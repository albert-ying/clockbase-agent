"""Find published clock analyses mentioning GEO series in the ClockBase atlas."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import re
import threading
import time
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

import pandas as pd
import requests

ATLAS = Path('/Entropy/kying/clockbase-verifier-v3-2026-09-27/final-2026-10-05/atlas-v3/atlas_v3_combined.csv')
OUTPUT = Path('/Entropy/kying/clockbase-benchmark-2026-10-05/retrieval')
REST = 'https://www.ebi.ac.uk/europepmc/webservices/rest'
ANNOTATIONS = 'https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds'
NCBI = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
TERMS = ('epigenetic clock', 'epigenetic age', 'DNAm age', 'DNA methylation age',
         'methylation age', 'epigenetic aging', 'age acceleration', 'PhenoAge',
         'GrimAge', 'DunedinPACE', 'Horvath clock', 'transcriptomic clock',
         'transcriptomic age', 'tAge', 'aging clock')
ACCESSION = re.compile(r'\b(GSE|GPL|GSM)(\d+)\b', re.I)
CLOCK = re.compile(r'epigenetic\s+(?:clock|age|aging|ageing)|(?:DNA\s*methylation|DNAm|methylation)\s*[- ]?age|'
                   r'age\s+acceleration|PhenoAge|GrimAge|Dunedin(?:PACE|PoAm)|Horvath|Hannum|'
                   r'skin[- ]and[- ]blood|transcriptomic\s+(?:clock|age)|\btAge\b|'
                   r'(?:aging|ageing|mortality|chronological)\s+clock|\bLin.{0,20}clock', re.I)
GROUP = re.compile(r'\b(?:controls?|versus|vs\.?|patients?|cases?|treated|treatment|untreated|'
                   r'intervention|disease|healthy|exposed|exposure|knockout|wild[- ]type|'
                   r'diet|restriction|obes(?:e|ity)|smok(?:er|ing)|HIV|infect(?:ed|ion)|'
                   r'cancer|tumou?r|diabet(?:es|ic)|syndrome)\b', re.I)
COMPARISON = re.compile(r'\b(?:compar(?:e[ds]?|ison|ing)|versus|vs\.?|differ(?:ence[sd]?|ent)|'
                        r'higher|lower|older|younger|unchanged|increas(?:e[ds]?|ing)|decreas(?:e[ds]?|ing)|'
                        r'accelerat(?:e[ds]?|ion)|decelerat(?:e[ds]?|ion)|reduc(?:e[ds]?|tion)|'
                        r'revers(?:e[ds]?|al))\b', re.I)
REPORT = re.compile(r'\b(?:p\s*(?:value)?\s*[=<>≤≥]|significant(?:ly)?|non[- ]significant|'
                    r'coefficient|effect\s+size|mean\s+difference|higher|lower|older|younger|unchanged|increased|'
                    r'decreased|accelerated|decelerated|reduced|reversed)\b', re.I)
APPLIED = re.compile(r'\b(?:calculat(?:e[ds]?|ing)|estimat(?:e[ds]?|ing)|predict(?:e[ds]?|ing|ion)|'
                     r'appl(?:y|ied)|assess(?:e[ds]?|ing)|measur(?:e[ds]?|ing)|evaluat(?:e[ds]?|ing)|'
                     r'analys(?:e[ds]?|is)|analyz(?:e[ds]?|ing)|observed|found|showed|revealed|'
                     r'accelerated|decelerated|higher|lower|older|younger|increased|decreased)\b', re.I)
EXCLUDE = re.compile(r'\b(?:training|train\s+set|trained|clock\s+development|chronological\s+age\s+prediction|'
                     r'age\s+correlation|correlation\s+with\s+(?:chronological\s+)?age)\b', re.I)


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    try:
        path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    except OSError as exc:
        raise RuntimeError(f'Cannot write {path}: {exc}') from exc


class CachedGET:
    """Cache response bytes and provenance; serialize starts to at most 2.5/s."""

    def __init__(self, output, offline=False):
        self.root = output / 'cache'
        self.root.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.lock = threading.Lock()
        self.last_start = 0.0
        self.network_requests = 0
        self.cache_hits = 0

    def get(self, url, params=None):
        if params:
            url += '?' + urlencode(sorted(params.items()), doseq=True)
        key = hashlib.sha256(url.encode()).hexdigest()
        meta_path = self.root / f'{key}.json'
        body_path = self.root / f'{key}.body'
        if meta_path.exists():
            try:
                meta = json.loads(meta_path.read_text())
                body = body_path.read_bytes()
                if hashlib.sha256(body).hexdigest() != meta['sha256']:
                    raise RuntimeError(f'Cache checksum mismatch: {body_path}')
                with self.lock:
                    self.cache_hits += 1
                return body, meta
            except (OSError, ValueError, KeyError) as exc:
                raise RuntimeError(f'Cannot read cached response {meta_path}: {exc}') from exc
        if self.offline:
            raise RuntimeError(f'Offline response missing: {url}')
        for attempt in range(4):
            with self.lock:
                time.sleep(max(0, .4 - (time.monotonic() - self.last_start)))
                self.last_start = time.monotonic()
                self.network_requests += 1
            try:
                response = requests.get(url, timeout=(15, 90), allow_redirects=False,
                                        headers={'User-Agent': 'ClockBase-literature-retrieval/2026-10-05',
                                                 'Accept': '*/*'})
                body = response.content
                meta = {'url': url, 'fetched_at': utcnow(), 'status': response.status_code,
                        'headers': dict(response.headers), 'sha256': hashlib.sha256(body).hexdigest(),
                        'bytes': len(body), 'cache_body': str(body_path), 'attempt': attempt + 1}
                if response.status_code in (429, 500, 502, 503, 504) and attempt < 3:
                    retry_path = self.root / f'{key}.attempt-{attempt + 1}'
                    retry_path.with_suffix(retry_path.suffix + '.body').write_bytes(body)
                    write_json(retry_path.with_suffix(retry_path.suffix + '.json'), meta)
                    delay = response.headers.get('Retry-After', '')
                    time.sleep(min(60, float(delay) if delay.isdigit() else 2 ** (attempt + 1)))
                    continue
                body_path.write_bytes(body)
                write_json(meta_path, meta)
                return body, meta
            except requests.RequestException as exc:
                logging.warning('GET attempt %s failed for %s: %s', attempt + 1, url, exc)
                write_json(self.root / f'{key}.error-{attempt + 1}.json',
                           {'url': url, 'attempt': attempt + 1, 'at': utcnow(), 'error': str(exc)})
                if attempt == 3:
                    raise RuntimeError(f'GET failed after four attempts: {url}') from exc
                time.sleep(2 ** (attempt + 1))
            except OSError as exc:
                raise RuntimeError(f'Cannot cache response for {url}: {exc}') from exc
        raise RuntimeError(f'GET exhausted retries: {url}')

    def json(self, url, params=None):
        body, meta = self.get(url, params)
        if meta['status'] != 200:
            raise RuntimeError(f"HTTP {meta['status']}: {meta['url']}")
        try:
            data = json.loads(body)
            if isinstance(data, dict) and (data.get('errCode') or data.get('error')):
                raise RuntimeError(f"API error at {meta['url']}: {data}")
            return data, meta
        except ValueError as exc:
            raise RuntimeError(f"Invalid JSON: {meta['url']}") from exc


def accession_ids(text):
    return sorted({prefix.upper() + number for prefix, number in ACCESSION.findall(text)})


def paper_key(record):
    return 'MED:' + record['pmid'] if record.get('pmid') else record['source'] + ':' + record['id']


def search(client, output):
    base = '(' + ' OR '.join('"' + term + '"' for term in TERMS) + ')'
    definitions = [('geo_annotation', base + ' AND ACCESSION_TYPE:geo'),
                   ('gse_literal', base + ' AND GSE'), ('gse_wildcard', base + ' AND GSE*')]
    papers, queries = {}, []
    for name, query in definitions:
        cursor, seen, pages = '*', set(), []
        hit_count = None
        while True:
            data, meta = client.json(REST + '/search', {'query': query, 'cursorMark': cursor,
                                                      'resultType': 'core', 'format': 'json', 'pageSize': 1000})
            if 'hitCount' not in data:
                raise RuntimeError(f'Unsuccessful search {query}: {data}')
            if hit_count is None:
                hit_count = data['hitCount']
            rows = data.get('resultList', {}).get('result', [])
            pages.append({'cursorMark': cursor, 'records': len(rows), 'response': meta})
            for record in rows:
                key = paper_key(record)
                seen.add(key)
                if key not in papers:
                    papers[key] = {'record': record, 'queries': []}
                papers[key]['queries'].append(name)
            new_cursor = data.get('nextCursorMark')
            if len(seen) >= hit_count or not rows:
                break
            if not new_cursor or new_cursor == cursor:
                raise RuntimeError(f'Paging stopped before hit count: {name} {len(seen)}/{hit_count}')
            cursor = new_cursor
        if len(seen) != hit_count:
            raise RuntimeError(f'Incomplete search: {name} {len(seen)}/{hit_count}')
        entry = {'name': name, 'query': query, 'hit_count': hit_count, 'retrieved_records': len(seen),
                 'date': pages[0]['response']['fetched_at'], 'pages': pages}
        queries.append(entry)
        logging.info('Search %s: %s hits, %s records', name, hit_count, len(seen))
    write_json(output / 'queries.json', {'terms': TERMS, 'queries': queries, 'unique_papers': len(papers),
                                       'fulltext_note': 'GSE* also matches accession tokens; bare GSE does not.'})
    write_json(output / 'search_records.json', papers)
    return papers, queries


def fetch_annotations(client, papers, output):
    annotations, errors = {}, []
    keys = sorted(papers)
    for offset in range(0, len(keys), 8):
        batch = keys[offset:offset + 8]
        ids = []
        lookup = {}
        for key in batch:
            record = papers[key]['record']
            article_id = 'PMC:' + record['pmcid'] if record.get('pmcid') else key
            ids.append(article_id)
            lookup[record.get('pmcid') or key] = key
            annotations[key] = {'annotations': [], 'status': 'not_returned'}
        try:
            data, meta = client.json(ANNOTATIONS, {'articleIds': ','.join(ids),
                                                  'type': 'Accession Numbers', 'format': 'JSON'})
            if not isinstance(data, list):
                raise RuntimeError(f'Unexpected annotation format: {data}')
            for result in data:
                key = lookup.get(result.get('pmcid')) or (result.get('source', '') + ':' + result.get('extId', ''))
                if key in annotations:
                    values = annotations[key]['annotations'] + [a for a in result.get('annotations', [])
                                                               if a.get('subType', '').upper() == 'GEO' or
                                                               any('geo:' in t.get('uri', '').lower()
                                                                   for t in a.get('tags', []))]
                    unique = {json.dumps(a, sort_keys=True): a for a in values}
                    annotations[key] = {'annotations': list(unique.values()), 'status': 'ok', 'response': meta}
        except RuntimeError as exc:
            errors.append({'stage': 'annotations', 'papers': batch, 'error': str(exc)})
            for key in batch:
                annotations[key] = {'annotations': [], 'status': 'failed', 'error': str(exc)}
            logging.warning('Annotation batch failed: %s', exc)
        if offset % 160 == 0:
            logging.info('Annotations: %s/%s', min(offset + 8, len(keys)), len(keys))
    write_json(output / 'annotations.json', annotations)
    return annotations, errors


def ncbi_accessions(record):
    values = set(accession_ids(str(record.get('accession', ''))))
    # GDS summaries sometimes identify the originating series in a numeric gse field.
    for number in re.findall(r'\b\d+\b', str(record.get('gse', ''))):
        if int(number) > 0:
            values.add('GSE' + number)
    return sorted(values)


def fetch_ncbi_links(client, papers, output):
    """Map each PMID separately to curated GEO links, then resolve GEO summaries."""
    links, errors, summaries = {}, [], {}
    keys = sorted(papers)
    pubmed = []
    for key in keys:
        record = papers[key]['record']
        if record.get('pmid'):
            pubmed.append(key)
            links[key] = {'status': 'not_returned', 'geo_uids': [], 'accession_evidence': []}
        else:
            links[key] = {'status': 'not_pubmed', 'geo_uids': [], 'accession_evidence': []}
    for offset in range(0, len(pubmed), 30):
        batch = pubmed[offset:offset + 30]
        lookup = {papers[k]['record']['pmid']: k for k in batch}
        try:
            data, meta = client.json(NCBI + '/elink.fcgi', {'dbfrom': 'pubmed', 'db': 'gds',
                                                         'id': list(lookup), 'retmode': 'json'})
            if 'linksets' not in data:
                raise RuntimeError(f'Unexpected NCBI link response: {data}')
            for result in data['linksets']:
                ids = [str(v) for v in result.get('ids', [])]
                if len(ids) != 1 or ids[0] not in lookup:
                    raise RuntimeError(f'Ambiguous PMID/GEO linkage: {ids}')
                key = lookup[ids[0]]
                uids = sorted({str(uid) for db in result.get('linksetdbs', []) if db.get('dbto') == 'gds'
                               for uid in db.get('links', [])})
                links[key] = {'status': 'ok', 'geo_uids': uids, 'accession_evidence': [], 'link_response': meta}
        except RuntimeError as exc:
            errors.append({'stage': 'ncbi_links', 'papers': batch, 'error': str(exc)})
            for key in batch:
                links[key]['status'] = 'failed'
            logging.warning('NCBI link batch failed: %s', exc)
        if offset % 300 == 0:
            logging.info('NCBI PMID/GEO cross-references: %s/%s', min(offset + 30, len(pubmed)), len(pubmed))
    uids = sorted({uid for entry in links.values() for uid in entry['geo_uids']})
    for offset in range(0, len(uids), 200):
        batch = uids[offset:offset + 200]
        try:
            data, meta = client.json(NCBI + '/esummary.fcgi', {'db': 'gds', 'id': ','.join(batch), 'retmode': 'json'})
            result = data.get('result', {})
            for uid in batch:
                if uid not in result or result[uid].get('error'):
                    raise RuntimeError(f'GEO summary missing/failed: {uid}')
                summaries[uid] = {'record': result[uid], 'response': meta}
        except RuntimeError as exc:
            errors.append({'stage': 'ncbi_geo_summary', 'geo_uids': batch, 'error': str(exc)})
            logging.warning('NCBI GEO summary batch failed: %s', exc)
    for key, entry in links.items():
        for uid in entry['geo_uids']:
            if uid not in summaries:
                entry['status'] = 'geo_summary_failed'
                continue
            item = summaries[uid]
            for accession in ncbi_accessions(item['record']):
                entry['accession_evidence'].append({'accession': accession, 'source': 'ncbi_pubmed_geo',
                                                   'geo_uid': uid, 'response_sha256': item['response']['sha256'],
                                                   'url': item['response']['url'],
                                                   'pubmed_link_response_sha256': entry['link_response']['sha256']})
    write_json(output / 'ncbi_crossreferences.json', {'links': links, 'geo_summaries': summaries})
    return links, errors


def normalized_text(element):
    return re.sub(r'\s+', ' ', ''.join(element.itertext())).strip()


def text_blocks(xml):
    """Keep article prose, table/figure text and section labels, excluding references."""
    root = ET.fromstring(xml)
    blocks = []

    def prose_parts(element):
        if element.text:
            yield element.text
        for child in element:
            child_tag = child.tag.rsplit('}', 1)[-1]
            if child_tag == 'xref':
                yield ' ' + normalized_text(child) + ' '
            elif child_tag not in ('fig', 'table-wrap', 'supplementary-material', 'table', 'ref-list'):
                yield from prose_parts(child)
            if child.tail:
                yield child.tail

    def embedded(element, section, label):
        for child in element:
            if child.tag.rsplit('}', 1)[-1] in ('fig', 'table-wrap', 'supplementary-material'):
                walk(child, section, label)
            else:
                embedded(child, section, label)

    def walk(element, section, label):
        tag = element.tag.rsplit('}', 1)[-1]
        if tag == 'article-meta':
            for abstract in element.findall('abstract'):
                walk(abstract, ['Abstract'], None)
            return
        if tag in ('ref-list', 'ref', 'journal-meta'):
            return
        if tag == 'sec':
            title = element.find('title')
            heading = normalized_text(title) if title is not None else element.get('sec-type', 'Section')
            section = section + [heading]
        elif tag in ('fig', 'table-wrap', 'supplementary-material'):
            child = element.find('label')
            label = normalized_text(child) if child is not None else element.get('id', tag)
        if tag in ('p', 'table', 'caption', 'abstract', 'title'):
            # Descend into structured abstracts/captions; emit their paragraph leaves once.
            if tag in ('abstract', 'caption') and any(c.tag in ('p', 'sec') for c in element):
                pass
            elif tag == 'title':
                return
            else:
                if tag == 'p':
                    text = re.sub(r'\s+', ' ', ''.join(prose_parts(element))).strip()
                else:
                    text = re.sub(r'\s+', ' ', ' '.join(element.itertext())).strip() if tag == 'table' else normalized_text(element)
                if text:
                    blocks.append({'text': text, 'section': ' / '.join(section) or 'Article',
                                   'label': label, 'element': tag, 'id': element.get('id'),
                                   'block_index': len(blocks)})
                if tag == 'p':
                    embedded(element, section, label)
                return
        for child in element:
            child_section = section
            child_tag = child.tag.rsplit('}', 1)[-1]
            if child_tag == 'abstract':
                child_section = ['Abstract']
            elif child_tag == 'body':
                child_section = ['Body']
            elif child_tag == 'back':
                child_section = ['Back matter']
            walk(child, child_section, label)

    walk(root, [], None)
    return blocks


def abstract_blocks(record):
    value = record.get('abstractText', '')
    try:
        text = normalized_text(ET.fromstring('<abstract>' + value + '</abstract>'))
    except ET.ParseError:
        text = re.sub(r'<[^>]*>', ' ', value)
        text = re.sub(r'\s+', ' ', text).strip()
    return [{'text': text, 'section': 'Abstract', 'label': None, 'element': 'abstract',
             'id': None, 'block_index': 0}] if text else []


def fetch_paper(client, key, item, annotation, ncbi_link, output):
    record = item['record']
    evidence, errors, responses = [], [], []
    blocks = abstract_blocks(record)
    fulltext_status = 'not_open_access' if record.get('isOpenAccess') != 'Y' else 'no_pmcid'
    if record.get('isOpenAccess') == 'Y' and record.get('pmcid'):
        try:
            body, meta = client.get(REST + '/' + record['pmcid'] + '/fullTextXML')
            responses.append({'kind': 'fulltext', **meta})
            if meta['status'] == 200:
                blocks = text_blocks(body)
                fulltext_status = 'ok'
                # Full-XML accession census includes references; triage uses non-reference prose.
                for accession in accession_ids(body.decode('utf-8', errors='replace')):
                    evidence.append({'accession': accession, 'source': 'fulltext_xml',
                                     'response_sha256': meta['sha256'], 'url': meta['url']})
            else:
                fulltext_status = 'http_' + str(meta['status'])
                errors.append({'stage': 'fulltext', 'paper': key, 'error': fulltext_status})
        except (RuntimeError, ET.ParseError) as exc:
            fulltext_status = 'failed'
            errors.append({'stage': 'fulltext', 'paper': key, 'error': str(exc)})
    for accession in accession_ids(json.dumps(record)):
        evidence.append({'accession': accession, 'source': 'core_record', 'paper': key})
    for entry in annotation.get('annotations', []):
        for accession in accession_ids(entry.get('exact', '') + ' ' + json.dumps(entry.get('tags', []))):
            evidence.append({'accession': accession, 'source': 'geo_annotation',
                             'section': entry.get('section'), 'annotation_id': entry.get('id'),
                             'context': entry.get('prefix', '') + entry.get('exact', '') + entry.get('postfix', ''),
                             'response_sha256': annotation['response']['sha256']})
    evidence.extend(ncbi_link['accession_evidence'])
    if ncbi_link.get('link_response'):
        responses.append({'kind': 'ncbi_pubmed_geo_links', **ncbi_link['link_response']})
    accessions = sorted({e['accession'] for e in evidence})
    paper = {'paper_id': key, 'ids': {k: record.get(k) for k in ('id', 'source', 'pmid', 'pmcid', 'doi')},
             'title': record.get('title'), 'journal': record.get('journalInfo', {}).get('journal', {}).get('title'),
             'year': record.get('pubYear'), 'first_publication_date': record.get('firstPublicationDate'),
             'oa': record.get('isOpenAccess') == 'Y', 'publication_types': record.get('pubTypeList', {}).get('pubType', []),
             'queries': sorted(set(item['queries'])), 'accessions': [a for a in accessions if a.startswith('GSE')],
             'gpl_mentions': [a for a in accessions if a.startswith('GPL')],
             'gsm_mentions': [a for a in accessions if a.startswith('GSM')], 'accession_evidence': evidence,
             'fulltext_status': fulltext_status, 'database_links_status': ncbi_link['status'],
             'europepmc_db_names': record.get('dbCrossReferenceList', {}).get('dbName', []),
             'annotations_status': annotation['status'], 'responses': responses, 'errors': errors}
    path = output / 'texts' / (hashlib.sha256(key.encode()).hexdigest() + '.json')
    write_json(path, {'paper_id': key, 'blocks': blocks, 'fulltext_status': fulltext_status})
    paper['text_path'] = str(path)
    return paper, blocks


def load_atlas(path):
    if not path.is_file():
        raise RuntimeError(f'Cannot proceed - required data file not found: {path}')
    frame = pd.read_csv(path, usecols=['modality', 'gse_id', 'clock'])
    atlas = {}
    for gse, group in frame.groupby('gse_id'):
        atlas[str(gse).strip().upper()] = {
            'modalities': sorted(group.modality.dropna().unique().tolist()), 'atlas_rows': len(group),
            'atlas_clocks': sorted(group.clock.dropna().unique().tolist()), 'atlas_clock_count': group.clock.nunique(),
            'atlas_by_modality': {modality: {'rows': len(sub), 'clocks': sorted(sub.clock.dropna().unique().tolist())}
                                  for modality, sub in group.groupby('modality')}}
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    provenance = {'path': str(path), 'sha256': digest.hexdigest(), 'bytes': path.stat().st_size,
                  'modified_at': datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
                  'rows': len(frame), 'series': len(atlas), 'rows_by_modality': frame.modality.value_counts().to_dict(),
                  'preprocessing': 'Read modality/gse_id/clock; uppercase and strip GSE keys; group without filtering rows.'}
    return atlas, provenance


def sentences(text):
    pattern = r'(?<!vs\.)(?<!Fig\.)(?<!fig\.)(?<!al\.)(?<!e\.g\.)(?<!i\.e\.)(?<=[.!?])\s+(?=[A-Z(])'
    return [s.strip() for s in re.split(pattern, text) if s.strip()]


def passage(block, sentence, role):
    return {'text': sentence, 'section': block['section'], 'label': block['label'],
            'element': block['element'], 'block_index': block['block_index'], 'role': role}


def triage(paper, gse, blocks, atlas_entry):
    """Require accession-local clock/comparison evidence for automatic LIKELY."""
    pattern = re.compile(r'\b' + re.escape(gse) + r'\b', re.I)
    anchors = [i for i, b in enumerate(blocks) if pattern.search(b['text'])]
    clock_blocks = [i for i, b in enumerate(blocks) if CLOCK.search(b['text'])]
    direct, near, remote, anchor_passages = [], [], [], []
    training = False
    for i in anchors:
        block = blocks[i]
        training = training or bool(EXCLUDE.search(block['text']))
        for sentence in sentences(block['text']):
            if pattern.search(sentence):
                anchor_passages.append(passage(block, sentence, 'series_mention'))
        # Require a clock result sentence within one sentence of the accession mention.
        parts = sentences(block['text'])
        if block['element'] != 'table' and not EXCLUDE.search(block['text']) and \
                not re.search(r'availability|references|methods?|introduction|supplement|associated data', block['section'], re.I):
            for index, sentence in enumerate(parts):
                if not (CLOCK.search(sentence) and COMPARISON.search(sentence) and REPORT.search(sentence)):
                    continue
                window = ' '.join(parts[max(0, index - 1):index + 2])
                if pattern.search(window) and GROUP.search(window) and APPLIED.search(window):
                    direct.append(i)
                    break
        for j in range(max(0, i - 1), min(len(blocks), i + 2)):
            if blocks[j]['section'] == block['section'] and CLOCK.search(blocks[j]['text']):
                near.append(j)
    for i in clock_blocks:
        if GROUP.search(blocks[i]['text']) and COMPARISON.search(blocks[i]['text']) and REPORT.search(blocks[i]['text']):
            remote.append(i)
    review = any('review' in p.lower() for p in paper['publication_types'])
    if direct and not review:
        classification = 'LIKELY'
        reason = 'A clock-result sentence and this GSE occur within one sentence in the same prose paragraph, with group/application language; manual confirmation needed.'
    elif review:
        classification = 'UNLIKELY'
        reason = 'Indexed as a review; accession/comparison mentions do not establish a primary clock analysis.'
    elif near or (anchors and remote) or (not anchors and clock_blocks):
        classification = 'POSSIBLE'
        reason = 'Clock analysis and this accession are present, but text rules cannot bind the reported comparison to this series.'
        if training:
            reason += ' Series also occurs in training/calibration context.'
        if paper['fulltext_status'] != 'ok':
            reason += ' Full text unavailable; abstract/annotation evidence only.'
    else:
        classification = 'UNLIKELY'
        reason = 'No accession-local clock comparison found in available prose; series may serve another analysis or a citation.'
    selected = []
    for i in list(dict.fromkeys(direct + near + remote))[:8]:
        for sentence in sentences(blocks[i]['text']):
            if CLOCK.search(sentence) or (i in direct and (GROUP.search(sentence) or REPORT.search(sentence))):
                selected.append(passage(blocks[i], sentence, 'clock_comparison' if i in direct else
                                        'nearby_clock_context' if i in near else 'unbound_clock_comparison'))
    for entry in paper['accession_evidence']:
        if entry['accession'] == gse and entry['source'] == 'geo_annotation' and not anchor_passages:
            anchor_passages.append({'text': entry['context'], 'section': entry.get('section'),
                                    'label': None, 'element': 'annotation_context', 'block_index': None,
                                    'role': 'series_mention'})
    unique = {}
    for p in anchor_passages[:10] + selected[:25]:
        unique[(p['block_index'], p['text'])] = p
    return {'paper_id': paper['paper_id'], 'ids': paper['ids'], 'title': paper['title'], 'year': paper['year'],
            'gse_id': gse, 'modality': '; '.join(atlas_entry['modalities']), **atlas_entry,
            'triage_class': classification, 'reason': reason, 'triage_method': 'deterministic_text_rules',
            'fulltext_status': paper['fulltext_status'], 'passages': list(unique.values()),
            'text_path': paper['text_path'], 'accession_evidence': [e for e in paper['accession_evidence'] if e['accession'] == gse],
            'ranking_score': 100 * len(direct) + 10 * len(set(near)) + len(remote),
            'rule_features': {'series_blocks': anchors, 'direct_blocks': direct, 'near_clock_blocks': sorted(set(near)),
                              'unbound_comparison_blocks': remote, 'training_context': training, 'review': review}}


def apply_manual_reviews(candidates, output):
    path = output / 'manual_reviews.jsonl'
    if not path.exists():
        return
    pairs = {(c['paper_id'], c['gse_id']): c for c in candidates}
    with path.open() as handle:
        for line in handle:
            review = json.loads(line)
            key = (review['paper_id'], review['gse_id'])
            if key not in pairs:
                raise RuntimeError(f'Manual review has no fetched atlas candidate: {key}')
            candidate = pairs[key]
            blocks = json.loads(Path(candidate['text_path']).read_text())['blocks']
            for p in review['passages']:
                index = p['block_index']
                if index is None or not 0 <= index < len(blocks) or \
                        p['text'] not in blocks[index]['text'] or p['section'] != blocks[index]['section'] or \
                        p.get('label') != blocks[index]['label']:
                    raise RuntimeError(f'Manual passage absent from cached text: {key}: {p["text"][:80]}')
            if review['triage_class'] not in ('LIKELY', 'POSSIBLE', 'UNLIKELY'):
                raise RuntimeError(f'Invalid manual triage class: {key}')
            candidate['automatic_triage_class'] = candidate['triage_class']
            candidate['automatic_reason'] = candidate['reason']
            candidate['triage_class'] = review['triage_class']
            candidate['reason'] = review['reason']
            candidate['triage_method'] = 'manual_reading_cached_text'
            candidate['manual_review'] = review
            candidate['passages'] = review['passages']
            candidate['ranking_score'] = review.get('ranking_score', candidate['ranking_score'])


def report(summary, top, output):
    counts = summary['triage_counts']
    lines = ['# Published aging-clock retrieval, 2026-10-05', '',
             f"Retrieved {summary['unique_papers']} distinct indexed papers; {summary['papers_with_gse']} mention a GSE.",
             f"The atlas intersection contains {summary['paper_gse_pairs']} paper/GSE pairs, "
             f"{summary['papers_in_atlas']} papers and {summary['series_in_atlas']} series.", '',
             '| Modality | LIKELY | POSSIBLE | UNLIKELY |', '|---|---:|---:|---:|']
    for modality, counter in sorted(summary['triage_by_modality'].items()):
        lines.append(f"| {modality} | {counter.get('LIKELY', 0)} | {counter.get('POSSIBLE', 0)} | {counter.get('UNLIKELY', 0)} |")
    lines += ['', f"Overall: {counts.get('LIKELY', 0)} LIKELY, {counts.get('POSSIBLE', 0)} POSSIBLE, "
              f"{counts.get('UNLIKELY', 0)} UNLIKELY. A pair in several modalities contributes to each modality row.", '',
              '## Method', '',
              'Europe PMC REST search used cursorMark paging, resultType=core, pageSize=1000, and the clock terms '
              'listed in queries.json. Three queries cover GEO annotation, literal GSE, and GSE*; the wildcard '
              'matches accession tokens that a bare GSE term misses. No publication-date or species filter was applied. '
              'The hit counts describe the index at each cached fetch date, not all papers that exist.', '',
              'Series were recovered from core records, Europe PMC GEO accession annotations (batches of eight), '
              'NCBI curated PubMed-to-GEO database cross-references (ELink with separate repeated PMID parameters '
              'and GEO ESummary), and a GSE/GPL/GSM regex over available '
              'open-access fullTextXML. GPL and GSM mentions are recorded separately. Reference-list accessions '
              'enter the accession census, but references are excluded from triage prose. Records are deduplicated '
              'by PMID when present, otherwise Europe PMC source/id. Publication versions without a shared PMID '
              'remain separate records.', '',
              'Atlas linkage is an exact GSE intersection. Rows and distinct clock names are counted from the '
              'provided CSV without filtering verifier classes. An accession match does not establish matching '
              'samples, contrasts, species, clocks or an eligible verified atlas row. Those checks belong to the '
              'statistic-extraction and matching steps.', '',
              'Deterministic LIKELY requires a sentence containing a clock, comparison language and a reported '
              'statistic or direction, with this accession, application and group language within one adjacent '
              'sentence in the same non-reference prose paragraph. '
              'Training/calibration context, methods, introductions, supplementary descriptions, '
              'data-availability sections, review articles and whole-table text '
              'cannot yield automatic LIKELY. Nearby or unbound clock comparisons yield POSSIBLE; other matches '
              'yield UNLIKELY. Manual review records in manual_reviews.jsonl override the rules with a reason and '
              'passages checked verbatim against cached text. The automatic class and reason remain on reviewed '
              'candidates. Supporting passages retain section, figure/table label and block index; unbound passages '
              'are explicitly labeled. These are retrieval classifications, not accepted benchmark ground truth.', '',
              '## Limits and provenance', '',
              'Paywalled full texts were not read. Full-text failures use the abstract and accession annotations. '
              'Figures are represented by captions and supplements by embedded descriptions; supplementary files '
              'and figure pixels were not downloaded or interpreted. A result only in a figure/supplement can remain '
              'POSSIBLE. Abstracts often omit GSE identifiers, and papers with neither searchable clock terms nor '
              'searchable accessions can be missed. The clock vocabulary is finite; indirect cohort names require '
              'manual review. Direction language is permitted as a reported comparison, but age-prediction '
              'performance alone is not a case-control result. Statistical values have not been extracted.', '',
              'Every HTTP response, including errors and transient retry responses, is preserved in cache/ with '
              'request URL, timestamp, status and SHA256. Request starts share a 0.4-second interval (at most 2.5/s), '
              'with four attempts and backoff for transport errors/429/5xx. Cached HTTP failures are reused; '
              'reruns do not automatically retry them. Offline mode fails on a missing cache entry. Search failures '
              'stop the run; per-paper source failures are counted in summary.json. texts/ holds normalized source '
              'blocks. search_records.json preserves core records; annotations.json preserves GEO annotations. '
              'ncbi_crossreferences.json preserves PMID-to-GEO UID mappings and summary records. Europe PMC '
              'databaseLinks does not support GEO as a database filter (cached API-error probes); core database '
              'names are retained, and NCBI supplies the curated GEO route.', '',
              f"Full-text status: `{json.dumps(summary['fulltext_status'], sort_keys=True)}`. "
              f"Source failures recorded: {len(summary['errors'])}.", '',
              f"Atlas: `{summary['atlas']['path']}`; SHA256 `{summary['atlas']['sha256']}`; "
              f"source modification timestamp `{summary['atlas']['modified_at']}`.", '',
              '## Rerun from the repository root', '', '```sh',
              'TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m benchmark.published_clock_2026_10_05.retrieve --offline',
              'TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m unittest benchmark.published_clock_2026_10_05.test_retrieve',
              'TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m benchmark.published_clock_2026_10_05.validate',
              '```', '', 'Omit --offline to fetch uncached URLs. --atlas and --output accept explicit alternative paths.', '',
              '## Most promising LIKELY pairs', '',
              'Ranked by manual-review priority when supplied, otherwise the number of accession-local and nearby '
              'clock/comparison passages. Several series from one paper are separate pairs. The complete supporting '
              'passages are in top_30_likely.jsonl and candidates.jsonl.', '',
              '| Rank | Paper | GSE | Atlas modality | Title |', '|---:|---|---|---|---|']
    for number, candidate in enumerate(top, 1):
        ids = candidate['ids']
        url = ('https://europepmc.org/articles/' + ids['pmcid']) if ids.get('pmcid') else \
              'https://europepmc.org/article/' + ids['source'] + '/' + ids['id']
        title = (candidate['title'] or 'Title unavailable in the core record').replace('|', '\\|')
        lines.append(f"| {number} | [{candidate['paper_id']}]({url}) | {candidate['gse_id']} | {candidate['modality']} | {title} |")
    lines += ['']
    try:
        (output / 'README.md').write_text('\n'.join(lines))
    except OSError as exc:
        raise RuntimeError(f'Cannot write retrieval README: {exc}') from exc


def write_jsonl(path, rows):
    try:
        with path.open('w') as handle:
            for row in rows:
                handle.write(json.dumps(row, ensure_ascii=False) + '\n')
    except OSError as exc:
        raise RuntimeError(f'Cannot write {path}: {exc}') from exc


def run(args):
    started = utcnow()
    atlas, provenance = load_atlas(args.atlas)
    client = CachedGET(args.output, args.offline)
    search_papers, queries = search(client, args.output)
    annotations, errors = fetch_annotations(client, search_papers, args.output)
    ncbi_links, ncbi_errors = fetch_ncbi_links(client, search_papers, args.output)
    errors.extend(ncbi_errors)
    (args.output / 'texts').mkdir(exist_ok=True)
    papers, candidates = [], []

    def process(key):
        return fetch_paper(client, key, search_papers[key], annotations[key], ncbi_links[key], args.output)

    with ThreadPoolExecutor(max_workers=3) as pool:
        for number, (paper, blocks) in enumerate(pool.map(process, sorted(search_papers)), 1):
            papers.append(paper)
            errors.extend(paper['errors'])
            for gse in paper['accessions']:
                if gse in atlas:
                    candidates.append(triage(paper, gse, blocks, atlas[gse]))
            if number % 50 == 0:
                logging.info('Papers %s/%s; atlas pairs %s', number, len(search_papers), len(candidates))
    apply_manual_reviews(candidates, args.output)
    candidates.sort(key=lambda c: (c['paper_id'], c['gse_id']))
    write_jsonl(args.output / 'papers.jsonl', papers)
    write_jsonl(args.output / 'candidates.jsonl', candidates)
    counts = Counter(c['triage_class'] for c in candidates)
    by_modality = defaultdict(Counter)
    for candidate in candidates:
        for modality in candidate['modalities']:
            by_modality[modality][candidate['triage_class']] += 1
    top = sorted((c for c in candidates if c['triage_class'] == 'LIKELY'),
                 key=lambda c: (-c['ranking_score'], c['paper_id'], c['gse_id']))[:30]
    write_jsonl(args.output / 'top_30_likely.jsonl', top)
    summary = {'started_at': started, 'completed_at': utcnow(), 'atlas': provenance,
               'source_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'query_hits': {q['name']: q['hit_count'] for q in queries}, 'unique_papers': len(papers),
               'papers_with_gse': sum(bool(p['accessions']) for p in papers),
               'distinct_gse_retrieved': len({a for p in papers for a in p['accessions']}),
               'fulltext_status': dict(Counter(p['fulltext_status'] for p in papers)),
               'annotation_status': dict(Counter(p['annotations_status'] for p in papers)),
               'database_link_status': dict(Counter(p['database_links_status'] for p in papers)),
               'papers_in_atlas': len({c['paper_id'] for c in candidates}),
               'series_in_atlas': len({c['gse_id'] for c in candidates}), 'paper_gse_pairs': len(candidates),
               'triage_counts': dict(counts), 'triage_by_modality': dict(by_modality),
               'triage_methods': dict(Counter(c['triage_method'] for c in candidates)),
               'accession_sources': dict(Counter(e['source'] for p in papers for e in p['accession_evidence'])),
               'top_likely_pairs': [{'paper_id': c['paper_id'], 'gse_id': c['gse_id'], 'title': c['title'],
                                     'modality': c['modality']} for c in top],
               'errors': errors, 'network_requests_this_run': client.network_requests,
               'cache_hits_this_run': client.cache_hits,
               'method': 'Literature retrieval and deterministic triage only; no published statistics extracted.'}
    write_json(args.output / 'summary.json', summary)
    report(summary, top, args.output)
    logging.info('Completed: %s papers, %s pairs, %s', len(papers), len(candidates), dict(counts))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--atlas', type=Path, default=ATLAS)
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    try:
        run(args)
    except (RuntimeError, OSError, ValueError) as exc:
        logging.error('%s', exc)
        raise SystemExit(1) from exc


if __name__ == '__main__':
    main()
