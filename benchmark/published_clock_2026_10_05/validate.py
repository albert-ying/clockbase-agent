"""Read back retrieval outputs against cached sources and the input atlas."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import logging
from pathlib import Path

import pandas as pd

from benchmark.published_clock_2026_10_05.retrieve import (
    ATLAS, OUTPUT, abstract_blocks, accession_ids, ncbi_accessions, paper_key, text_blocks, write_json,
)


def read_jsonl(path):
    with path.open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def validate(output, atlas_path):
    summary = json.loads((output / 'summary.json').read_text())
    queries = json.loads((output / 'queries.json').read_text())['queries']
    records = json.loads((output / 'search_records.json').read_text())
    annotations = json.loads((output / 'annotations.json').read_text())
    ncbi = json.loads((output / 'ncbi_crossreferences.json').read_text())
    papers = read_jsonl(output / 'papers.jsonl')
    candidates = read_jsonl(output / 'candidates.jsonl')
    top = read_jsonl(output / 'top_30_likely.jsonl')
    checks = Counter()

    def require(condition, message, category):
        if not condition:
            raise RuntimeError(message)
        checks[category] += 1

    metadata = {}
    for path in sorted((output / 'cache').glob('*.json')):
        meta = json.loads(path.read_text())
        if 'sha256' not in meta:
            continue
        body = path.with_suffix('.body').read_bytes()
        require(hashlib.sha256(body).hexdigest() == meta['sha256'], f'Cache checksum mismatch: {path}', 'cache_checksums')
        require(len(body) == meta['bytes'], f'Cache length mismatch: {path}', 'cache_lengths')
        metadata[meta['sha256']] = meta
    parsed_responses = {}

    def cached_json(digest):
        if digest not in parsed_responses:
            parsed_responses[digest] = json.loads(Path(metadata[digest]['cache_body']).read_text())
        return parsed_responses[digest]

    fetched = {}
    for query in queries:
        keys = set()
        for page in query['pages']:
            data = json.loads(Path(page['response']['cache_body']).read_text())
            require(data['request']['queryString'] == query['query'], f'Query mismatch: {query["name"]}', 'query_provenance')
            rows = data['resultList']['result']
            for record in rows:
                key = paper_key(record)
                keys.add(key)
                fetched.setdefault(key, record)
        require(len(keys) == query['retrieved_records'] == query['hit_count'],
                f'Incomplete search: {query["name"]}', 'paging_completeness')
    require(set(fetched) == set(records), 'Search/core-record census differs', 'paper_census')
    paper_map = {p['paper_id']: p for p in papers}
    require(len(paper_map) == len(papers) == len(records) == summary['unique_papers'], 'Paper census mismatch', 'paper_census')
    source_accessions, block_cache = {}, {}
    for paper in papers:
        key = paper['paper_id']
        record = records[key]['record']
        require(record == fetched[key], f'Core record differs from cached response: {key}', 'core_provenance')
        require(paper['title'] == record.get('title'), f'Title not from core: {key}', 'core_provenance')
        annotation = annotations[key]
        if annotation['status'] == 'ok':
            data = cached_json(annotation['response']['sha256'])
            articles = [item for item in data if
                        (record.get('pmcid') and item.get('pmcid') == record['pmcid']) or
                        item.get('source', '') + ':' + item.get('extId', '') == key]
            require(bool(articles), f'Annotation article absent from fetched response: {key}', 'annotation_provenance')
            expected = [a for article in articles for a in article.get('annotations', [])
                        if a.get('subType', '').upper() == 'GEO' or
                        any('geo:' in t.get('uri', '').lower() for t in a.get('tags', []))]
            expected = list({json.dumps(a, sort_keys=True): a for a in expected}.values())
            require(annotation['annotations'] == expected, f'GEO annotations differ from fetched response: {key}', 'annotation_provenance')
        blocks = json.loads(Path(paper['text_path']).read_text())['blocks']
        if paper['fulltext_status'] == 'ok':
            response = next(r for r in paper['responses'] if r['kind'] == 'fulltext')
            expected_blocks = text_blocks(Path(response['cache_body']).read_bytes())
        else:
            expected_blocks = abstract_blocks(record)
        require(blocks == expected_blocks, f'Normalized text differs from fetched source: {key}', 'text_provenance')
        block_cache[key] = blocks
        for evidence in paper['accession_evidence']:
            accession = evidence['accession']
            source = evidence['source']
            if source == 'core_record':
                values = accession_ids(json.dumps(record))
            elif source == 'geo_annotation':
                values = {a for entry in annotations[key]['annotations']
                          for a in accession_ids(entry.get('exact', '') + ' ' + json.dumps(entry.get('tags', [])))}
                require(evidence['response_sha256'] in metadata, f'Annotation response missing: {key}', 'accession_source')
            elif source == 'ncbi_pubmed_geo':
                digest = evidence['response_sha256']
                require(digest in metadata and evidence['pubmed_link_response_sha256'] in metadata,
                        f'NCBI response missing: {key}', 'accession_source')
                link_data = cached_json(evidence['pubmed_link_response_sha256'])
                uid = evidence['geo_uid']
                require(any([str(v) for v in s.get('ids', [])] == [record['pmid']] and
                            any(uid in [str(v) for v in db.get('links', [])] for db in s.get('linksetdbs', [])
                                if db.get('dbto') == 'gds') for s in link_data['linksets']),
                        f'GEO UID is not linked to this PMID: {key} {uid}', 'ncbi_link_provenance')
                summary_data = cached_json(digest)['result'][uid]
                require(summary_data == ncbi['geo_summaries'][uid]['record'], f'GEO summary differs: {uid}', 'ncbi_link_provenance')
                values = ncbi_accessions(summary_data)
            else:
                digest = evidence['response_sha256']
                require(digest in metadata, f'Accession response missing: {key} {accession}', 'accession_source')
                if digest not in source_accessions:
                    body = Path(metadata[digest]['cache_body']).read_text(errors='replace')
                    source_accessions[digest] = set(accession_ids(body))
                values = source_accessions[digest]
            require(accession in values, f'Accession not in fetched source: {key} {accession}', 'accession_source')
        require(set(paper['accessions']) == {e['accession'] for e in paper['accession_evidence']
                                            if e['accession'].startswith('GSE')},
                f'GSE accession census differs: {key}', 'accession_census')
        expected_evidence = {('core_record', a) for a in accession_ids(json.dumps(record))}
        expected_evidence.update(('geo_annotation', a) for entry in annotation['annotations']
                                 for a in accession_ids(entry.get('exact', '') + ' ' + json.dumps(entry.get('tags', []))))
        expected_evidence.update(('ncbi_pubmed_geo', a) for uid in ncbi['links'][key]['geo_uids']
                                 if uid in ncbi['geo_summaries']
                                 for a in ncbi_accessions(ncbi['geo_summaries'][uid]['record']))
        if paper['fulltext_status'] == 'ok':
            expected_evidence.update(('fulltext_xml', a) for a in
                                     accession_ids(Path(response['cache_body']).read_text(errors='replace')))
        require(expected_evidence == {(e['source'], e['accession']) for e in paper['accession_evidence']},
                f'Accession extraction incomplete: {key}', 'accession_census')
    if not atlas_path.is_file():
        raise RuntimeError(f'Cannot proceed - required data file not found: {atlas_path}')
    digest = hashlib.sha256()
    with atlas_path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    require(digest.hexdigest() == summary['atlas']['sha256'], 'Atlas changed since retrieval', 'atlas_provenance')
    frame = pd.read_csv(atlas_path, usecols=['gse_id', 'modality', 'clock'])
    frame['gse_id'] = frame.gse_id.str.strip().str.upper()
    row_counts = frame.groupby('gse_id').size().to_dict()
    clocks = frame.groupby('gse_id').clock.agg(lambda x: sorted(set(x.dropna()))).to_dict()
    modalities = frame.groupby('gse_id').modality.agg(lambda x: sorted(set(x.dropna()))).to_dict()
    expected_pairs = {(p['paper_id'], gse) for p in papers for gse in p['accessions'] if gse in row_counts}
    candidate_map = {(c['paper_id'], c['gse_id']): c for c in candidates}
    require(set(candidate_map) == expected_pairs and len(candidate_map) == len(candidates),
            'Atlas intersection incomplete or duplicated', 'intersection')
    require(len(candidates) == summary['paper_gse_pairs'], 'Candidate total differs', 'counts')
    counts = Counter()
    by_modality = defaultdict(Counter)
    for candidate in candidates:
        key, gse = candidate['paper_id'], candidate['gse_id']
        require(candidate['atlas_rows'] == row_counts[gse] and candidate['atlas_clocks'] == clocks[gse] and
                candidate['atlas_clock_count'] == len(clocks[gse]) and candidate['modalities'] == modalities[gse],
                f'Atlas counts differ: {key} {gse}', 'atlas_counts')
        require(candidate['accession_evidence'] == [e for e in paper_map[key]['accession_evidence'] if e['accession'] == gse],
                f'Candidate source differs from paper: {key} {gse}', 'candidate_provenance')
        if key not in block_cache:
            block_cache[key] = json.loads(Path(candidate['text_path']).read_text())['blocks']
        blocks = block_cache[key]
        for passage in candidate['passages']:
            index = passage['block_index']
            if index is None:
                condition = any(passage['text'] == e.get('context') for e in candidate['accession_evidence'])
            else:
                condition = 0 <= index < len(blocks) and passage['text'] in blocks[index]['text'] and \
                            passage['section'] == blocks[index]['section'] and passage['label'] == blocks[index]['label']
            require(condition, f'Passage absent from source: {key} {gse}', 'passage_provenance')
        require(candidate['triage_class'] in ('LIKELY', 'POSSIBLE', 'UNLIKELY'), 'Invalid triage', 'triage')
        counts[candidate['triage_class']] += 1
        for modality in candidate['modalities']:
            by_modality[modality][candidate['triage_class']] += 1
    require(dict(counts) == summary['triage_counts'], 'Triage counts differ', 'counts')
    require(dict(by_modality) == summary['triage_by_modality'], 'Modality counts differ', 'counts')
    require(len(top) == min(30, counts['LIKELY']), 'Top-pair count differs', 'top_candidates')
    for candidate in top:
        require(candidate['triage_class'] == 'LIKELY' and
                candidate == candidate_map[(candidate['paper_id'], candidate['gse_id'])],
                'Top candidate differs from fetched candidate', 'top_candidates')
    result = {'status': 'PASS', 'checks': dict(checks), 'total_checks': sum(checks.values()),
              'papers': len(papers), 'candidates': len(candidates), 'top_pairs': len(top),
              'limits': 'Validates source linkage/counts; triage is a literature-screening judgment, not benchmark acceptance.'}
    write_json(output / 'validation.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--atlas', type=Path, default=ATLAS)
    args = parser.parse_args()
    try:
        result = validate(args.output, args.atlas)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        logging.error('Retrieval validation failed: %s', exc)
        raise SystemExit(1) from exc
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
