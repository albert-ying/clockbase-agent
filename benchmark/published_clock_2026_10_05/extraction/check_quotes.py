"""Check extractor A's schema, literal cached quotations and copied statistics."""
import argparse
import hashlib
import json
import logging
import re
import xml.etree.ElementTree as ET
from collections import Counter
from decimal import Decimal
from html.parser import HTMLParser
from pathlib import Path


ROOT = Path('/Entropy/kying/clockbase-benchmark-2026-10-05')
FIELDS = {'id', 'tier', 'source', 'gse', 'modality', 'tissue', 'clock_reported',
          'clock_clockbase', 'measure', 'case_label', 'control_label', 'case_definition',
          'control_definition', 'n_case', 'n_control', 'direction', 'significant',
          'threshold', 'p_value', 'effect', 'covariates', 'extraction'}


def normalized(text):
    return re.sub(r'\s+', ' ', text).strip()


class DocumentParser(HTMLParser):
    """Keep visible paragraph and cell text in cached HTML."""
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.root = ET.Element('document')
        self.stack = [self.root]

    def handle_starttag(self, tag, attrs):
        node = ET.SubElement(self.stack[-1], tag, dict(attrs))
        if tag not in {'area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input',
                       'link', 'meta', 'param', 'source', 'track', 'wbr'}:
            self.stack.append(node)

    def handle_startendtag(self, tag, attrs):
        ET.SubElement(self.stack[-1], tag, dict(attrs))

    def handle_endtag(self, tag):
        for i in range(len(self.stack) - 1, 0, -1):
            if self.stack[i].tag == tag:
                self.stack = self.stack[:i]
                break

    def handle_data(self, data):
        node = self.stack[-1]
        if len(node):
            node[-1].tail = (node[-1].tail or '') + data
        else:
            node.text = (node.text or '') + data


def document(path):
    if path.suffix == '.xml':
        return ET.parse(path).getroot()
    parser = DocumentParser()
    parser.feed(path.read_text())
    return parser.root


def text_of(node):
    return normalized(''.join(node.itertext()))


def cached_texts(path):
    """Expose cached text blocks, paragraphs and whitespace-delimited table rows."""
    if path.suffix == '.json':
        data = json.loads(path.read_text())
        return [normalized(b['text']) for b in data['blocks']]
    if path.suffix in {'.xml', '.html'}:
        tree = document(path)
        texts = [text_of(n) for n in tree.iter() if n.tag in {'p', 'td', 'th', 'caption'}]
        texts += [' '.join(text_of(c) for c in row if c.tag in {'td', 'th'})
                  for row in tree.findall('.//tr')]
        return texts
    return [normalized(path.read_text())]


def numbers(text):
    """Read decimal and scientific notation without rounding reported values."""
    text = text.translate(str.maketrans('−–⁻⁰¹²³⁴⁵⁶⁷⁸⁹', '--'+ '-0123456789'))
    text = text.replace('\u00ad', '').replace('×', 'x')
    text = re.sub(r'(?<=\d),(?=\d{3}(?:\D|$))', '', text)
    pattern = (r'(?<![\w.])([+-]?\s*(?:\d+(?:\.\d*)?|\.\d+))'
               r'(?:\s*(?:[eE]\s*([+-]?\s*\d+)|[xX]\s*10\s*([+-]?\s*\d+)))?'
               r'(?![\w.])')
    values = []
    for match in re.finditer(pattern, text):
        mantissa, exponent, power = match.groups()
        value = Decimal(re.sub(r'\s', '', mantissa))
        if exponent is not None or power is not None:
            value *= Decimal(10) ** int(re.sub(r'\s', '', exponent or power))
        elif mantissa.strip() == '10':
            # Superscript exponents in visible XML become 10-8.
            tail = text[match.end():]
            exp = re.match(r'([+-]\d+)(?!\w)', tail)
            if exp:
                value = Decimal(10) ** int(exp[1])
        values.append(value)
    return values


def check(root, files):
    """Verify quotes against their cache hashes and check every schema statistic."""
    evidence = json.loads((root / 'evidence_manifest_a.json').read_text())
    sources = {}
    ids = set()
    rows = []
    failures = []
    fragments = 0
    statistics = 0
    for name in files:
        for line, payload in enumerate((root / name).read_text().splitlines(), 1):
            try:
                row = json.loads(payload)
                if set(row) != FIELDS:
                    raise ValueError('Record fields differ from GROUND_TRUTH_SCHEMA.md')
                if row['id'] in ids:
                    raise ValueError('Duplicate record ID')
                ids.add(row['id'])
                if row['tier'] not in {'T1', 'T2'} or row['direction'] not in {'higher', 'lower', 'none'}:
                    raise ValueError('Invalid tier or direction')
                if row['measure'] not in {'age', 'acceleration_residual', 'acceleration_difference', 'rate', 'other'}:
                    raise ValueError('Invalid measure')
                if set(row['source']) != {'doi', 'pmid', 'pmcid', 'citation', 'location', 'quote', 'url', 'retrieved'}:
                    raise ValueError('Source fields differ from the schema')
                if set(row['extraction']) != {'extractors', 'agreement', 'notes'}:
                    raise ValueError('Extraction fields differ from the schema')
                entries = evidence['records'][row['id']]
                if normalized(row['source']['quote']) != normalized('\n'.join(e['text'] for e in entries)):
                    raise ValueError('Source quote differs from its evidence fragments')
                for entry in entries:
                    path = Path(entry['file'])
                    if path not in sources:
                        digest = hashlib.sha256(path.read_bytes()).hexdigest()
                        sources[path] = digest, cached_texts(path)
                    digest, texts = sources[path]
                    if digest != entry['sha256']:
                        raise ValueError(f'Cache identity changed: {path}')
                    quote = normalized(entry['text'])
                    if not quote or not any(quote in text for text in texts):
                        raise ValueError(f'Quote absent from {path.name}: {quote[:100]}')
                    fragments += 1
                values = numbers(row['source']['quote'])
                copied = [(key, row[key]) for key in ['n_case', 'n_control', 'p_value'] if row[key] is not None]
                if row['effect'] is not None:
                    if set(row['effect']) != {'value', 'unit', 'type'}:
                        raise ValueError('Effect fields differ from the schema')
                    copied.append(('effect.value', row['effect']['value']))
                for key, value in copied:
                    if isinstance(value, bool) or not isinstance(value, (int, float)):
                        raise ValueError(f'Invalid numeric value: {key}')
                    if Decimal(str(value)) not in values:
                        raise ValueError(f'{key}={value} is absent from the quote')
                    statistics += 1
                for key in ['n_case', 'n_control']:
                    if row[key] is not None and (not isinstance(row[key], int) or row[key] < 1):
                        raise ValueError(f'Invalid reported count: {key}')
                if row['p_value'] is not None and not 0 <= row['p_value'] <= 1:
                    raise ValueError('P value outside [0, 1]')
                if row['tier'] == 'T2' and (copied or row['significant'] is not None or row['threshold'] is not None):
                    raise ValueError('Expected direction acquired a reported statistic')
                rows.append(row)
            except (OSError, ValueError, KeyError, TypeError, ET.ParseError) as exc:
                failures.append(f'{name}:{line}: {exc}')
    if failures:
        raise ValueError('\n'.join(failures))
    return {'result': 'passed', 'records': len(rows), 'quote_fragments': fragments,
            'statistics_checked': statistics, 'cache_files': len(sources),
            'tiers': dict(Counter(r['tier'] for r in rows)),
            'modalities': dict(Counter(r['modality'] for r in rows)),
            'clock_clockbase_set': sum(r['clock_clockbase'] is not None for r in rows)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT / 'extraction')
    parser.add_argument('--files', nargs='+', default=['extractor_a.jsonl', 'curated_normalized.jsonl'])
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    try:
        print(json.dumps(check(args.root, args.files), indent=2))
    except (OSError, ValueError, KeyError) as exc:
        logging.error('Quote verification failed: %s', exc)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
