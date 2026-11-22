"""Load GEO descriptions and construct prompts without published results."""
import csv
import gzip
import hashlib
import logging
import re
from pathlib import Path

from audit_system import paths


MODALITIES = {
    'human DNAm': ('gse_me_hs.csv', 'GEO_hs_dnam'),
    'human RNA-seq': ('gse_rna_hs.csv', 'GEO_hs_rna'),
    'mouse RNA-seq': ('gse_rna_mm.csv', 'GEO_mm_rna'),
}
CLOCKS = {
    'Horvath pan-tissue age': 'Estimates chronological age from DNA methylation across tissues.',
    'Hannum': 'Estimates chronological age from blood DNA methylation.',
    'Consensus Hannum and Horvath': 'Combines Hannum and Horvath DNA methylation age estimates.',
    'PhenoAge': 'Estimates a phenotypic age target from DNA methylation.',
    'GrimAge': 'Estimates a mortality-associated age score from DNA methylation.',
    'GrimAge2': 'Estimates a mortality-associated age score from DNA methylation using the second version.',
    'DunedinPACE': 'Estimates the pace of biological aging from DNA methylation.',
    'Skin and blood': 'Estimates chronological age from DNA methylation in skin and blood tissues.',
    'Lin': 'Estimates chronological age from DNA methylation.',
    'Zhang': 'Estimates chronological age from DNA methylation.',
    'tAge mouse': 'Estimates chronological age from mouse gene expression.',
    'tAge rodent': 'Estimates chronological age from rodent gene expression.',
    'tAge multi-species': 'Estimates chronological age from gene expression across species.',
    'tAge mouse mortality': 'Estimates a mortality-associated score from mouse gene expression.',
    'tAge rodent mortality': 'Estimates a mortality-associated score from rodent gene expression.',
    'tAge multi-species mortality': 'Estimates a mortality-associated score from gene expression across species.',
}
MEASURES = {
    'age': 'Clock output.',
    'acceleration_residual': 'Residual of clock output regressed on chronological age.',
    'acceleration_difference': 'Clock output minus chronological age.',
    'rate': 'Pace of biological aging.',
}
RESULT_WORDS = re.compile(
    r'\b(?:accelerat\w*|decelerat\w*|rejuvenat\w*|higher|lower|older|younger|'
    r'significan\w*|non[- ]significant|no difference|not different|(?:did |does |do )?not differ|age reversal|'
    r'increas\w*|decreas\w*|elevat\w*|reduc\w*|unchanged|similar)\b|'
    r'\b(?:slower|faster|progeria|premature ag[ei]ng|no (?:association|effect|change))\b|'
    r'\b(?:p|p[- ]value|fdr)\s*(?:[<>=≤≥]|of\s+\d)', re.I)
TASK = """Predict the case group relative to the control group for the specified clock measure in this series.
Use the study description and your prior knowledge only. Do not access sample data, tools or external sources.
Predict higher, lower, or none (no difference), and whether the difference would be significant at P < 0.05.
Confidence is your probability (0 to 1) that the direction prediction is correct.
Say whether you recognize a publication reporting this clock result for this GEO series; if yes, identify it.
Return JSON only:
{"prediction": "higher|lower|none", "significance": true,
 "confidence": 0.5, "recognizes_publication": false, "publication": null}
Use a boolean for significance and recognition. If recognition is true, publication must identify the paper.
Choose each value independently; the example values do not describe this study."""


class BlindingError(ValueError):
    def __init__(self, message, flags):
        super().__init__(message)
        self.flags = flags


def _normal(text):
    return re.sub(r'\s+', ' ', str(text)).strip().casefold()


def _receipt(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return {'path': paths.portable(path), 'sha256': digest.hexdigest(), 'mtime_ns': path.stat().st_mtime_ns}


def clock_description(record):
    """Map reported names to neutral names, omitting result wording in clock labels."""
    text = ' '.join(str(record.get(key) or '') for key in ('clock_clockbase', 'clock_reported')).lower()
    compact = re.sub(r'[^a-z0-9]', '', text)
    if 'consensus' in text and 'hannum' in text and 'horvath' in text:
        name = 'Consensus Hannum and Horvath'
    elif 'tage' in compact or ('mortality' in text and 'rna' in str(record.get('modality', '')).lower()):
        if 'multi' in text or 'mammal' in text:
            species = 'multi-species'
        elif 'rodent' in text:
            species = 'rodent'
        elif 'mouse' in text or 'murine' in text or record.get('modality') == 'mouse RNA-seq':
            species = 'mouse'
        else:
            raise ValueError(f'No neutral tAge species description for {record.get("clock_reported")!r}')
        name = 'tAge ' + species + (' mortality' if 'mort' in text or 'lifespan' in text else '')
    elif 'dunedinpace' in compact:
        name = 'DunedinPACE'
    elif re.search(r'grimage(?:v)?2(?!\d)', compact):
        name = 'GrimAge2'
    elif 'grimage' in compact:
        name = 'GrimAge'
    elif 'phenoage' in compact:
        name = 'PhenoAge'
    elif ('skin' in text and 'blood' in text) or re.search(r'horvath(?:v)?2(?!\d)', compact):
        name = 'Skin and blood'
    elif 'horvath' in text:
        name = 'Horvath pan-tissue age'
    else:
        name = next((key for key in ('Hannum', 'Lin', 'Zhang')
                     if re.search(rf'\b{key.lower()}(?:\b|\d)', text)), None)
    if name is None:
        raise ValueError(f'No neutral clock description for {record.get("clock_reported")!r}')
    return name, CLOCKS[name]


class StudyDescriptions:
    def __init__(self, geo=None):
        self.geo = Path(geo) if geo is not None else paths.GEO
        self.tables = {}

    def get(self, record):
        """Read the modality index, falling back to the series matrix header."""
        modality, gse = record['modality'], record['gse']
        if modality not in MODALITIES:
            raise ValueError(f'Unsupported S4 modality: {modality}')
        if not re.fullmatch(r'GSE\d+', gse):
            raise ValueError(f'Invalid GEO accession: {gse!r}')
        index, folder = MODALITIES[modality]
        table_path = self.geo / index
        if modality not in self.tables:
            rows = {}
            if table_path.exists():
                try:
                    with table_path.open(newline='') as handle:
                        for row in csv.DictReader(handle):
                            rows[row.get('accession')] = row
                    provenance = _receipt(table_path)
                except (OSError, csv.Error) as error:
                    raise RuntimeError(f'Cannot read GEO description index {table_path}: {error}') from error
            else:
                provenance = None
            self.tables[modality] = rows, provenance
        rows, provenance = self.tables[modality]
        row = rows.get(gse, {})
        study = {'title': row.get('title') or '', 'summary': row.get('description') or '',
                 'overall_design': '', 'n_samples': _count(row.get('samples')), 'sources': []}
        if row:
            study['sources'].append(provenance)
        matrix_dir = self.geo / folder / 'series_matrix'
        if not study['title'] or not study['summary'] or study['n_samples'] is None:
            candidates = sorted(set(matrix_dir.glob(f'{gse}_series_matrix.txt*')) |
                                set(matrix_dir.glob(f'{gse}-*_series_matrix.txt*')))
            sample_ids = set()
            for path in candidates:
                header, gsms = _matrix_header(path)
                sample_ids.update(gsms)
                for key in ('title', 'summary', 'overall_design'):
                    if not study[key] and header[key]:
                        study[key] = header[key]
                study['sources'].append(_receipt(path))
            if study['n_samples'] is None and sample_ids:
                study['n_samples'] = len(sample_ids)
        if not study['title'] or not (study['summary'] or study['overall_design']):
            raise RuntimeError(f'Cannot proceed - required data file not found: {table_path} '
                               f'(description for {gse}; fallback {matrix_dir}/{gse}_series_matrix.txt.gz)')
        return study


def _count(value):
    try:
        count = int(value)
    except (TypeError, ValueError):
        return None
    return count if count > 0 else None


def _matrix_header(path):
    header, samples = {'title': [], 'summary': [], 'overall_design': []}, set()
    keys = {'!Series_title': 'title', '!Series_summary': 'summary', '!Series_overall_design': 'overall_design'}
    opener = gzip.open if path.suffix == '.gz' else open
    try:
        with opener(path, 'rt') as handle:
            for line in handle:
                if line.startswith('!series_matrix_table_begin'):
                    break
                fields = next(csv.reader([line], delimiter='\t'))
                if fields and fields[0] in keys:
                    header[keys[fields[0]]].extend(fields[1:])
                elif fields and fields[0] == '!Sample_geo_accession':
                    samples.update(fields[1:])
    except (OSError, EOFError, csv.Error) as error:
        raise RuntimeError(f'Cannot read GEO series header {path}: {error}') from error
    return {key: ' '.join(values) for key, values in header.items()}, samples


def forbidden_fields(record):
    """Source text and result values forbidden in study-specific prompt content."""
    source = record.get('source') or {}
    fields = {'source.quote': source.get('quote'), 'source.title': source.get('title'),
              'source.paper_title': source.get('paper_title'), 'source.citation': source.get('citation'),
              'threshold': record.get('threshold')}
    for key in ('p_value', 'statistic', 'reported_statistic', 'result', 'result_text'):
        fields[key] = record.get(key)
    effect = record.get('effect') or {}
    if isinstance(effect, dict):
        fields['effect.value'] = effect.get('value')
    return {key: _normal(value) for key, value in fields.items() if value is not None and str(value).strip()}


def check_blinding(stimulus, record):
    """Check stimulus, excluding the fixed task's answer vocabulary and P < 0.05 target."""
    content = _normal(stimulus)
    issues = []
    for key, value in forbidden_fields(record).items():
        pattern = rf'(?<!\w){re.escape(value)}(?!\w)'
        if re.search(pattern, content):
            issues.append(key)
    # Direction/significance booleans are categorical answer options, not unique source text.
    # Detect their result wording in the study content instead of banning the response schema.
    if RESULT_WORDS.search(stimulus):
        issues.append('result_wording')
    return sorted(set(issues))


def _blind_text(text, record):
    parts = re.split(r'(?<=[.!?])\s+|\n+', str(text))
    kept = [part for part in parts if part.strip() and not check_blinding(part, record)]
    return ' '.join(kept), len(kept) != len([part for part in parts if part.strip()])


def build_prompt(record, study):
    """Use explicit input fields only; keep revealing group definitions with a flag."""
    flags = []
    metadata = {}
    for key in ('title', 'summary', 'overall_design'):
        metadata[key], redacted = _blind_text(study.get(key) or '', record)
        if redacted:
            flags.append(f'geo_{key}_redacted')
    name, description = clock_description(record)
    measure = MEASURES.get(record['measure'])
    if measure is None:
        raise ValueError(f'No neutral measure description for {record["measure"]!r}')
    n_case, n_control = _count(record.get('n_case')), _count(record.get('n_control'))
    n_samples = study.get('n_samples')
    if n_samples is None:
        flags.append('series_sample_count_unavailable')
    case, control = record['case_definition'], record['control_definition']
    group_issues = check_blinding(f'{case}\n{control}', record)
    if group_issues:
        flags.append('case_control_answer_leakage')
    neutral = (f'GEO series: {record["gse"]}\nModality: {record["modality"]}\n'
               f'Title: {metadata["title"] or "[withheld]"}\n'
               f'Summary: {metadata["summary"] or "[withheld]"}\n'
               f'Overall design: {metadata["overall_design"] or "[not supplied]"}\n'
               f'Clock: {name}. {description}\nMeasure: {measure}')
    issues = check_blinding(neutral, record)
    # A short numeric effect can also be a sample count. Omit that count rather than leaking its value.
    counts = []
    for label, value in (('Series samples', n_samples), ('Case samples', n_case), ('Control samples', n_control)):
        line = f'{label}: {value if value is not None else "not reported"}'
        if check_blinding(line, record):
            line = f'{label}: [withheld]'
            flags.append('sample_count_result_collision')
        counts.append(line)
    stimulus = neutral + '\n' + '\n'.join(counts) + f'\nCase definition: {case}\nControl definition: {control}'
    # Source quotes/statistics in definitions cannot be sent; retain a failed row for the benchmark.
    hard_group_issues = [issue for issue in group_issues if issue != 'result_wording']
    issues.extend(hard_group_issues)
    if issues:
        raise BlindingError(f'Prompt blinding failed for {record["id"]}: '
                            f'{", ".join(sorted(set(issues)))}', sorted(set(flags)))
    if not metadata['title'] and not metadata['summary'] and not metadata['overall_design']:
        flags.append('geo_description_fully_redacted')
    logging.info('Prepared %s (%s)', record['id'], ', '.join(flags) or 'no prompt flags')
    prompt = stimulus + '\n\n' + TASK
    return {'prompt': prompt, 'stimulus': stimulus, 'flags': sorted(set(flags)),
            'blinding_issues': group_issues, 'sources': study.get('sources', []),
            'estimated_prompt_tokens': (len(prompt) + 3) // 4}
