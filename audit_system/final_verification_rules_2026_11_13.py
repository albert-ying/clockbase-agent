"""Deterministic output, single-cell unit and title-pair checks on recovered arms.

Only deposited GEO headers, GEO sample tables and empirical predictions are read.
Software fixtures in the tests are synthetic and are never analysis inputs.
"""
import csv
import gzip
import re
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system import design_estimand_2026_10_10 as de, paths
from audit_system.membership_recovery_2026_09_27 import MODALITIES, sha256

POLICY = Path(__file__).parent / 'policy'
UNITS = POLICY / 'clock_units_2026_11_13.csv'
HEADERS = POLICY / 'geo_series_headers_2026_11_13.csv'
OUTPUT_REASON = 'clock outputs outside the plausible range'
CELL_REASON = 'samples are single cells; individuals not recorded'
OUTPUT_CODE = 'CLOCK_OUTPUTS_OUTSIDE_PLAUSIBLE_RANGE'
CELL_CODE = 'SINGLE_CELLS_AS_INDEPENDENT_SAMPLES'
LOWER_YEARS, UPPER_YEARS, FLOOR_SHARE = -100., 250., .9
AGE_ERROR_YEARS, AGE_COVERAGE = 100., .5
# biolearn/model.py anti_trafo: (1 + adult_age) * exp(x) - 1 for x < 0.
# The same transform is used by Horvathv1, Horvathv2 and PEDBE. Its limit is -1.
FLOOR_ATOL = 1e-8
SINGLE = re.compile(r'\bsingle[ -]?(?:cells?|nucle(?:us|i|ar))\b|\bscRNA[ -]?(?:seq)?\b|\bsnRNA[ -]?seq\b|\bSmart[ -]?seq\w*\b|\b10x\b|\bDrop[ -]?seq\b|\bCEL[ -]?seq\w*\b|\bMARS[ -]?seq\b|\bFluidigm\s+C1\b', re.I)
BULK = re.compile(r'\bbulk\b|\bpooled\s+(?:cells|nuclei)\b|\bstranded\s+(?:mRNA|(?:total\s+)?RNA)|\bSMART[ -]?Seq\s+Stranded\b|\bTruSeq\b|\bwhole\s+(?:preimplantation\s+)?embryos?\b', re.I)
CELL_TITLE = re.compile(r'\bsingle[ -]?cell\b|\bcell[ _-]*\d+\b|(?:^|[_ -])[ACGT]{10,}(?:-\d+)?(?:$|[_ ])|(?:SS2|Smart[ -]?seq|SNK\d+|plate).*[_ -][A-P](?:0?[1-9]|1\d|2[0-4])$', re.I)
WELL_TITLE = re.compile(r'[_ -][A-P](?:0?[1-9]|1\d|2[0-4])$|\b(?:BC|barcode)[_-]?\d+\b', re.I)
PLATE_PROTOCOL = re.compile(r'\b(?:Smart[ -]?seq(?:2)?|SS2|Fluidigm\s+C1|C1\s+(?:system|platform)|CEL[ -]?seq\w*|MARS[ -]?seq\w*)\b', re.I)
CELL_INDEX = re.compile(r'\bsingle[ -]+(?:cells?|nucle(?:us|i))\b|\b(?:cell|nucleus|nuclei)[ -]*(?:id|index|#)?[ -]*\d+\b|\b(?:cell|nucleus|nuclei)[ -]+[A-Z]\d+\b', re.I)
CELL_FIELD = re.compile(r'^(?:single[ _-]*)?(?:cell|nucleus|nuclei)[ _-]*(?:id|index|barcode)$|^(?:plate[ _-]*)?well(?:[ _-]*(?:id|position))?$', re.I)
MULTICELL_LIBRARY = re.compile(r'\b(?:10\s*x|InDrop|Drop[ -]?seq|droplet|Visium|spatial|whole\s+(?:preimplantation\s+)?embryo|bulk|SMART[ -]?Seq\s+Stranded)\b', re.I)
BEFORE = re.compile(r'\b(?:pre[ -]?(?:training|treatment|intervention)?|before|baseline)\b', re.I)
AFTER = re.compile(r'\b(?:post[ -]?(?:training|treatment|intervention)?|after|follow[ -]?up)\b', re.I)
SUBJECT_TOKEN = re.compile(r'^(\d+)[A-Za-z](?=[ _-])|\b((?:subject|patient|donor|participant|individual|volunteer|animal|mouse|mice|rat|person|P)[ _#-]*\d+)(?=[ _-]|$)', re.I)


@lru_cache(maxsize=None)
def policy_table(path, index):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {path}')
    try:
        frame = pd.read_csv(path, dtype=str).fillna('')
    except (OSError, ValueError) as error:
        raise ValueError(f'Cannot read verification policy {path}: {error}') from error
    return frame.set_index(index).to_dict('index')


FLOORS = {clock: float(row['lower_transform_bound'])
          for clock, row in policy_table(str(UNITS), 'name').items() if row.get('lower_transform_bound')}


def metadata(geo, modality, gse):
    """Sample titles and public series header evidence, without analysis traces."""
    base = Path(geo) / MODALITIES[modality]['folder']
    files, frames = [], []
    for directory in ('processed_metadata_preharmonize', 'processed_metadata'):
        path = base / directory / f'{gse}.csv'
        frame = de.read_frame(path)
        if frame is not None:
            frames.append(frame)
            files.append(str(path))
    headers = dict(policy_table(str(HEADERS), 'gse_id').get(gse, {}))
    titles = {}
    fields = {}
    characteristics, sample_protocols = {}, {}
    age_fields = {}
    for frame in frames:
        if 'title' in frame:
            for gsm, value in frame['title'].items():
                if de.clean(value):
                    titles.setdefault(str(gsm), str(value))
        for column in frame:
            if str(column) == 'harmonized.age' or de.AGE_FIELD.fullmatch(str(column).strip().lower()):
                parsed = {}
                for gsm, value in frame[column].items():
                    match = de.YEARS.fullmatch(str(value).strip().lower())
                    if match:
                        parsed[str(gsm)] = float(match.group(1))
                if parsed:
                    age_fields.setdefault(str(column), {}).update(parsed)
            if de.SUBJECT_NAME.match(str(column).strip().lower()):
                fields.setdefault(str(column), {str(g): de.clean(v) for g, v in frame[column].items() if de.clean(v)})
            if column != 'title':
                target = sample_protocols if re.search(r'protocol|library|extract', str(column), re.I) else characteristics
                for gsm, value in frame[column].items():
                    if de.clean(value):
                        target.setdefault(str(gsm), {}).setdefault(str(column), str(value))
    # Local series matrices take precedence over the deposited header snapshot.
    for path in sorted((base / 'series_matrix').glob(f'{gse}*_series_matrix*.txt.gz')):
        try:
            with gzip.open(path, 'rt', errors='replace') as handle:
                matrix = {}
                for row in csv.reader(handle, delimiter='\t'):
                    if row and row[0] == '!series_matrix_table_begin':
                        break
                    if row and row[0].startswith(('!Series_', '!Sample_')):
                        matrix.setdefault(row[0], []).append(row[1:])
            for key, name in (('title', 'title'), ('summary', 'summary'), ('overall_design', 'design')):
                if '!Series_' + key in matrix:
                    headers[name] = ' '.join(v for row in matrix['!Series_' + key] for v in row)
            gsms = matrix.get('!Sample_geo_accession', [[]])[0]
            for gsm, title in zip(gsms, matrix.get('!Sample_title', [[]])[0]):
                titles.setdefault(gsm, title)
            files.append(str(path))
        except (OSError, ValueError) as error:
            raise ValueError(f'Cannot read GEO header {path}: {error}') from error
    return {'headers': headers, 'titles': titles, 'individual_fields': fields,
            'characteristics': characteristics, 'sample_protocols': sample_protocols,
            'age_fields_years': age_fields, 'files': files}


def cell_evidence(meta, case, control):
    titles, headers = meta.get('titles', {}), meta.get('headers', {})
    gsms = case + control
    hits = {name: sorted(set(m.group() for m in SINGLE.finditer(headers.get(name, ''))))
            for name in ('title', 'summary', 'design')}
    assay_hits = {}
    for name in hits:
        text = headers.get(name, '')
        accepted = []
        for sentence in re.split(r'(?<=[.!?])\s+', text):
            for match in SINGLE.finditer(sentence):
                if re.match(r'\s+(?:suspensions?|resolution|phenotyping)\b', sentence[match.end():], re.I):
                    continue
                if name == 'summary':
                    if re.search(r'\b(?:recently|previously|prior studies|previous studies|published studies)\b', sentence, re.I):
                        continue
                    if re.match(r'single', match.group(), re.I) and not re.search(
                            r'RNA|seq|transcriptom|gene.expression', sentence[max(0,match.start()-40):match.end()+60], re.I):
                        continue
                accepted.append(match.group())
        assay_hits[name] = sorted(set(accepted))
    sample_hits = []
    sample_evidence = {}
    plate_context = bool(PLATE_PROTOCOL.search(' '.join(headers.values())))
    library_context = bool(MULTICELL_LIBRARY.search(headers.get('title', '') + ' ' + headers.get('design', '')))
    library_context = library_context or bool(re.search(
        r'\b(?:10\s*x|InDrop|Drop[ -]?seq|droplet|Visium|spatial)\b', headers.get('summary', ''), re.I))
    bulk = bool(BULK.search(headers.get('title', '') + ' ' + headers.get('design', '')))
    for gsm in gsms:
        title = titles.get(gsm, '')
        normalized = title.replace('_', ' ')
        characteristics = meta.get('characteristics', {}).get(gsm, {})
        protocol = ' '.join(meta.get('sample_protocols', {}).get(gsm, {}).values())
        cell_fields = {k: v for k, v in characteristics.items() if CELL_FIELD.fullmatch(k) and de.clean(v)}
        singular = bool(re.search(r'\bsingle[ -]+(?:cell|nucleus)\b', normalized, re.I))
        explicit = singular or (bool(CELL_INDEX.search(normalized)) and any(assay_hits.values()))
        direct_index = bool(re.search(r'\b(?:cell|nucleus)[ -]*(?:[A-Z])?\d+\b', normalized, re.I))
        # Underscores are word separators in deposited sample names. An assay
        # acronym or a library identifier alone does not identify a cell.
        well = bool(WELL_TITLE.search(title))
        plate = bool(PLATE_PROTOCOL.search(normalized + ' ' + protocol)) or plate_context
        barcode = bool(re.search(r'(?:^|[ _-])[ACGT]{10,}(?:-\d+)?(?:$|[ _])', title, re.I))
        indexed_barcode = bool(re.search(r'(?:^|[_ -])(?:MSCB|BC|barcode)[_-]\d+(?:$|[_ -])|\bBC\d+[_-]\d+$', title, re.I))
        dual_index = bool(re.search(r'\bN\d+[_-]N\d+\b', title, re.I))
        nuclear_well = bool(re.search(r'(?:^|[_ -])SNK\d+', title, re.I)) and any(
            re.search(r'\b\d+\s+wel', v, re.I) for v in characteristics.values())
        char_explicit = any(re.fullmatch(r'single[ -]+(?:cell|nucleus)(?:[ -]*\d+)?', str(v).strip(), re.I)
                            for k, v in characteristics.items() if re.search(r'cell|nucle|source|description', k, re.I))
        multicell = bool(MULTICELL_LIBRARY.search(normalized)) or bool(re.search(
            r'\b(?:suspension|pooled|clone|clonal|derived|spheroid|material)\b|\bsingle cells\b', normalized, re.I))
        multicell = multicell or any(MULTICELL_LIBRARY.search(value) for field, value in characteristics.items()
                                     if re.search(r'^source_name$|platform|library|sample type', field, re.I))
        for field, value in characteristics.items():
            if re.fullmatch(r'(?:number[ _-]*of[ _-]*)?cells?|cell[ _-]*(?:number|count)', field, re.I):
                count = re.match(r'\s*(\d+(?:\.\d+)?)', value)
                if count and float(count[1]) > 1:
                    multicell = True
        # Plate coordinates are also used for bulk wells and indexed library
        # preparation. A donor-library, bulk or droplet design cannot establish
        # one-cell wells merely by naming a plate-based amplification kit.
        if bulk or library_context:
            explicit = singular and direct_index
        library_replicate = bool(re.search(r'\brep(?:licate)?[ -]*\d+\b', normalized, re.I))
        donor_library = any(all(g in mapping for g in arm) and len(set(mapping[g] for g in arm)) == len(arm)
                            for mapping in meta.get('individual_fields', {}).values() for arm in (case, control) if arm)
        if library_replicate and not singular:
            explicit = False
        if (library_replicate or donor_library) and not (singular and direct_index):
            explicit = False
        assay_title = bool(re.search(r'\bsingle[ -]+cell[ -]+(?:m?RNA|sequencing|transcriptom)', normalized, re.I))
        indexed_title = bool(re.search(r'(?:^|[_ -])[A-Za-z0-9]+[_-]\d+(?:$|[ ,:_-])', title))
        if assay_title and not (direct_index or indexed_title):
            explicit = False
        if re.search(r'\bembryo\b', normalized, re.I) and not singular:
            multicell = True
        if not (plate_context or any(assay_hits.values())) or library_context or bulk:
            cell_fields = {}
        evidence = []
        if explicit and not multicell:
            evidence.append('sample title identifies a single cell or nucleus or cell index')
        if cell_fields and not multicell:
            evidence.append('sample characteristic records a cell identifier or plate well')
        if char_explicit and not multicell and not library_context:
            evidence.append('sample characteristic identifies a single cell or nucleus')
        coordinate = bool(re.search(r'[_ -][A-P](?:0?[1-9]|1\d|2[0-4])$', title, re.I))
        plate_structure = bool(re.search(r'SS2|SNK\d+|plate', title, re.I)) or (
            bool(re.search(r'\bsingle[ -]+(?:cell|nucle)', ' '.join(headers.values()), re.I))
            and (len(gsms) >= 96 or bool(re.search(r'\b(?:plate|wells?)\b', headers.get('design', ''), re.I))))
        if coordinate and plate and plate_structure and not multicell and not library_context and not bulk:
            evidence.append('plate well with plate-based single-cell protocol')
        if well and nuclear_well and any(assay_hits.values()) and not library_context:
            evidence.append('nucleus plate well with recorded well-plate condition and nuclear assay')
        if dual_index and plate and not library_context:
            evidence.append('dual-index sample structure with plate-based single-cell protocol')
        if indexed_barcode and any(assay_hits.values()) and not library_context:
            evidence.append('indexed cell barcode with a single-cell assay')
        if barcode and (any(assay_hits.values()) or plate):
            evidence.append('sample title is a cell barcode with a single-cell assay')
        if multicell:
            evidence = []
        if evidence:
            sample_hits.append(gsm)
            sample_evidence[gsm] = {'title': title, 'basis': evidence, 'cell_fields': cell_fields}
    wells = [g for g in gsms if WELL_TITLE.search(titles.get(g, ''))]
    explicit = any(assay_hits.values())
    bulk = bool(BULK.search(headers.get('title', '') + ' ' + headers.get('design', '')))
    majority = bool(gsms) and len(sample_hits) > len(gsms) / 2
    summary_only = bool(assay_hits['summary']) and not (assay_hits['title'] or assay_hits['design'])
    ambiguity_reasons = []
    if summary_only and re.search(r'\[RNA[ -]?seq\]',headers.get('title',''),re.I):
        ambiguity_reasons.append('RNA-seq subseries title; single-cell method only in summary')
    if summary_only and re.search(r'RNA[ -]?seq(?:uencing)?\s+of\s+FACS[ -]sorted',headers.get('design',''),re.I):
        ambiguity_reasons.append('design describes RNA-seq of FACS-sorted populations')
    replicates = sum(bool(re.search(r'\brep(?:licate)?[ _-]*\d+|[_ -]rep(?:licate)?[ _-]*\d+',titles.get(g,''),re.I)) for g in gsms)
    if summary_only and gsms and replicates > len(gsms)/2:
        ambiguity_reasons.append('replicate library titles; single-cell method only in summary')
    if any(hits.values()) and not explicit:
        ambiguity_reasons.append('single-cell wording does not identify a sequencing assay of these samples')
    if bulk:
        ambiguity_reasons.append('bulk or whole-embryo sequencing design or library protocol')
    # A series can deposit both bulk and cell measurements. A bulk design wins
    # over a summary mentioning an accompanying single-cell experiment, unless
    # the compared sample titles themselves identify cells.
    # Assay-like plate/barcode titles can resolve mixed series at comparison level.
    strong_titles = sum(bool(re.search(r'SS2|Smart[ -]?seq|SNK\d+|\bsingle[ -]?cell\b|[ACGT]{10,}',titles.get(g,''),re.I)) for g in gsms)
    strong_majority = bool(gsms) and strong_titles > len(gsms)/2
    # Require evidence in each arm. Series text can support a sample-level
    # plate/well interpretation but cannot by itself trigger reclassification.
    arm_evidence = {name: {'samples': len(arm), 'cell_samples': sum(g in sample_evidence for g in arm)}
                    for name, arm in (('case', case), ('control', control))}
    detected = all(v['samples'] and v['cell_samples'] > v['samples'] / 2 for v in arm_evidence.values())
    if not detected and any(hits.values()):
        ambiguity_reasons.append('no majority of per-cell sample identifiers in both comparison arms')
    ambiguous = not detected and (any(hits.values()) or bool(wells) or bool(ambiguity_reasons))
    return {'detected': detected, 'ambiguous': ambiguous, 'series_signals': hits,
            'assay_signals': assay_hits, 'ambiguity_reasons': ambiguity_reasons,
            'series_contexts': {name:[sentence for sentence in re.split(r'(?<=[.!?])\s+',headers.get(name,''))
                                      if SINGLE.search(sentence) or BULK.search(sentence)]
                                for name in ('title','summary','design')},
            'bulk_design_signal': bulk, 'cell_title_samples': len(sample_hits),
            'well_title_samples': len(wells), 'samples': len(gsms),
            'arm_evidence': arm_evidence,
            'sample_evidence': {g: sample_evidence[g] for arm in (case, control)
                                for g in [s for s in arm if s in sample_evidence][:4]},
            'example_titles': [titles[g] for g in (sample_hits or wells or gsms)[:4] if g in titles]}


def individuals(meta, gsms, series):
    """Require a recorded person/animal for every compared cell; never invent units."""
    candidates = []
    for field, mapping in meta.get('individual_fields', {}).items():
        if all(g in mapping for g in gsms):
            candidates.append((len(set(mapping[g] for g in gsms)), field, mapping))
    mapping = series.get('subjects') or {}
    if mapping and all(g in mapping for g in gsms):
        candidates.append((len(set(mapping[g] for g in gsms)), series.get('subject_source') or 'GEO sample field', mapping))
    if candidates:
        # A recorded complete individual field takes precedence over a token
        # inferred from a title; group prefixes must not split a known donor.
        _, field, mapping = max(candidates, key=lambda x: (x[0], x[1]))
        return field, {g: mapping[g] for g in gsms}
    # An explicit mouse/donor token in a cell title records its individual.
    # Retain the preceding group text so independently numbered animals in
    # different age/genotype cohorts cannot become spurious matched pairs.
    title_subjects = {}
    for gsm in gsms:
        title = meta.get('titles', {}).get(gsm, '')
        token = re.search(r'(?:^|[_ -])(?:donor|patient|mouse|animal|subject|m)[_ -]?\d+(?=[_ -]|$)|^P\d+(?=[_ -]|$)', title, re.I)
        if token:
            title_subjects[gsm] = title[:token.end()].strip().casefold()
    if title_subjects and all(g in title_subjects for g in gsms):
        candidates.append((len(set(title_subjects.values())), 'explicit animal or patient token in GEO title', title_subjects))
    animal_libraries = animal_library_individuals(meta)
    if animal_libraries and all(g in animal_libraries for g in gsms):
        candidates.append((len(set(animal_libraries.values())), 'GEO design identifies one animal per library', animal_libraries))
    if not candidates:
        return None, {}
    _, field, mapping = max(candidates, key=lambda x: (x[0], x[1]))
    return field, {g: mapping[g] for g in gsms}


def animal_library_individuals(meta):
    """A series design can explicitly enumerate the animals represented by its libraries.

    Require a complete partition of all sample titles with the recorded group
    sizes; a numbered title alone never establishes an animal identifier.
    """
    titles=meta.get('titles',{})
    if not titles:
        return {}
    design=meta.get('headers',{}).get('design','')
    sex_counts=re.search(r'\b(\d+)\s+females?\s+and\s+(\d+)\s+males?\s*/\s*(?:experimental\s+)?condition',design,re.I)
    animal_codes={g:re.match(r'^(\d+)([fm])\s+([A-Za-z][A-Za-z0-9-]*)[, ]',t,re.I) for g,t in titles.items()}
    if sex_counts and all(animal_codes.values()):
        groups={m[3].casefold() for m in animal_codes.values()}
        tokens={g:''.join(m.groups()).casefold() for g,m in animal_codes.items()}
        counts=Counter((m[3].casefold(),m[2].casefold()) for m in animal_codes.values())
        if len(set(tokens.values()))==len(titles) and all(
                counts[(group,sex)]==int(sex_counts[position])
                for group in groups for sex,position in (('f',1),('m',2))):
            return tokens
    for sentence in re.split(r'(?<=[.!?])\s+',design):
        if not re.search(r'\b(?:animals|mice)\b',sentence,re.I) or re.search(r'\bpool',sentence,re.I):
            continue
        groups=[]
        for part in sentence.strip().rstrip('.').split(','):
            match=re.fullmatch(r'\s*(\d+)\s+([A-Za-z+ -]+?)\s*(?:animals|mice)?\s*',part,re.I)
            if not match:
                groups=[]
                break
            words=[w.casefold() for w in re.findall(r'[A-Za-z]+',match[2]) if w.casefold() not in ('animals','mice','and')]
            groups.append((int(match[1]),words))
        if not groups or sum(n for n,_ in groups)!=len(titles):
            continue
        assigned=[]
        for n,words in groups:
            members=[g for g,t in titles.items() if all(re.search(r'\b'+re.escape(w)+r'\b',t.replace('_',' '),re.I) for w in words)]
            if not words or len(members)!=n:
                assigned=[]
                break
            assigned.extend(members)
        if len(assigned)==len(titles) and len(set(assigned))==len(titles):
            return {g:'animal library: '+titles[g] for g in titles}
    return {}


def title_pairs(meta, case, control):
    """Exact subject bijection and opposite temporal labels in every sample title."""
    if not case or not control or len(case) != len(control) or set(case) & set(control):
        return None
    titles = meta.get('titles', {})
    maps, sides = [], []
    for arm in (case, control):
        mapping, side = {}, []
        for gsm in arm:
            title = titles.get(gsm, '').replace('_', ' ')
            match = SUBJECT_TOKEN.search(title)
            before, after = bool(BEFORE.search(title)), bool(AFTER.search(title))
            if match is None or before == after:
                return None
            token = re.sub(r'\s+', '', next(v for v in match.groups() if v)).casefold()
            mapping[gsm] = token
            side.append('before' if before else 'after')
        if len(set(mapping.values())) != len(arm) or len(set(side)) != 1:
            return None
        maps.append(mapping)
        sides.append(side[0])
    if sides[0] == sides[1] or set(maps[0].values()) != set(maps[1].values()):
        return None
    subjects={**maps[0],**maps[1]}
    conflicts={}
    for field,recorded in sorted(meta.get('individual_fields',{}).items()):
        if not all(g in recorded for g in case+control):
            continue
        inconsistent={token:sorted(set(recorded[g] for g in case+control if subjects[g]==token))
                      for token in sorted(set(subjects.values()))}
        inconsistent={token:values for token,values in inconsistent.items() if len(values)>1}
        if inconsistent:
            conflicts[field]=inconsistent
    return {'subjects': subjects, 'tokens': sorted(maps[0].values()),
            'conflicting_individual_fields':conflicts,
            'case_time': sides[0], 'control_time': sides[1],
            'titles': {g: titles[g] for g in case + control}}


def output_evidence(clock, values, case, control, series=None):
    info = policy_table(str(UNITS), 'name').get(clock)
    if not info or info['output_unit'] != 'years':
        return None
    floor = FLOORS.get(clock)
    ages = dict((series or {}).get('ages', {}))
    meta = (series or {}).get('verification_metadata', {})
    age_fields = meta.get('age_fields_years', {})
    evidence = {'unit': info['output_unit'], 'bounds_without_ages': [LOWER_YEARS, UPPER_YEARS],
                'age_error_years': AGE_ERROR_YEARS, 'age_coverage_threshold': AGE_COVERAGE,
                'age_field': (series or {}).get('age_field'),
                'lower_transform_bound': floor, 'floor_tolerance': FLOOR_ATOL,
                'floor_share_threshold': FLOOR_SHARE, 'arms': {}}
    for name, arm in (('case', case), ('control', control)):
        a = np.asarray([values[g] for g in arm if g in values.index and pd.notna(values[g])], dtype=float)
        if not len(a):
            continue
        mean = float(a[0]) if np.all(a == a[0]) else float(a.mean())
        n_floor = int(np.count_nonzero(np.isclose(a, floor, rtol=0, atol=FLOOR_ATOL))) if floor is not None else 0
        # Coverage is evaluated within each arm, independently of the stricter
        # whole-series coverage required for age-adjusted inferential tests.
        candidates = [(sum(g in mapping for g in arm), field, mapping)
                      for field, mapping in age_fields.items()]
        candidates.append((sum(g in ages for g in arm), (series or {}).get('age_field') or '', ages))
        _, arm_age_field, arm_ages = max(candidates, key=lambda item: (item[0], item[1] == 'harmonized.age', item[1]))
        recorded = [float(arm_ages[g]) for g in arm if g in arm_ages and np.isfinite(arm_ages[g])]
        coverage = len(recorded) / len(arm)
        age_mean = float(np.mean(recorded)) if recorded else None
        use_age = coverage >= AGE_COVERAGE
        outside = abs(mean - age_mean) > AGE_ERROR_YEARS if use_age else not LOWER_YEARS <= mean <= UPPER_YEARS
        evidence['arms'][name] = {'mean': mean, 'samples': len(a), 'samples_at_floor': n_floor,
                                 'recorded_age_samples': len(recorded), 'age_coverage': coverage,
                                 'recorded_age_field': arm_age_field,
                                 'mean_recorded_age_years': age_mean, 'uses_recorded_age': use_age,
                                 'screen': 'mean deviation from recorded age' if use_age else 'absolute range without sufficient recorded ages',
                                 'floor_share': n_floor / len(a),
                                 'outside': outside,
                                 'floor_saturated': floor is not None and n_floor / len(a) >= FLOOR_SHARE}
    evidence['detected'] = any(v['outside'] or v['floor_saturated'] for v in evidence['arms'].values())
    return evidence


def _set_result(record, estimate, series):
    if 'recomputed_before_verification_rules' not in record:
        record['recomputed_before_verification_rules'] = dict(record.get('recomputed', {}))
    record['recomputed'] = estimate['recomputed']
    record['estimand'] = estimate['estimand']
    record['design'] = {**estimate['design'], 'metadata': [paths.portable(p) for p in series.get('files', [])]}
    record['soft_flags'] = sorted(set(record.get('soft_flags', [])) | set(estimate.get('soft', [])))
    record.pop('pending_reason', None)
    if record['recomputed'].get('p') is None:
        record['final_class'] = 'DESCRIPTIVE'
        if record['recomputed'].get('undefined_reason'):
            record['pending_reason'] = [record['recomputed']['undefined_reason']]
    else:
        record['final_class'] = 'PROVISIONAL' if record['soft_flags'] else 'VERIFIED'


def apply(record, values, series, modality, contrast_types=()):
    """Apply before duplicate selection and before the deposited arm decisions/BH."""
    if values is None or not record.get('inputs') or record['final_class'] not in ('VERIFIED', 'PROVISIONAL', 'DESCRIPTIVE'):
        return
    case, control = (record['inputs'][k] for k in ('case_gsms', 'control_gsms'))
    if not case or not control:
        return
    previous_class = record['final_class']
    meta = (series or {}).get('verification_metadata', {})
    evidence = {}
    output = output_evidence(record['clock'], values, case, control, series)
    if output and output['detected']:
        evidence['implausible_outputs'] = output
    cell = cell_evidence(meta, case, control) if modality != 'dnam' else None
    if cell and (cell['detected'] or cell['ambiguous']):
        evidence['single_cells'] = cell
    if cell and cell['detected']:
        field, mapping = individuals(meta, case + control, series or {})
        cell.update(individual_field=field, individuals=mapping)
        record['reasons'] = sorted(set(record['reasons']) | {CELL_CODE})
        if mapping:
            individual_series = {**series, 'subjects': mapping, 'subject_source': field, 'force_individuals': True}
            estimate = de.estimate(values, case, control, individual_series, modality, record['clock'],
                                   (record['comparison_name'], record['condition_name']), contrast_types)
            _set_result(record, estimate, individual_series)
            cell['decision'] = 'per-individual test'
        else:
            summary = {k: v for k, v in record['recomputed'].items() if k in ('mean_case', 'mean_control', 'mean_diff')}
            # Independent units are unknown; the sample counts remain in evidence.
            summary.update(n_case=0, n_control=0, undefined_reason=CELL_REASON)
            estimate = {'recomputed': summary, 'estimand': CELL_REASON, 'soft': [],
                        'design': {'method': 'SINGLE_CELL_INDIVIDUALS_UNRECORDED',
                                   'samples_case': len(case), 'samples_control': len(control)}}
            _set_result(record, estimate, series or {})
            record['pending_reason'] = [CELL_REASON]
            cell['decision'] = CELL_REASON
        record['actions'].append({'action': 'SINGLE_CELL_INDIVIDUAL_UNITS', 'individual_field': field})
    else:
        paired = title_pairs(meta, case, control)
        if paired and paired['conflicting_individual_fields']:
            evidence['title_pair_conflict']=paired
        elif paired:
            evidence['title_pairs'] = paired
            paired_series = {**series, 'subjects': paired['subjects'], 'subject_source': 'one-to-one GEO title subject tokens'}
            estimate = de.estimate(values, case, control, paired_series, modality, record['clock'],
                                   (record['comparison_name'], record['condition_name']), contrast_types)
            if estimate is not None:
                _set_result(record, estimate, paired_series)
                record['actions'].append({'action': 'PAIR_SUBJECT_TOKENS_IN_TITLES', 'tokens': paired['tokens']})
    if output and output['detected']:
        record['reasons'] = sorted(set(record['reasons']) | {OUTPUT_CODE})
        record['final_class'] = 'PENDING_INPUT'
        record['pending_reason'] = [OUTPUT_REASON]
        record['actions'].append({'action': 'WITHHOLD_IMPLAUSIBLE_CLOCK_OUTPUTS', 'reason': OUTPUT_REASON})
        for name in ('p', 't', 'fdr_family', 'cohen_d'):
            record['recomputed'].pop(name, None)
    if evidence:
        record['final_class_before_verification_rules'] = previous_class
        record['verification_rules'] = evidence
        clear_nondefect_cell_flag(record)


def clear_nondefect_cell_flag(record):
    """A recorded animal per GSM is already an independent unit, not a defect.

    This normalization only changes the historical defect flag. It cannot
    change estimates, classes, memberships or multiplicity families.
    """
    cell=record.get('verification_rules',{}).get('single_cells',{})
    mapping=cell.get('individuals') or {}
    inputs=record.get('inputs') or {}
    arms=[inputs.get(name,[]) for name in ('case_gsms','control_gsms')]
    if mapping and all(arm and all(g in mapping for g in arm) and len(set(mapping[g] for g in arm))==len(arm) for arm in arms):
        record['reasons']=[code for code in record['reasons'] if code!=CELL_CODE]


def provenance():
    return {str(path.name): {'sha256': sha256(path), 'mtime_ns': path.stat().st_mtime_ns}
            for path in (UNITS, HEADERS, Path(__file__))}
