"""Outcome-blinded study-design packets; no scoring or scientific approval."""
import json
import argparse
import csv
from collections import Counter
from fractions import Fraction
import hashlib
import io
import math
import os
import numpy as np
from pathlib import Path
import platform
import re
import shutil
import statistics
import sys
import tempfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone


SAMPLE_FIELDS = {
    'gse_id', 'gsm_id', 'title', 'source_name', 'raw.title',
    'raw.source_name_ch1', 'raw.organism_ch1', 'raw.molecule_ch1',
    'raw.platform_id', 'raw_characteristics',
    'harmonized.age', 'harmonized.sex', 'harmonized.subject',
    'harmonized.genetic_information', 'harmonized.treatment',
    'harmonized.compound', 'harmonized.disease', 'harmonized.disease_stage',
    'harmonized.specimen_part', 'harmonized.cell_line',
    'harmonized.developmental_stage', 'harmonized.time',
}
COMPARISON_FIELDS = {
    'comparison_name', 'condition_name', 'condition_category', 'covariate',
    'statistical_test', 'clock',
}

CONTEXT_VERSION = 'study-context-2026-09-16-v1'
SOURCE_TABLE_VERSION = 'clockbase-source-table-v1'
SAMPLE_METHOD_VERSION = 'geo-sample-methods-2026-09-19-v1'
SAMPLE_METHOD_FIELDS = frozenset({
    'treatment_protocol_ch1', 'growth_protocol_ch1', 'extract_protocol_ch1',
    'label_ch1', 'label_protocol_ch1', 'hyb_protocol', 'scan_protocol',
    'data_processing', 'library_strategy', 'library_source', 'library_selection',
})


def is_sample_method_field(field):
    """Accept native GEO method names and numbered repetitions, not substrings."""
    base, sep, suffix = field.partition('.')
    return base in SAMPLE_METHOD_FIELDS and (not sep or re.fullmatch(r'[0-9]+', suffix) is not None)


def method_provenance(samples, *, prefix=''):
    return {
        'version': SAMPLE_METHOD_VERSION,
        'fields_present': sorted({key for row in samples for key in row
                                 if key.startswith(prefix) and is_sample_method_field(key[len(prefix):])}),
        'blinding_scope': 'Field-level withholding of structured outcomes and prior verdicts only. '
                          'Source method narratives may contain published findings; retaining a field '
                          'does not establish outcome blinding, method validity or scientific acceptance.',
    }


def source_json_bytes(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def compact_source_packet(packet):
    """Encode repeated source strings without dropping or normalizing a field."""
    if not isinstance(packet, dict) or 'sample_table' in packet:
        raise ValueError('Source packet is not an unencoded mapping')
    rows = packet.get('samples')
    if not isinstance(rows, list) or not rows or any(not isinstance(row, dict) for row in rows):
        raise ValueError('Source packet has no sample records')
    ids = [row.get('sample_id') for row in rows]
    if any(not isinstance(sid, str) or not sid for sid in ids) or len(set(ids)) != len(ids):
        raise ValueError('Source sample identities are missing or duplicated')
    if any(not isinstance(key, str) or not key or not isinstance(value, str)
           for row in rows for key, value in row.items()):
        raise ValueError('Source sample fields must contain exact strings')
    columns = {}
    for field in sorted(set().union(*(row.keys() for row in rows)) - {'sample_id'}):
        values = [row.get(field) for row in rows]
        present = [value for value in values if value is not None]
        candidates = []
        if len(present) == len(rows) and len(set(present)) == 1:
            candidates.append({'constant': present[0]})
        for prefix in dict.fromkeys(('', os.path.commonprefix(present))):
            suffixes = [value[len(prefix):] if value is not None else None for value in values]
            dictionary = list(dict.fromkeys(value for value in suffixes if value is not None))
            lookup = {value: i for i, value in enumerate(dictionary)}
            common = {'prefix': prefix} if prefix else {}
            candidates.extend([
                {**common, 'values': suffixes},
                {**common, 'dictionary': dictionary,
                 'indices': [lookup[value] if value is not None else -1 for value in suffixes]}])
        columns[field] = min(candidates, key=lambda value: len(source_json_bytes(value)))
    encoded = {key: value for key, value in packet.items() if key != 'samples'}
    encoded['sample_table'] = {'version': SOURCE_TABLE_VERSION, 'sample_ids': ids, 'columns': columns,
                              'packet_sha256': hashlib.sha256(source_json_bytes(packet)).hexdigest()}
    if expand_source_packet(encoded) != packet:
        raise ValueError('Source compaction changed the packet')
    return encoded


def expand_source_packet(encoded):
    """Strictly reconstruct the original packet and verify its declared digest."""
    if not isinstance(encoded, dict) or 'samples' in encoded or not isinstance(encoded.get('sample_table'), dict):
        raise ValueError('Invalid encoded source packet')
    table = encoded['sample_table']
    if (set(table) != {'version', 'sample_ids', 'columns', 'packet_sha256'}
            or table['version'] != SOURCE_TABLE_VERSION):
        raise ValueError('Unknown source table schema')
    ids, columns = table['sample_ids'], table['columns']
    if (not isinstance(ids, list) or not ids or any(not isinstance(sid, str) or not sid for sid in ids)
            or len(set(ids)) != len(ids) or not isinstance(columns, dict)):
        raise ValueError('Invalid source table identities or columns')
    rows = [{'sample_id': sid} for sid in ids]
    for field, column in columns.items():
        if not isinstance(field, str) or not field or field == 'sample_id' or not isinstance(column, dict):
            raise ValueError('Invalid source table field')
        prefix = column.get('prefix', '')
        if not isinstance(prefix, str): raise ValueError('Source prefix must be a string')
        keys = set(column) - {'prefix'}
        if keys == {'constant'} and 'prefix' not in column:
            if not isinstance(column['constant'], str): raise ValueError('Invalid source constant')
            values = [column['constant']] * len(ids)
        elif keys == {'values'}:
            values = column['values']
        elif keys == {'dictionary', 'indices'}:
            dictionary, indices = column['dictionary'], column['indices']
            if (not isinstance(dictionary, list) or not dictionary
                    or any(not isinstance(v, str) for v in dictionary) or len(set(dictionary)) != len(dictionary)
                    or not isinstance(indices, list)
                    or any(type(i) is not int or not -1 <= i < len(dictionary) for i in indices)):
                raise ValueError('Invalid source dictionary or indices')
            values = [dictionary[i] if i != -1 else None for i in indices]
        else:
            raise ValueError('Unknown source column encoding')
        if (not isinstance(values, list) or len(values) != len(ids)
                or any(value is not None and not isinstance(value, str) for value in values)):
            raise ValueError('Source column length or value type differs')
        for row, value in zip(rows, values):
            if value is not None: row[field] = prefix + value
    packet = {key: value for key, value in encoded.items() if key != 'sample_table'}
    packet['samples'] = rows
    if hashlib.sha256(source_json_bytes(packet)).hexdigest() != table['packet_sha256']:
        raise ValueError('Reconstructed source packet digest differs')
    return packet


METHODS_TITLES = {'methods', 'materials and methods', 'material and methods',
                  'experimental procedures', 'experimental methods', 'star methods'}

# Only explicit source labels establish units; titles and replicate numbering
# do not. Unknown labels remain available in the packet for source review.
DESIGN_LABELS = {
    'unit': {'subject id', 'donor id', 'patient id', 'animal id', 'mouse id',
             'individual id'},
    'treatment': {'treatment', 'compound', 'drug'},
    'age': {'age', 'age at collection', 'age at harvest', 'age y', 'age yr', 'age yrs', 'age years'},
    'sex': {'sex', 'gender'},
    'genotype': {'genotype', 'genetic background'},
    'tissue': {'tissue', 'tissue type', 'cell type'},
    'batch': {'batch', 'sequencing batch', 'library batch', 'array batch'},
}
CONTROL_VALUES = {'control', 'untreated', 'vehicle', 'none', 'no treatment',
                  'dmso', 'mock', 'placebo', 'pbs'}
UNRESOLVED_SOURCE_VALUES = {'', 'na', 'n a', 'unknown', 'missing',
                            'not available', 'not reported', 'not applicable'}


def normalized_source_label(value):
    return re.sub(r'[^a-z0-9]+', ' ', value.casefold()).strip()


def structured_characteristics(sample):
    """Parse labeled source characteristics without interpreting sample titles."""
    fields = {}
    for field, text in sample.items():
        if field.startswith('characteristics_ch1') and isinstance(text, str):
            label, separator, value = text.partition(':')
            if separator and value.strip():
                label = normalized_source_label(label)
                fields.setdefault(label, []).append({'field': field, 'value': value.strip(), 'quote': text})
    return fields


def source_metadata_summary(packet, sample_ids=None, *, detailed=False):
    """Calculate descriptive source quantities without clocks, outcomes or verdicts.

    Counts describe assay records, not verified biological units. Numeric values
    are grouped by their literal units; no unit conversion or imputation occurs.
    High-cardinality frequency tables may be omitted from this derived prompt
    summary; the original source packet always retains every value.
    """
    samples = packet['samples']
    by_id = {sample['sample_id']: sample for sample in samples}
    if len(by_id) != len(samples):
        raise ValueError('Metadata summary requires unique source identities')
    ids = list(by_id) if sample_ids is None else list(sample_ids)
    if len(set(ids)) != len(ids) or set(ids) - set(by_id):
        raise ValueError('Metadata summary selection duplicates or invents source identities')
    parsed, labels = {}, set()
    for sid in ids:
        fields = {}
        for field, text in by_id[sid].items():
            if not isinstance(text, str):
                continue
            if field.startswith('characteristics_ch1'):
                label, separator, value = text.partition(':')
                label = normalized_source_label(label) if separator else field
                fields.setdefault(label, []).append((field, value.strip() if separator else text))
            elif field in {'organism_ch1', 'platform_id', 'molecule_ch1'}:
                fields[field] = [(field, text)]
        parsed[sid] = fields
        labels.update(fields)
    summaries = {}
    missing_tokens = {'', 'na', 'n/a', 'n a', 'not available', 'unknown', 'missing', 'not reported'}
    for label in sorted(labels):
        values, absent, missing, conflicts, source_fields = [], 0, 0, 0, set()
        for sid in ids:
            declarations = parsed[sid].get(label, [])
            source_fields.update(field for field, _ in declarations)
            unique = {value for _, value in declarations}
            if not unique:
                absent += 1
            elif len(unique) != 1:
                conflicts += 1
            elif next(iter(unique)).casefold() in missing_tokens:
                missing += 1
            else:
                values.append(next(iter(unique)))
        counts = Counter(values)
        item = {'source_fields': sorted(source_fields), 'observed': len(values),
                'absent': absent, 'declared_missing': missing, 'conflicting': conflicts,
                'distinct_values': len(counts)}
        if detailed or len(counts) <= 12:
            item['counts'] = dict(sorted(counts.items()))
        else:
            item['frequency_table_omitted_from_summary'] = True
        if re.match(r'^(age|time|duration)( |$)', label):
            numeric = {}
            for value in values:
                match = re.fullmatch(r'([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*([^\d]*)', value)
                if match and math.isfinite(float(match[1])):
                    numeric.setdefault(match[2].strip(), []).append(float(match[1]))
            if numeric:
                item['numeric_by_literal_unit'] = [
                    {'unit': unit, 'n': len(v), 'minimum': min(v), 'median': statistics.median(v),
                     'maximum': max(v)} for unit, v in sorted(numeric.items())]
                item['nonnumeric_observed'] = len(values) - sum(len(v) for v in numeric.values())
        summaries[label] = item
    return {'schema': 'source-metadata-summary-v1', 'assay_samples': len(ids),
            'selected_source_sha256': hashlib.sha256(source_json_bytes([by_id[sid] for sid in ids])).hexdigest(),
            'fields': summaries,
            'scope': 'Descriptive source metadata only; counts do not establish biological independence, covariate adjustment or reference truth'}


def design_matrix_diagnostic(design, column_names, sample_ids, contrasts, *,
                             rank_rtol=1e-12, estimability_tolerance=1e-10):
    """Diagnose algebraic design support without making a scientific judgment.

    ``design`` is an n-by-p finite numeric matrix and ``contrasts`` maps named
    contrasts to p-dimensional coefficient vectors. Rank is determined from
    singular values using ``rank_rtol * largest_singular_value``. A contrast
    is estimable when its Euclidean residual after projection onto the row
    space is no larger than ``estimability_tolerance``. The returned
    scientific fields are deliberately unresolved: algebraic estimability
    cannot establish independent biological units or measurement applicability.
    """
    def identities(values, label):
        try:
            values = list(values)
        except TypeError as error:
            raise ValueError(f'{label} must be a finite sequence of identities') from error
        if (not values or any(not isinstance(value, str) or not value.strip() for value in values)
                or len(values) != len(set(values))):
            raise ValueError(f'{label} must contain unique, nonempty identities')
        return values

    if isinstance(rank_rtol, bool) or not isinstance(rank_rtol, (int, float)):
        raise ValueError('rank_rtol must be a finite positive number')
    if isinstance(estimability_tolerance, bool) or not isinstance(estimability_tolerance, (int, float)):
        raise ValueError('estimability_tolerance must be a finite nonnegative number')
    rank_rtol, estimability_tolerance = float(rank_rtol), float(estimability_tolerance)
    if not math.isfinite(rank_rtol) or rank_rtol <= 0:
        raise ValueError('rank_rtol must be a finite positive number')
    if not math.isfinite(estimability_tolerance) or estimability_tolerance < 0:
        raise ValueError('estimability_tolerance must be a finite nonnegative number')

    try:
        matrix = np.asarray(design, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError('Design matrix must be numeric and rectangular') from error
    if matrix.ndim != 2 or not matrix.shape[0] or not matrix.shape[1]:
        raise ValueError('Design matrix must be a nonempty two-dimensional matrix')
    if not np.isfinite(matrix).all():
        raise ValueError('Design matrix must contain only finite values')
    names = identities(column_names, 'Column names')
    samples = identities(sample_ids, 'Sample IDs')
    if len(names) != matrix.shape[1]:
        raise ValueError('Column-name count does not match design matrix columns')
    if len(samples) != matrix.shape[0]:
        raise ValueError('Sample-ID count does not match design matrix rows')
    if not isinstance(contrasts, dict) or not contrasts:
        raise ValueError('Contrasts must be a nonempty mapping of names to vectors')

    _, singular_values, right_vectors = np.linalg.svd(matrix, full_matrices=False)
    largest = float(singular_values[0])
    threshold = rank_rtol * largest
    rank = int(np.count_nonzero(singular_values > threshold)) if largest else 0
    row_basis = right_vectors[:rank]
    contrast_results = {}
    for name, values in contrasts.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError('Contrast names must be nonempty strings')
        try:
            vector = np.asarray(values, dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError(f'Contrast {name!r} must be numeric') from error
        if vector.ndim != 1 or vector.shape[0] != matrix.shape[1]:
            raise ValueError(f'Contrast {name!r} dimension does not match design matrix columns')
        if not np.isfinite(vector).all():
            raise ValueError(f'Contrast {name!r} must contain only finite values')
        projection = row_basis.T @ (row_basis @ vector) if rank else np.zeros_like(vector)
        residual = float(np.linalg.norm(vector - projection))
        contrast_results[name] = {
            'vector': vector.tolist(),
            'row_space_residual_norm': residual,
            'estimability_tolerance': estimability_tolerance,
            'estimable': residual <= estimability_tolerance,
            'status': 'ESTIMABLE' if residual <= estimability_tolerance else 'NOT_ESTIMABLE',
        }
    return {
        'schema': 'clockbase-design-matrix-diagnostic-v1',
        'row_count': int(matrix.shape[0]),
        'column_count': int(matrix.shape[1]),
        'column_names': names,
        'sample_ids': samples,
        'singular_values': singular_values.tolist(),
        'rank_rtol': rank_rtol,
        'rank_threshold': threshold,
        'rank': rank,
        'residual_df': int(matrix.shape[0] - rank),
        'contrasts': contrast_results,
        'scientific_status': 'CANNOT VERIFY',
        'biological_independence_status': 'CANNOT VERIFY',
        'measurement_applicability_status': 'CANNOT VERIFY',
        'scope': 'Algebraic design-matrix diagnostics only; no biological independence, covariate adjustment or measurement applicability is established',
    }


SOURCE_DIAGNOSTIC_CURRENT = 'source-consistency-v3'
SOURCE_DIAGNOSTIC_LEGACY = 'source-consistency-v1'
# Reviewed pre-migration source implementation; no archived code is executed.
SOURCE_DIAGNOSTIC_LEGACY_IMPLEMENTATION = '8701221f9a8fef489baf58e4916c50baed7f05e7f66f81b661fa0e0cc87a3baf'


def source_diagnostic_profile(profile):
    if profile not in {SOURCE_DIAGNOSTIC_CURRENT, SOURCE_DIAGNOSTIC_LEGACY}:
        raise ValueError('Unknown source diagnostic profile')
    return profile


def _verify_requested_comparison_v1(packet, case_ids, control_ids):
    """Bind a two-arm request to complete, exact source-labeled categories.

    Supported syntax is '<source case value> versus <source control value>',
    optionally prefixed by 'Compare', with an optional 'within label: value'
    conjunction of exact source strata. Full characteristic quotations are also
    accepted. Prose and inferred synonyms cannot pass through a substring match.
    """
    result = {'status': 'CANNOT VERIFY', 'profile': 'exact-source-comparison-v1',
              'scope': 'Requested source category labels and complete directional assay membership only',
              'failures': [], 'unresolved': [], 'bindings': [], 'reference_judgment': False}
    if packet is None:
        result['unresolved'].append('Source packet is unavailable')
        return result
    target = packet.get('target', '')
    result['requested_target'] = target
    target_scope = re.split(r'\s+within\s+', target.strip(), maxsplit=1, flags=re.I)
    parts = re.split(r'\s+(?:versus|vs\.?)\s+', re.sub(r'^compare\s+', '', target_scope[0], flags=re.I), flags=re.I)
    if len(parts) != 2 or not all(parts):
        result['unresolved'].append('Requested prose has no explicit two-arm source binding')
        return result
    normalize = lambda text: ' '.join(text.casefold().split())
    strata = []
    if len(target_scope) == 2:
        for clause in re.split(r'\s+and\s+', target_scope[1], flags=re.I):
            label, separator, value = clause.partition(':')
            if not separator or not label.strip() or not value.strip():
                result['unresolved'].append('Requested stratum requires an explicit source label and value')
                return result
            strata.append((normalized_source_label(label), normalize(value)))
        if len({label for label, _ in strata}) != len(strata):
            result['unresolved'].append('Requested stratum labels are repeated or conflicting')
            return result
    result['strata'] = [{'label': label, 'value': value} for label, value in strata]
    case_label, control_label = map(normalize, parts)
    if case_label == control_label:
        result['unresolved'].append('Requested case and control labels are identical')
        return result
    samples = packet.get('samples', [])
    by_id = {sample['sample_id']: sample for sample in samples}
    members = case_ids + control_ids
    if len(by_id) != len(samples) or len(set(members)) != len(members) or set(members) - set(by_id):
        result.update(status='FAIL', failures=['Requested comparison has invalid sample membership'])
        return result
    parsed = {sid: structured_characteristics(sample) for sid, sample in by_id.items()}
    labels = sorted({label for fields in parsed.values() for label in fields})
    for label in labels:
        selected = {'case': [], 'control': []}
        quotes = []
        ambiguous = False
        for sid, fields in parsed.items():
            declarations = fields.get(label, [])
            values = {normalize(item['value']) for item in declarations}
            terms = values | {normalize(item['quote']) for item in declarations}
            terms |= {normalize(label + ': ' + item['value']) for item in declarations}
            matched = [arm for arm, term in (('case', case_label), ('control', control_label)) if term in terms]
            if matched and (len(values) != 1 or len(matched) != 1):
                ambiguous = True
            elif matched:
                source_strata = [{normalize(item['value']) for item in fields.get(key, [])}
                                 for key, _ in strata]
                if any(len(values) != 1 for values in source_strata):
                    ambiguous = True
                    continue
                if any(value not in values for (_, value), values in zip(strata, source_strata)):
                    continue
                selected[matched[0]].append(sid)
                quotes.extend({'sample_id': sid, 'field': item['field'], 'quote': item['quote']}
                              for item in declarations)
                quotes.extend({'sample_id': sid, 'field': item['field'], 'quote': item['quote']}
                              for key, _ in strata for item in fields[key])
        if not ambiguous and selected['case'] and selected['control']:
            result['bindings'].append({'source_label': label, 'case_ids': selected['case'],
                                       'control_ids': selected['control'], 'evidence': quotes})
    if not result['bindings']:
        result['unresolved'].append('Requested arm labels have no unambiguous matching source characteristic')
        return result
    memberships = {(tuple(sorted(item['case_ids'])), tuple(sorted(item['control_ids']))) for item in result['bindings']}
    if len(memberships) != 1:
        result['unresolved'].append('Requested labels map to different source groupings')
        return result
    expected_case, expected_control = next(iter(memberships))
    if set(case_ids) != set(expected_case) or set(control_ids) != set(expected_control):
        result['status'] = 'FAIL'
        result['failures'].append('Assigned arms omit, swap or add samples relative to the source-defined request')
    else:
        result['status'] = 'PASS'
    return result


def verify_requested_comparison(packet, case_ids, control_ids, *, diagnostic_profile=SOURCE_DIAGNOSTIC_CURRENT):
    """Bind a two-arm request to complete, exact source-labeled categories.

    Supported syntax is '<source case value> versus <source control value>',
    optionally prefixed by 'Compare', with an optional 'within label: value'
    conjunction of exact source strata. Full characteristic quotations are also
    accepted. Prose and inferred synonyms cannot pass through a substring match.
    """
    if source_diagnostic_profile(diagnostic_profile) == SOURCE_DIAGNOSTIC_LEGACY:
        return _verify_requested_comparison_v1(packet, case_ids, control_ids)
    result = {'status': 'CANNOT VERIFY', 'profile': 'exact-source-comparison-v3',
              'scope': 'Requested source category labels and complete directional assay membership only',
              'failures': [], 'unresolved': [], 'bindings': [], 'reference_judgment': False}
    if packet is None:
        result['unresolved'].append('Source packet is unavailable')
        return result
    target = packet.get('target', '')
    result['requested_target'] = target
    target_scope = re.split(r'\s+within\s+', target.strip(), maxsplit=1, flags=re.I)
    parts = re.split(r'\s+(?:versus|vs\.?)\s+', re.sub(r'^compare\s+', '', target_scope[0], flags=re.I), flags=re.I)
    if len(parts) != 2 or not all(parts):
        result['unresolved'].append('Requested prose has no explicit two-arm source binding')
        return result
    normalize = lambda text: ' '.join(text.casefold().split())
    strata = []
    if len(target_scope) == 2:
        for clause in re.split(r'\s+and\s+', target_scope[1], flags=re.I):
            label, separator, value = clause.partition(':')
            if not separator or not label.strip() or not value.strip():
                result['unresolved'].append('Requested stratum requires an explicit source label and value')
                return result
            strata.append((normalized_source_label(label), normalize(value)))
        if len({label for label, _ in strata}) != len(strata):
            result['unresolved'].append('Requested stratum labels are repeated or conflicting')
            return result
    result['strata'] = [{'label': label, 'value': value} for label, value in strata]
    case_label, control_label = map(normalize, parts)
    if case_label == control_label:
        result['unresolved'].append('Requested case and control labels are identical')
        return result
    samples = packet.get('samples', [])
    by_id = {sample['sample_id']: sample for sample in samples}
    members = case_ids + control_ids
    if len(by_id) != len(samples) or len(set(members)) != len(members) or set(members) - set(by_id):
        result.update(status='FAIL', failures=['Requested comparison has invalid sample membership'])
        return result
    parsed = {sid: structured_characteristics(sample) for sid, sample in by_id.items()}
    # An explicit source-key assertion can contradict an assigned record even
    # when both requested categories do not exist in the requested stratum.
    # This does not interpret titles, synonyms, biological units or arm intent.
    declared = [part.partition(':') for part in parts]
    requested_values = [value if separator else label for label, separator, value in declared]
    target_unknown = any(normalized_source_label(value) in UNRESOLVED_SOURCE_VALUES
                         for value in requested_values + [value for _, value in strata])
    explicit_label = None
    if (all(separator and label.strip() and value.strip()
            for label, separator, value in declared)
            and normalized_source_label(declared[0][0]) == normalized_source_label(declared[1][0])):
        explicit_label = normalized_source_label(declared[0][0])
        contradictions = []
        for (arm, ids), (_, _, requested) in zip((('case', case_ids), ('control', control_ids)), declared):
            for sid in ids:
                for label, expected in [(explicit_label, normalize(requested))] + strata:
                    declarations = parsed[sid].get(label, [])
                    values = {normalize(item['value']) for item in declarations}
                    if (len(values) == 1
                            and not {normalized_source_label(value) for value in values} & UNRESOLVED_SOURCE_VALUES
                            and normalized_source_label(expected) not in UNRESOLVED_SOURCE_VALUES
                            and expected not in values):
                        contradictions.append({'sample_id': sid, 'assigned_arm': arm,
                            'source_label': label, 'requested_value': expected,
                            'source_values': sorted(values), 'evidence': [
                                {'sample_id': sid, 'field': item['field'], 'quote': item['quote']}
                                for item in declarations]})
        if contradictions:
            result.update(status='FAIL', explicit_source_contradictions=contradictions)
            result['failures'].append('Assigned samples contradict explicitly named source arm or stratum values')
            return result
    if target_unknown:
        result['unresolved'].append('Requested arm or stratum value is an unresolved source declaration')
        return result
    labels = [explicit_label] if explicit_label is not None else sorted({label for fields in parsed.values() for label in fields})
    uncertain_labels = []
    for label in labels:
        selected = {'case': [], 'control': []}
        quotes = []
        ambiguous = False
        recognizable = label == explicit_label
        for sid, fields in parsed.items():
            declarations = fields.get(label, [])
            values = {normalize(item['value']) for item in declarations}
            terms = values | {normalize(item['quote']) for item in declarations}
            terms |= {normalize(label + ': ' + item['value']) for item in declarations}
            matched = [arm for arm, term in (('case', case_label), ('control', control_label)) if term in terms]
            recognizable |= bool(matched)
            source_strata = [{normalize(item['value']) for item in fields.get(key, [])}
                             for key, _ in strata]
            def known(declarations):
                return len(declarations) == 1 and not {
                    normalized_source_label(value) for value in declarations} & UNRESOLVED_SOURCE_VALUES
            # A known nonmatching stratum proves this record ineligible even
            # when another declaration is unknown. Otherwise unknown eligibility
            # cannot be silently converted into absence from a category.
            if any(known(found) and expected not in found
                   for (_, expected), found in zip(strata, source_strata)):
                continue
            if not known(values):
                ambiguous = True
            elif len(matched) > 1:
                ambiguous = True
            elif matched:
                if not all(known(found) for found in source_strata):
                    ambiguous = True
                    continue
                selected[matched[0]].append(sid)
                quotes.extend({'sample_id': sid, 'field': item['field'], 'quote': item['quote']}
                              for item in declarations)
                quotes.extend({'sample_id': sid, 'field': item['field'], 'quote': item['quote']}
                              for key, _ in strata for item in fields[key])
        if label == explicit_label:
            omissions = [{'sample_id': sid, 'source_arm': arm,
                          'evidence': [item for item in quotes if item['sample_id'] == sid]}
                         for arm, supplied in (('case', case_ids), ('control', control_ids))
                         for sid in selected[arm] if sid not in supplied]
            if omissions:
                result.update(status='FAIL', known_membership_disagreements=omissions)
                result['failures'].append('Assigned arms omit or misplace known eligible samples under the explicit source request')
                return result
        if recognizable and ambiguous:
            uncertain_labels.append(label)
        if not ambiguous and selected['case'] and selected['control']:
            result['bindings'].append({'source_label': label, 'case_ids': selected['case'],
                                       'control_ids': selected['control'], 'evidence': quotes})
    if uncertain_labels:
        result['unresolved'].append('Missing, unknown or conflicting declarations prevent complete source membership: ' + ', '.join(uncertain_labels))
        return result
    if not result['bindings']:
        result['unresolved'].append('Requested arm labels have no unambiguous matching source characteristic')
        return result
    memberships = {(tuple(sorted(item['case_ids'])), tuple(sorted(item['control_ids']))) for item in result['bindings']}
    if len(memberships) != 1:
        result['unresolved'].append('Requested labels map to different source groupings')
        return result
    expected_case, expected_control = next(iter(memberships))
    if set(case_ids) != set(expected_case) or set(control_ids) != set(expected_control):
        result['status'] = 'FAIL'
        result['failures'].append('Assigned arms omit, swap or add samples relative to the source-defined request')
    else:
        result['status'] = 'PASS'
    return result


def propose_source_comparisons(packet, case_ids, control_ids, *, diagnostic_profile=SOURCE_DIAGNOSTIC_CURRENT):
    """Propose exact source labels/strata for existing groups, without approval.

    Each proposal must reproduce complete membership when read back by the
    comparison verifier. Shared strata are reduced by deterministic deletion;
    this defines labels, not a causal model or an independently validated repair.
    """
    source_diagnostic_profile(diagnostic_profile)
    if packet is None or not case_ids or not control_ids:
        return []
    members = case_ids + control_ids
    samples = packet.get('samples', [])
    by_id = {row['sample_id']: structured_characteristics(row) for row in samples}
    if len(by_id) != len(samples) or len(set(members)) != len(members) or set(members) - set(by_id):
        return []
    labels = sorted({label for row in by_id.values() for label in row})
    def values(label, ids):
        found = [{item['value'] for item in by_id[sid].get(label, [])} for sid in ids]
        return set().union(*found) if found and all(len(v) == 1 for v in found) else set()
    common = {label: next(iter(v)) for label in labels if len(v := values(label, members)) == 1}
    proposals = []
    for label in labels:
        cases, controls = values(label, case_ids), values(label, control_ids)
        if len(cases) != 1 or len(controls) != 1 or cases == controls:
            continue
        case, control = next(iter(cases)), next(iter(controls))
        base = f'{label}: {case} versus {label}: {control}'
        strata = [(key, value) for key, value in common.items() if key != label]
        def bind(selected):
            target = base + (' within ' + ' and '.join(f'{key}: {value}' for key, value in selected) if selected else '')
            return target, verify_requested_comparison({**packet, 'target': target}, case_ids, control_ids,
                                                        diagnostic_profile=diagnostic_profile)
        target, bound = bind(strata)
        if bound['status'] != 'PASS':
            continue
        for clause in list(strata):
            shorter = [item for item in strata if item != clause]
            trial_target, trial = bind(shorter)
            if trial['status'] == 'PASS':
                strata, target, bound = shorter, trial_target, trial
        proposals.append({'status': 'SOURCE_LABEL_PROPOSAL', 'original_target': packet.get('target'),
                          'proposed_target': target, 'source_binding': bound,
                          'scientific_status': 'CANNOT VERIFY', 'applied': False,
                          'scope': 'Source labels and sample scope only; biological interpretation and repair approval remain separate'})
    return proposals


def source_field_by_arm(packet, case_ids, control_ids, field='platform_id'):
    """Describe a supplied direct source field by assigned arm, without inference.

    Counts refer to assay records. A separated field is conditional on these
    groups and does not prove their meaning, a field effect or unit independence.
    """
    result = {'status': 'CANNOT VERIFY', 'field': field, 'pattern': 'UNRESOLVED',
              'source_metadata_sha256': packet.get('source_sha256') if packet else None,
              'counts': {}, 'evidence': [], 'unresolved_ids': [],
              'scientific_status': 'CANNOT VERIFY', 'biological_independence_status': 'CANNOT VERIFY',
              'reference_judgment': False, 'counts_are': 'Source assay records; not independent biological units'}
    if packet is None:
        return result
    samples = packet.get('samples', [])
    by_id = {sample['sample_id']: sample for sample in samples}
    members = case_ids + control_ids
    if len(by_id) != len(samples) or len(set(members)) != len(members) or set(members) - set(by_id):
        result.update(status='FAIL', pattern='INVALID_MEMBERSHIP')
        return result
    if not case_ids or not control_ids:
        return result
    for arm, ids in (('case', case_ids), ('control', control_ids)):
        counts = Counter()
        for sid in ids:
            value = by_id[sid].get(field)
            if not isinstance(value, str) or normalized_source_label(value) in UNRESOLVED_SOURCE_VALUES:
                result['unresolved_ids'].append(sid)
                continue
            counts[value.strip()] += 1
            result['evidence'].append({'sample_id': sid, 'assigned_arm': arm, 'field': field, 'quote': value})
        result['counts'][arm] = dict(sorted(counts.items()))
    if result['unresolved_ids']:
        return result
    case, control = result['counts']['case'], result['counts']['control']
    result['status'] = 'PASS'
    if len(case) == len(control) == 1:
        result['pattern'] = ('SAME_CONSTANT' if set(case) == set(control) else 'PERFECT_ARM_SEPARATION')
    else:
        result['pattern'] = 'VARIES_WITHIN_ARM'
    return result


def source_grouping_diagnostics(packet, case_ids, control_ids, *, diagnostic_profile=SOURCE_DIAGNOSTIC_CURRENT):
    """Keep target interpretation separate from conditional source membership.

    A witness describes the supplied groups, not which groups should have been
    selected. Even reversed or out-of-scope groups can have exact source labels.
    Nothing here applies a relabeling or changes a scientific decision.
    """
    legacy = source_diagnostic_profile(diagnostic_profile) == SOURCE_DIAGNOSTIC_LEGACY
    proposals = propose_source_comparisons(packet, case_ids, control_ids, diagnostic_profile=diagnostic_profile)
    return {
        'schema': 'clockbase-source-grouping-diagnostics-v1' if legacy else 'clockbase-source-grouping-diagnostics-v2',
        'source_metadata_sha256': packet.get('source_sha256') if packet else None,
        'scientific_status': 'CANNOT VERIFY', 'reference_judgment': False,
        'target_binding': verify_requested_comparison(packet, case_ids, control_ids, diagnostic_profile=diagnostic_profile),
        **({} if legacy else {'platform_by_arm': source_field_by_arm(packet, case_ids, control_ids)}),
        'source_membership_witnesses': [
            {**proposal, 'case_ids': list(case_ids), 'control_ids': list(control_ids),
             'source_metadata_sha256': packet.get('source_sha256'),
             'reference_judgment': False} for proposal in proposals],
        'scope': 'Exact source-label witnesses conditional on supplied groups; no target interpretation or control-role approval',
        'limitations': ['Witnesses do not establish that the supplied groups answer the requested question',
                        'No biological-unit identity, independence, covariate adjustment or clock applicability is established',
                        'All source-label proposals remain unapplied; target-binding failures and uncertainties are retained']}


def verify_source_design(packet, case_ids, control_ids, *, diagnostic_profile=SOURCE_DIAGNOSTIC_CURRENT):
    """Check the explicit-unit, covariate-matched drug-contrast profile.

    PASS supports this bounded unadjusted comparison under source declarations.
    It does not certify unreported covariates, randomization, causal treatment
    effects or a reference truth label. Missing declarations stay unresolved.
    """
    legacy = source_diagnostic_profile(diagnostic_profile) == SOURCE_DIAGNOSTIC_LEGACY
    scope = 'Unadjusted independent drug-versus-control contrast with explicit biological units and matched source covariates'
    result = {'status': 'CANNOT VERIFY', 'profile': 'explicit-unit-drug-contrast-v2',
              'scope': scope, 'failures': [], 'unresolved': [], 'sample_bindings': {},
              'reference_judgment': False}
    if packet is None:
        result['unresolved'].append('Source packet is unavailable')
        return result
    samples = packet.get('samples', [])
    by_id = {sample['sample_id']: sample for sample in samples}
    members = case_ids + control_ids
    if len(by_id) != len(samples) or len(set(members)) != len(members) or set(members) - set(by_id):
        result.update(status='FAIL', failures=['Missing, duplicated or invented source sample identity'])
        return result
    if not legacy:
        result['platform_by_arm'] = source_field_by_arm(packet, case_ids, control_ids)
    if min(len(case_ids), len(control_ids)) < 2:
        result['unresolved'].append('Two biological units per arm are required to estimate variance')
        return result
    result['target_binding'] = verify_requested_comparison(packet, case_ids, control_ids, diagnostic_profile=diagnostic_profile)
    result['comparison_proposals'] = (propose_source_comparisons(packet, case_ids, control_ids, diagnostic_profile=diagnostic_profile)
                                      if result['target_binding']['status'] != 'PASS' else [])
    result['failures'].extend(result['target_binding']['failures'])
    result['unresolved'].extend(result['target_binding']['unresolved'])
    treatment_binding = (result['target_binding']['status'] == 'PASS' and any(
        item['source_label'] in DESIGN_LABELS['treatment'] for item in result['target_binding']['bindings']))
    if result['target_binding']['status'] == 'PASS' and not treatment_binding:
        result['unresolved'].append('The drug-contrast profile requires a target bound to a source treatment/compound/drug label')
    missing_values = ({'', 'na', 'n a', 'not available', 'unknown', 'missing', 'not reported'}
                      if legacy else UNRESOLVED_SOURCE_VALUES)
    for sid in members:
        source = structured_characteristics(by_id[sid])
        bindings = {}
        for semantic, labels in DESIGN_LABELS.items():
            matches = [item for label, items in source.items() if label in labels for item in items]
            values = {item['value'].strip() if semantic == 'unit' else normalized_source_label(item['value'])
                      for item in matches}
            unknown = {normalized_source_label(value) for value in values} & missing_values
            if semantic == 'unit':
                unknown |= {normalized_source_label(value) for value in values} & {'none', 'not applicable'}
            if len(values) != 1 or unknown:
                result['unresolved'].append(f'{sid}: missing or conflicting explicit {semantic} declaration')
            else:
                bindings[semantic] = {'value': next(iter(values)), 'evidence': matches}
        result['sample_bindings'][sid] = bindings
    bindings = result['sample_bindings']
    if all('unit' in bindings[sid] for sid in members):
        units = [bindings[sid]['unit']['value'] for sid in members]
        duplicates = sorted(unit for unit, count in Counter(units).items() if count > 1)
        if duplicates:
            result['unresolved'].append('Repeated source unit identifiers require pairing/replication or namespace resolution: ' + ', '.join(duplicates))
    drug_contrast = False
    if all('treatment' in bindings[sid] for sid in members):
        case = {bindings[sid]['treatment']['value'] for sid in case_ids}
        control = {bindings[sid]['treatment']['value'] for sid in control_ids}
        result['treatment_values'] = {'case': sorted(case), 'control': sorted(control)}
        if len(case) != 1 or len(control) != 1:
            result['unresolved'].append('Multiple treatment values per arm require a different design profile')
        elif case == control:
            result['unresolved'].append('Treatment labels do not define this contrast; another comparison dimension may apply')
        elif case & CONTROL_VALUES and not control & CONTROL_VALUES:
            result['failures'].append('Case/control direction contradicts the source intervention/control labels')
        elif not control & CONTROL_VALUES:
            result['unresolved'].append('Control condition is outside the declared control vocabulary')
        else:
            drug_contrast = treatment_binding
    # Required covariates are compared as exact distributions, not by a
    # nonsignificant imbalance test. Differing scales/labels need adjudication.
    for covariate in ('age', 'sex', 'genotype', 'tissue', 'batch'):
        if not all(covariate in bindings[sid] for sid in members):
            continue
        a = Counter(bindings[sid][covariate]['value'] for sid in case_ids)
        b = Counter(bindings[sid][covariate]['value'] for sid in control_ids)
        if drug_contrast and len(a) == len(b) == 1 and set(a).isdisjoint(b):
            result['failures'].append('Treatment arm is confounded with source ' + covariate)
        elif {key: Fraction(value, len(case_ids)) for key, value in a.items()} != {
                key: Fraction(value, len(control_ids)) for key, value in b.items()}:
            result['unresolved'].append('Source ' + covariate + ' is not exactly balanced; adjustment/profile review required')
    if legacy:
        for field in ('organism_ch1', 'platform_id'):
            values = [by_id[sid].get(field, '').strip() for sid in members]
            if any(not value for value in values):
                result['unresolved'].append('Missing source ' + field)
            elif len(set(values)) != 1:
                result['unresolved'].append('Mixed source ' + field + ' requires stratification or adjustment')
    else:
        for field in ('organism_ch1', 'platform_id'):
            values = [by_id[sid].get(field, '') for sid in members]
            if ((field == 'platform_id' and result['platform_by_arm']['status'] != 'PASS')
                    or any(not isinstance(value, str) or normalized_source_label(value) in UNRESOLVED_SOURCE_VALUES
                           for value in values)):
                result['unresolved'].append('Missing or unresolved source ' + field)
            elif field == 'platform_id' and result['platform_by_arm']['pattern'] == 'PERFECT_ARM_SEPARATION':
                result['unresolved'].append('Assigned arm and source platform_id have perfect categorical separation; '
                                            'these records cannot identify separate arm and platform terms without additional assumptions')
            elif len({value.strip() for value in values}) != 1:
                result['unresolved'].append('Mixed source ' + field + ' requires stratification or adjustment')
    result['status'] = ('FAIL' if result['failures'] else
                        'CANNOT VERIFY' if result['unresolved'] else 'PASS')
    result['limitations'] = ['Source declarations are not independently established biological truth',
                            'Only listed covariates and explicit unit identities are checked',
                            'No causal or organismal aging claim is authorized']
    return result


def series_design(payload, study_id):
    """Retain GEO design fields, not its abstract or previous analysis outputs."""
    fields = {}
    allowed = {'title', 'overall_design', 'type', 'pubmed_id', 'sample_id'}
    identities = []
    for line in payload.decode('utf-8').splitlines():
        key, sep, value = line.partition(' = ')
        if key == '^SERIES' and sep:
            identities.append(value.strip())
        elif sep and key.startswith('!Series_') and key[8:] in allowed:
            fields.setdefault(key[8:], []).append(value.strip())
    if identities != [study_id]:
        raise ValueError('Study-context GEO identity mismatch')
    return {key: '\n'.join(values) for key, values in fields.items()}


def article_methods(payload, study_id, pmcid, geo_pubmed_ids, *, hide_statistics=False):
    """Extract only Methods sections from linked JATS, without forwarding results.

    Association is checked against the GEO citation or an exact accession in
    the article body. That link does not establish which experimental cohort
    corresponds to each GEO library. Consumers must assess that separately.
    """
    if not re.fullmatch(r'PMC\d+', pmcid):
        raise ValueError('Invalid publication archive identity')
    if b'<!ENTITY' in payload.upper():
        raise ValueError('XML entity declarations are not supported')
    try:
        root = ET.fromstring(payload)
    except ET.ParseError as error:
        raise ValueError('Cannot parse publication JATS: ' + str(error)) from error
    if root.tag != 'article':
        raise ValueError('Publication response is not JATS article XML')
    ids = {}
    for item in root.findall('./front/article-meta/article-id'):
        ids.setdefault(item.get('pub-id-type'), []).append(''.join(item.itertext()).strip())
    if ids.get('pmcid', ids.get('pmc', [])) not in ([pmcid], [pmcid[3:]]):
        raise ValueError('Downloaded publication identity differs from request')
    body = root.find('body')
    if body is None:
        raise ValueError('Publication has no article body')
    linked = bool(set(ids.get('pmid', [])) & set(geo_pubmed_ids))
    if geo_pubmed_ids and not linked:
        raise ValueError('Article does not match a GEO-cited publication; accession mentions may be secondary reuse')
    # An accession in the bibliography alone cannot establish this link.
    body_text = ' '.join(body.itertext())
    if not linked and not re.search(r'(?<![A-Za-z0-9])' + re.escape(study_id) + r'(?![A-Za-z0-9])', body_text):
        raise ValueError('Publication has no verified citation/accession link to this study')
    basis = 'GEO_PUBMED_CITATION' if linked else 'EXACT_ACCESSION_IN_ARTICLE_BODY'
    documents = []
    source_number = 0

    def heading(section):
        title = section.find('title')
        return ' '.join(''.join(title.itertext()).split()) if title is not None else ''

    def walk(section, in_methods=False, blocked=False, statistical=False):
        nonlocal source_number
        title = heading(section)
        normalized = re.sub(r'[^a-z0-9]+', ' ', title.lower()).strip()
        sec_type = section.get('sec-type', '').lower()
        statistical = statistical or bool(re.search(r'\bstatistic\w*\b', normalized)) or bool(
            re.search(r'\bstatistic\w*\b', sec_type))
        # A combined Results/Methods section does not meet the blinded contract.
        blocked = blocked or bool(re.search(r'\b(results?|discussion|conclusions?)\b', normalized)) or sec_type in {
            'results', 'discussion', 'conclusions', 'results-discussion'}
        selected = not blocked and (in_methods or normalized in METHODS_TITLES or sec_type == 'methods')
        if selected:
            paragraphs = []
            for child in section:
                if child.tag in {'p', 'list', 'def-list', 'table-wrap'}:
                    paragraphs.append(' '.join(''.join(child.itertext()).split()))
            text = '\n\n'.join(p for p in paragraphs if p)
            if text:
                source_number += 1
                if not hide_statistics or not statistical:
                    documents.append({'source_id': f'{pmcid}:methods:{source_number}',
                                      'source_kind': 'publication_methods',
                                      'fields': {'title': title, 'text': text},
                                      'link_basis': basis,
                                      'cohort_correspondence': 'CANNOT VERIFY'})
        for child in section.findall('sec'):
            walk(child, selected, blocked, statistical)

    for section in body.findall('sec'):
        walk(section)
    return documents


def context_documents(study_id, sources, directory, *, hide_statistics=False):
    """Reconstruct every quoted field from hashed source files on each load."""
    directory = Path(directory).resolve()
    if not isinstance(sources, list) or not sources:
        raise ValueError('Study context has no source declarations')
    payloads = []
    seen = set()
    for source in sources:
        if set(source) != {'kind', 'id', 'path', 'sha256', 'url'} or source['id'] in seen:
            raise ValueError('Study-context source schema/identity mismatch')
        seen.add(source['id'])
        relative = Path(source['path'])
        path = (directory / relative).resolve()
        if relative.is_absolute() or not path.is_relative_to(directory):
            raise ValueError('Study-context sources must be inside the context directory')
        try:
            payload = path.read_bytes()
        except OSError as error:
            raise ValueError('Cannot read study-context source: ' + str(path)) from error
        if hashlib.sha256(payload).hexdigest() != source['sha256']:
            from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
            raise PipelineIntegrityError('Study-context source hash changed: ' + source['id'])
        if not source['url'].startswith('https://'):
            raise ValueError('Study-context source URL must be HTTPS')
        payloads.append((source, payload))
    series = [(s, p) for s, p in payloads if s['kind'] == 'geo_series']
    if len(series) != 1 or series[0][0]['id'] != study_id:
        raise ValueError('Study context requires one matching GEO series')
    source, payload = series[0]
    fields = series_design(payload, study_id)
    pubmed_ids = fields.pop('pubmed_id', '').splitlines()
    sample_ids = fields.pop('sample_id', '').splitlines()
    if not sample_ids or len(sample_ids) != len(set(sample_ids)):
        raise ValueError('Study-context GEO sample membership is missing/duplicated')
    documents = [{'source_id': 'GEO:' + study_id, 'source_kind': 'geo_series_design',
                  'fields': fields, 'sample_ids': sample_ids,
                  'source_url': source['url'], 'source_sha256': source['sha256']}]
    for source, payload in payloads:
        if source['kind'] == 'geo_series':
            continue
        if source['kind'] != 'publication_jats':
            raise ValueError('Unsupported study-context source type')
        for document in article_methods(payload, study_id, source['id'], pubmed_ids,
                                        hide_statistics=hide_statistics):
            documents.append({**document, 'source_url': source['url'], 'source_sha256': source['sha256']})
    return documents


def read_study_context(path, study_id, *, blinded=False):
    path = Path(path)
    try:
        context = json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise ValueError('Cannot read study context: ' + str(path)) from error
    if context.get('version') != CONTEXT_VERSION or context.get('study_id') != study_id:
        raise ValueError('Study-context version/identity mismatch')
    documents = context_documents(study_id, context['sources'], path.parent)
    if context.get('documents') != documents:
        from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
        raise PipelineIntegrityError('Study-context quotations differ from source extraction')
    if blinded:
        # Validate the preserved unfiltered extraction first. Design filtering
        # must not alter archived quotation identity or rewrite old contexts.
        documents = context_documents(study_id, context['sources'], path.parent, hide_statistics=True)
    return documents, context.get('acquisition', [])


def build_packet(study_id, samples, comparisons):
    """Preserve samples separately; omit predictions and previous model labels."""
    if not isinstance(study_id, str) or not study_id.strip():
        raise ValueError('A study identifier is required')
    selected, seen = [], set()
    for row in samples:
        if row.get('gse_id') != study_id:
            raise ValueError('Sample belongs to a different study')
        gsm = row.get('gsm_id')
        if not isinstance(gsm, str) or not gsm.strip() or gsm in seen:
            raise ValueError('Missing or duplicate sample identifier')
        seen.add(gsm)
        selected.append({k: v for k, v in row.items() if k in SAMPLE_FIELDS or
                         (k.startswith('raw.') and is_sample_method_field(k[4:]))})
    contrasts = []
    for index, row in enumerate(comparisons, 1):
        if row.get('gse_id', study_id) != study_id:
            raise ValueError('Comparison belongs to a different study')
        contrasts.append({'source_row': index, **{
            k: v for k, v in row.items() if k in COMPARISON_FIELDS
        }})
    return {
        'study_id': study_id,
        'scope': 'Study-design review; structured numerical outcomes withheld',
        'evidence_note': 'Raw metadata and harmonized annotations are distinct; '
                         'neither sample titles nor model annotations establish '
                         'independent biological replication.',
        'samples': selected,
        'method_provenance': method_provenance(selected, prefix='raw.'),
        'comparisons': contrasts,
        'scientific_verdict': 'NOT_EVALUATED',
    }


def parse_review(raw):
    """Reject malformed replies; never salvage text into a passing decision."""
    try:
        value = json.loads(raw)
    except (json.JSONDecodeError, TypeError) as exc:
        raise ValueError('Review must be a complete JSON object') from exc
    if not isinstance(value, dict) or set(value) != {'status', 'findings', 'unresolved'}:
        raise ValueError('Review schema mismatch')
    if value['status'] not in {'NEEDS_REVIEW', 'ISSUES_FOUND', 'NO_DESIGN_ISSUE_IDENTIFIED'}:
        raise ValueError('Unsupported review status; model output cannot approve a claim')
    for field in ('findings', 'unresolved'):
        if not isinstance(value[field], list) or any(
            not isinstance(item, str) or not item.strip() for item in value[field]
        ):
            raise ValueError('Findings and unresolved items must be nonempty strings')
    if value['status'] == 'ISSUES_FOUND' and not value['findings']:
        raise ValueError('An issue decision requires findings')
    if value['status'] == 'NEEDS_REVIEW' and not value['unresolved']:
        raise ValueError('An unresolved decision requires unresolved items')
    if value['status'] == 'NO_DESIGN_ISSUE_IDENTIFIED' and (
        value['findings'] or value['unresolved']
    ):
        raise ValueError('A no-issue decision cannot contain issues or unresolved items')
    return value


def read_source_csv(path):
    """Read without pandas' silent duplicate-header renaming or row truncation."""
    try:
        payload = Path(path).read_bytes()
        if payload.startswith(b'version https://git-lfs.github.com/spec/v1'):
            raise ValueError(f'LFS pointer has no source payload: {path}')
        records = list(csv.reader(io.StringIO(payload.decode('utf-8-sig'))))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f'Cannot read source CSV {path}: {exc}') from exc
    if len(records) < 2 or not records[0]:
        raise ValueError(f'Empty or header-only source CSV: {path}')
    header = records[0]
    if any(not column.strip() for column in header):
        raise ValueError(f'Blank column name in {path}; supply an explicit schema')
    if any(len(row) != len(header) for row in records[1:]):
        raise ValueError(f'Ragged CSV records in {path}')
    repeated = {key for key in header if header.count(key) > 1}
    for key in repeated:
        if key not in {'sample_id', 'gsm_id', 'GSM'}:
            raise ValueError(f'Ambiguous duplicate column {key!r} in {path}')
        indices = [i for i, column in enumerate(header) if column == key]
        if any(len({row[i] for i in indices}) != 1 for row in records[1:]):
            raise ValueError(f'Conflicting duplicate sample columns in {path}')
    return [dict(zip(header, row)) for row in records[1:]], sorted(repeated)


def metadata_packet(path, study_id):
    """Adapt one complete metadata CSV; do not read results or comparisons."""
    rows, collapsed = read_source_csv(path)
    id_columns = [key for key in ('gsm_id', 'GSM', 'sample_id', 'geo_accession')
                  if key in rows[0]]
    if not id_columns:
        raise ValueError('No explicit sample identifier column')
    raw_fields = {'age', 'sex', 'gender', 'batch', 'subject', 'donor', 'patient',
                  'pair', 'timepoint', 'time', 'tissue', 'strain', 'genotype',
                  'treatment', 'condition', 'replicate', 'strain background',
                  'strain/background', 'rna-seq run', 'sequencing run'}
    samples = []
    for row in rows:
        if len({row[key] for key in id_columns}) != 1:
            raise ValueError('Conflicting sample identifier aliases')
        if row.get('gse_id', study_id) != study_id:
            raise ValueError('Mixed-study metadata is not accepted')
        selected = {key: value for key, value in row.items() if key in SAMPLE_FIELDS or
                    (key.startswith('raw.') and is_sample_method_field(key[4:]))}
        selected.update(gse_id=study_id, gsm_id=row[id_columns[0]])
        for native, target in [('title', 'raw.title'),
                               ('source_name_ch1', 'raw.source_name_ch1'),
                               ('organism_ch1', 'raw.organism_ch1'),
                               ('molecule_ch1', 'raw.molecule_ch1'),
                               ('platform_id', 'raw.platform_id')]:
            if native in row:
                selected[target] = row[native]
        for key, value in row.items():
            if is_sample_method_field(key):
                target = 'raw.' + key
                if target in selected and selected[target] != value:
                    raise ValueError('Conflicting native and raw method field: ' + key)
                selected[target] = value
        extra = {key: value for key, value in row.items()
                 if key.lower() in raw_fields or key.startswith('characteristics_ch1')}
        if extra:
            selected['raw_characteristics'] = json.dumps(
                {'original': row.get('raw_characteristics', ''), 'fields': extra},
                sort_keys=True)
        samples.append(selected)
    packet = build_packet(study_id, samples, [])
    packet['adapter'] = {'collapsed_identical_id_columns': collapsed,
                         'scope': 'Field-level withholding only; source text and '
                                  'prior harmonized annotations require review.'}
    return packet


def replay_saved_outputs(samples_path, contract_path):
    """Bounded development fixture replay, never scientific acceptance."""
    import numpy
    import scipy
    from scipy import stats

    rows, collapsed = read_source_csv(samples_path)
    try:
        contract = json.loads(Path(contract_path).read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f'Cannot read replay contract {contract_path}: {exc}') from exc
    if contract.get('test') != 'two-sided Welch independent-samples t-test':
        raise ValueError('Unsupported test; paired or adjusted designs need a separate implementation')
    if contract.get('multiplicity_adjustment') != 'none':
        raise ValueError('This replay cannot reconstruct an accepted FDR family')
    expected_effect = ('case mean minus control mean divided by pooled within-group '
                       'sample standard deviation, ddof=1')
    if contract.get('effect_definition') != expected_effect:
        raise ValueError('Unsupported effect definition')
    # This also rejects duplicate IDs and mixed-study records.
    build_packet(contract['study_id'], rows, [])
    by_id = {row['gsm_id']: row for row in rows}
    groups = [contract[key] for key in ('case_gsm_ids', 'control_gsm_ids',
                                      'background_gsm_ids_not_in_drug_test')]
    if any(not isinstance(group, list) or any(not isinstance(x, str) or not x for x in group)
           for group in groups):
        raise ValueError('Each group must explicitly list sample identifiers')
    members = [sample for group in groups for sample in group]
    if len(members) != len(set(members)) or set(members) != set(by_id):
        raise ValueError('Groups must partition the complete input sample set exactly')
    if len(rows) != contract['sample_count_in_fixture']:
        raise ValueError('Sample count differs from the declared contract')
    if min(len(groups[0]), len(groups[1])) < 2:
        raise ValueError('Singleton groups cannot support this variance estimate')
    rtol = float(contract['relative_comparison_tolerance'])
    atol = float(contract['absolute_comparison_tolerance'])
    if not all(math.isfinite(x) and x >= 0 for x in (rtol, atol)):
        raise ValueError('Invalid comparison tolerances')
    references = contract['expected_from_saved_analysis']
    if not isinstance(references, dict) or not references:
        raise ValueError('At least one explicit clock reference is required')
    results = []
    for clock, expected in references.items():
        # Exact column names only; clock aliases and duplicate version labels are not guessed.
        case = [float(by_id[sample][clock]) for sample in groups[0]]
        control = [float(by_id[sample][clock]) for sample in groups[1]]
        if not all(math.isfinite(x) for x in case + control):
            raise ValueError(f'Nonfinite input for {clock}')
        pooled = math.sqrt(((len(case)-1)*statistics.variance(case) +
                            (len(control)-1)*statistics.variance(control)) /
                           (len(case)+len(control)-2))
        if pooled == 0:
            raise ValueError(f'Undefined standardized effect for {clock}')
        p = float(stats.ttest_ind(case, control, equal_var=False).pvalue)
        effect = (statistics.mean(case)-statistics.mean(control))/pooled
        targets = (float(expected['p_value']), float(expected['effect_size']))
        if not all(math.isfinite(x) for x in (p, effect, *targets)):
            raise ValueError(f'Nonfinite calculation or reference for {clock}')
        if not 0 <= targets[0] <= 1:
            raise ValueError(f'Invalid reference P value for {clock}')
        match = math.isclose(p, targets[0], rel_tol=rtol, abs_tol=atol) and \
            math.isclose(effect, targets[1], rel_tol=rtol, abs_tol=atol)
        results.append({'clock': clock, 'p_value': p, 'effect_size': effect,
                        'n_case': len(case), 'n_control': len(control),
                        'saved_output_check': 'PASS' if match else 'FAIL'})
    return {'study_id': contract['study_id'], 'stage': 'saved_output_arithmetic',
            'evaluation_role': 'development_only', 'results': results,
            'status': 'PASS' if all(r['saved_output_check'] == 'PASS' for r in results) else 'FAIL',
            'scientific_status': 'CANNOT VERIFY',
            'unresolved': ['Biological units and study-design adjudication',
                           'Upstream clock input and implementation validity',
                           'Independent reference judgment and accepted testing family'],
            'collapsed_identical_id_columns': collapsed,
            'versions': {'numpy': numpy.__version__, 'scipy': scipy.__version__},
            'model_calls': 0, 'paid_api_cost_usd': 0}


def run_local(mode, inputs, output, *, study_id=None, resume=False):
    """Preserve inputs and publish a complete run atomically; refuse stale resume."""
    output = Path(output)
    paths = {name: Path(path).resolve() for name, path in inputs.items()}
    paths['implementation.py'] = Path(__file__).resolve()
    try:
        payloads = {name: path.read_bytes() for name, path in paths.items()}
        config = {'mode': mode, 'study_id': study_id, 'python': platform.python_version(),
                  'sources': {name: str(path) for name, path in paths.items()}}
        if mode == 'replay':
            import numpy
            import scipy
            config['numpy'] = numpy.__version__
            config['scipy'] = scipy.__version__
        if output.exists():
            if not resume:
                raise ValueError(f'Output exists; use --resume to verify unchanged inputs: {output}')
            manifest = json.loads((output/'manifest.json').read_text())
            if manifest['config'] != config or any(
                (output/'inputs'/name).read_bytes() != value for name, value in payloads.items()
            ):
                raise ValueError('Resume rejected: input, implementation or environment changed')
            result = json.loads((output/'result.json').read_text())
            if result != manifest['result']:
                raise ValueError('Resume rejected: saved result differs from completion manifest')
            return result
        if mode == 'prepare':
            result = metadata_packet(paths['metadata.csv'], study_id)
        elif mode == 'replay':
            result = replay_saved_outputs(paths['samples.csv'], paths['contract.json'])
        else:
            raise ValueError('Unsupported local run mode')
        if any(path.read_bytes() != payloads[name] for name, path in paths.items()):
            raise ValueError('Source changed during execution; no completed output written')
        manifest = {'config': config, 'result': result,
                    'completed_at_utc': datetime.now(timezone.utc).isoformat(),
                    'inputs': {name: {'bytes': len(payloads[name]),
                                     'mtime_ns': path.stat().st_mtime_ns}
                               for name, path in paths.items()}}
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=output.name+'.partial-', dir=output.parent))
        try:
            (staging/'inputs').mkdir()
            for name, payload in payloads.items():
                (staging/'inputs'/name).write_bytes(payload)
            for name, value in [('result.json', result), ('manifest.json', manifest)]:
                (staging/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
            # Existing completed output directories are nonempty and cannot be replaced.
            if output.exists():
                raise ValueError('Another run created the output directory')
            os.rename(staging, output)
        except Exception:
            shutil.rmtree(staging)
            raise
        return result
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f'Local {mode} failed for {output}: {exc}') from exc


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='mode', required=True)
    prepare = commands.add_parser('prepare', help='Build an outcome-withheld metadata packet')
    prepare.add_argument('--metadata', required=True)
    prepare.add_argument('--study', required=True)
    replay = commands.add_parser('replay', help='Replay an explicit development fixture')
    replay.add_argument('--samples', required=True)
    replay.add_argument('--contract', required=True)
    for command in (prepare, replay):
        command.add_argument('--output', required=True)
        command.add_argument('--resume', action='store_true')
    args = parser.parse_args(argv)
    inputs = ({'metadata.csv': args.metadata} if args.mode == 'prepare' else
              {'samples.csv': args.samples, 'contract.json': args.contract})
    try:
        result = run_local(args.mode, inputs, args.output,
                           study_id=getattr(args, 'study', None), resume=args.resume)
        print(json.dumps({'output': args.output,
                          'status': result.get('status', 'PACKET_PREPARED'),
                          'scientific_status': result.get('scientific_status', 'CANNOT VERIFY')}))
        return 1 if result.get('status') == 'FAIL' else 0
    except (ValueError, ImportError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
