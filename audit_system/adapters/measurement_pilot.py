"""Translate a measurement-pilot report into immutable audit records."""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

from ..hashing import ArtifactIdentityError, FileIdentity, inspect_file
from ..models import (
    ActorRecord,
    ActorRole,
    ActorType,
    ArtifactRecord,
    AuditRecord,
    CheckRecord,
    CheckStatus,
    EvidenceClassification,
    IdentityAssurance,
)


BLOOD_CLOCK_PROFILES = {
    'Hannum': {'age_years': (19, 101), 'publication': 'https://doi.org/10.1016/j.molcel.2012.10.016'},
    'PhenoAge': {'age_years': (21, 100), 'publication': 'https://doi.org/10.18632/aging.101414'},
}


# Public model bytes and metadata bridge were independently checked against
# tAge 0dba58fba356fecfbbb7c6f0cb27ced59ee6f99f and Zenodo 18763485.
MOUSE_RNA_MODELS = {
    'Hazard_Mouse_All_All_scaleddiff_EN': 'c201a0f8a46f4d5ddcb09892ff3193450099a46088ec724c970abf86a50a7551',
    'Hazard_Rodents_All_All_scaleddiff_EN': '60aee54a2cf54acabe416f7bb0dabc10318e73fb51a951505a328bcf2031e2dc',
}
MOUSE_RNA_ASSETS = {
    'gene_mapping': '4b9a38152a64ee9cdcfba69074ca73d1ebc285f04e4d1b4afbe18fe738faf313',
    'feature_registry': 'ef3f9e30fe4869a40b2c4915b8e1710e94933191b1e553a6deb0c199ca22af9e',
}


def _rna_sample_field(key):
    return bool(key in {'title', 'source_name_ch1', 'organism_ch1', 'molecule_ch1', 'platform_id'}
        or re.fullmatch(r'(?:characteristics_ch1|(?:treatment|growth|extract)_protocol_ch1|library_(?:strategy|source|selection))(?:\.\d+)?', key))


def _rna_primary_samples(path):
    """Parse primary fields independently of the acquisition/packet code."""
    samples, current = {}, None
    for line in Path(path).read_text().splitlines():
        key, sep, value = line.partition(' = ')
        if key == '^SAMPLE' and sep:
            sid = value.strip()
            if not re.fullmatch(r'GSM\d+', sid) or sid in samples:
                raise ValueError('Primary RNA sample identity is malformed or duplicated')
            current = samples[sid] = {}
        elif current is not None and sep and key.startswith('!Sample_'):
            field = key[len('!Sample_'):]
            if not _rna_sample_field(field): continue
            suffix = 0
            while (field if suffix == 0 else field + '.' + str(suffix)) in current: suffix += 1
            current[field if suffix == 0 else field + '.' + str(suffix)] = value.strip()
    if not samples: raise ValueError('Primary RNA source contains no samples')
    return samples


def _rna_domain(samples):
    """Positive, bounded bulk-adipose evidence; absence is not an error label."""
    from ..source_packet import structured_characteristics, normalized_source_label
    unresolved, witnesses = [], []
    for sid, sample in samples.items():
        def values(field):
            return [value for key, value in sample.items()
                    if key == field or re.fullmatch(re.escape(field) + r'\.\d+', key)]
        expected = {'organism_ch1': 'Mus musculus', 'molecule_ch1': 'total RNA',
                    'library_strategy': 'RNA-Seq', 'library_source': 'transcriptomic'}
        for field, value in expected.items():
            if set(values(field)) != {value}:
                unresolved.append(sid + ': source ' + field + ' is outside or absent from the declared profile')
        fields = structured_characteristics(sample)
        tissue = {normalized_source_label(row['value']) for key in ('tissue', 'tissue type') for row in fields.get(key, [])}
        if tissue != {'adipose tissue'}: unresolved.append(sid + ': one explicit adipose-tissue stratum is required')
        protocols = [(key, value) for key, value in sample.items()
                     if re.fullmatch(r'extract_protocol_ch1(?:\.\d+)?', key)]
        # These sentence forms are deliberately bounded. Other explicit specimen
        # protocols need a reviewed profile, not inference from a tissue word.
        specimen = [(k, v) for k, v in protocols if re.search(
            r'\bRNA was (?:isolated|extracted)\b[^.\n]*\bfrom homogenized adipose tissue\b', v, re.I)]
        library = [(k, v) for k, v in protocols if re.search(
            r'\bLibraries were prepared using\b[^.\n]*\bmRNA Library Prep(?:aration)? Kit\b', v, re.I)]
        exclusions = re.search(r'\b(?:single[ -](?:cell|nucleus)|metacells?|cultured cells?|cell lines?|microarrays?)\b',
                               ' '.join(sample.values()), re.I)
        negated = any(re.search(r'\b(?:not|never|without)\b', v, re.I) for _, v in specimen + library)
        if not specimen or not library or exclusions or negated:
            unresolved.append(sid + ': positive homogenized-bulk-adipose/mRNA-library evidence is unavailable or conflicting')
        for field, quote in specimen + library:
            witnesses.append({'sample_id': sid, 'source_field': field, 'quote': quote})
    return unresolved, witnesses


def _rna_matrix(path):
    import csv
    import pandas as pd
    with Path(path).open() as handle: header = next(csv.reader(handle), None)
    if not header or len(header) < 2 or len(set(header[1:])) != len(header[1:]) or any(not item for item in header[1:]):
        raise ValueError('RNA matrix has absent or duplicate column identities')
    frame = pd.read_csv(path, index_col=0)
    if frame.empty or frame.index.has_duplicates or frame.index.isna().any():
        raise ValueError('RNA matrix has absent or duplicate row identities')
    frame.index = frame.index.astype(str)
    if frame.index.has_duplicates: raise ValueError('RNA row identity conversion is ambiguous')
    return frame


def _rna_source_calculation(raw, mapping, registry):
    """Reconstruct the public recipe without importing producer helpers or R."""
    import numpy as np
    import pandas as pd
    counts = raw.to_numpy(dtype=float)
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError('Source RNA counts must be finite and nonnegative')
    keep = (counts >= 10).sum(axis=1) >= counts.shape[1] * .2
    accumulated = {}
    for gene, row in zip(raw.index[keep], counts[keep]):
        target = mapping.get(gene)
        if target is not None:
            accumulated.setdefault(target, np.zeros(counts.shape[1], dtype=float))[:] += row
    if len(accumulated) < 2: raise ValueError('Mapped RNA universe cannot support sample SD scaling')
    mapped = pd.DataFrame([accumulated[g] for g in sorted(accumulated)], index=sorted(accumulated), columns=raw.columns)
    array = mapped.to_numpy()
    sizes = array.sum(axis=0)
    positive = (array > 0).all(axis=1)
    if not positive.any() or not np.isfinite(sizes).all() or (sizes <= 0).any():
        raise ValueError('RLE reference or library size is undefined')
    geometric = np.exp(np.log(array[positive]).mean(axis=1))
    relative = np.median(array[positive] / geometric[:, None], axis=0) / sizes
    factors = relative / np.exp(np.log(relative).mean())
    multipliers = 1e7 / (sizes * factors)
    logged = np.log10(array * multipliers[None, :] + 1)
    sd = logged.std(axis=0, ddof=1)
    if not np.isfinite(logged).all() or not np.isfinite(sd).all() or (sd <= 0).any():
        raise ValueError('RNA sample scaling is undefined')
    standardized = (logged - logged.mean(axis=0)) / sd
    scaled = pd.DataFrame(standardized, index=mapped.index, columns=mapped.columns).reindex(registry).T
    centered = scaled - scaled.median(axis=0)
    normalization = pd.DataFrame({'library_size': sizes, 'norm_factor': factors, 'count_multiplier': multipliers}, index=raw.columns)
    return mapped, scaled, centered, normalization


def _verify_mouse_rna(packet, clocks, clock, selected_ids, check, result):
    import csv
    import joblib
    import numpy as np
    import pandas as pd
    import sklearn
    result.update(profile='published-mouse-bulk-adipose-rle-scaleddiff-v1',
        profile_source='https://doi.org/10.1038/s41586-026-10542-3',
        application_code_commit='0dba58fba356fecfbbb7c6f0cb27ced59ee6f99f',
        scope='Cohort-centered log10 expected-mortality model score in source-declared mouse bulk adipose; no age, survival, efficacy or aging-reversal inference',
        biological_status='CANNOT VERIFY', coverage_threshold=None,
        numerical_checks=[], source_witnesses=[])
    if check.get('normalization_profile') != 'published_mouse_rle':
        result['unresolved'].append('This RNA measurement profile requires the explicit published_mouse_rle recipe')
        return result
    required = {'raw_counts', 'gene_mapping', 'feature_registry', 'source_packet', 'sample_source', 'metadata_source',
                'preprocessing', 'transformed_input', 'model', 'feature_coverage',
                'preprocessing_python', 'preprocessing_r', 'ortholog_r'}
    declarations = check.get('measurement_sources', {})
    if not isinstance(declarations, dict): declarations = {}
    missing = sorted(required - set(declarations))
    implementation_labels = {'producer', 'prediction_python', 'parameter_verifier'}
    implementations = clocks.get('prediction_implementation_sources', {})
    if (missing or not packet or not selected_ids or not clocks.get('prediction_sha256')
            or not isinstance(implementations, dict) or set(implementations) != implementation_labels
            or not clocks.get('implementation_sha256')):
        result['unresolved'].append('RNA source/contrast/execution provenance is incomplete: ' + ', '.join(missing))
        return result
    if sklearn.__version__ != '1.3.2':
        result['unresolved'].append('RNA fitted-model inspection requires the declared sklearn 1.3.2 runtime')
        return result

    def register(label, declaration, *, implementation=False):
        if not isinstance(declaration, dict) or not declaration.get('path') or not declaration.get('sha256'):
            raise KeyError('Source declaration absent or incomplete: ' + label)
        if implementation and (not declaration.get('executed_path')
                or Path(declaration['path']).resolve() == Path(declaration['executed_path']).resolve()):
            result['unresolved'].append('Immutable implementation snapshot is unavailable: ' + label)
            return Path(declaration['path'])
        if not Path(declaration['path']).is_file():
            raise FileNotFoundError(2, 'Required RNA artifact is unavailable', declaration['path'])
        try:
            identity = inspect_file(Path(declaration['path']))
        except ArtifactIdentityError as error:
            raise OSError('RNA artifact cannot be read: ' + declaration['path']) from error
        if identity.sha256 != declaration['sha256']:
            raise ArtifactIdentityError('RNA bound artifact identity differs: ' + label)
        result['evidence_files'].append({'kind': label, 'path': str(identity.path), 'sha256': identity.sha256})
        return identity.path

    def compare(label, observed, expected, *, exact=False):
        if list(observed.index) != list(expected.index) or list(observed.columns) != list(expected.columns):
            raise ValueError('RNA axes differ: ' + label)
        a, b = observed.to_numpy(dtype=float), expected.to_numpy(dtype=float)
        if not np.array_equal(np.isnan(a), np.isnan(b)) or np.isinf(a).any():
            raise ValueError('RNA missing-feature mask differs: ' + label)
        okay = np.array_equal(a, b, equal_nan=True) if exact else np.allclose(a, b, equal_nan=True, rtol=1e-10, atol=1e-12)
        if not okay: raise ValueError('RNA numerical values differ: ' + label)
        finite = np.isfinite(b)
        result['numerical_checks'].append({'step': label, 'status': 'PASS',
            'max_abs_difference': float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else None})

    try:
        files = {label: register(label, declarations[label], implementation=label in
                 {'preprocessing_python', 'preprocessing_r', 'ortholog_r'}) for label in sorted(required)}
        if any(not isinstance(item, dict) or not isinstance(item.get('path'), str) or not item['path']
               or not isinstance(item.get('sha256'), str) or not item['sha256'] for item in implementations.values()):
            raise KeyError('A labeled prediction implementation declaration is unavailable')
        declared_implementations = {item['path']: item['sha256'] for item in implementations.values()}
        if len(declared_implementations) != len(implementation_labels) or declared_implementations != clocks['implementation_sha256']:
            raise KeyError('Complete prediction implementation bindings are unavailable')
        for label, declaration in implementations.items():
            register('prediction_implementation:' + label, declaration, implementation=True)
        prediction_path = register('clock_predictions', {'path': clocks['path'], 'sha256': clocks['prediction_sha256']})
        for label, digest in {**MOUSE_RNA_ASSETS, 'model': MOUSE_RNA_MODELS[clock.removesuffix('_tAge')]}.items():
            if declarations[label]['sha256'] != digest:
                result['unresolved'].append('RNA asset is outside the pinned public application profile: ' + label)
        if result['unresolved']: return result
        bound = json.loads(files['source_packet'].read_text())
        for field in ('study_id', 'samples', 'source_sha256', 'sample_source', 'metadata_source'):
            if packet.get(field) != bound.get(field): raise ValueError('RNA source packet binding differs: ' + field)
        # Targets, summaries and design proposals are deliberately not compared.
        for label in ('sample_source', 'metadata_source'):
            if packet[label] != declarations[label]: raise ValueError('RNA primary source declaration differs: ' + label)
        if packet['source_sha256'] != declarations['metadata_source']['sha256']:
            raise ValueError('RNA metadata file identity differs from source packet')
        source = _rna_primary_samples(files['sample_source'])
        ids = [sample['sample_id'] for sample in packet['samples']]
        if len(ids) != len(set(ids)) or set(ids) != set(source) or len(selected_ids) != len(set(selected_ids)) or set(selected_ids) - set(ids):
            raise ValueError('RNA source/contrast sample membership is absent or duplicated')
        with files['metadata_source'].open() as handle: metadata = list(csv.DictReader(handle))
        meta = {}
        for row in metadata:
            aliases = [row[k] for k in ('GSM', 'gsm_id', 'geo_accession') if row.get(k)]
            if not aliases or len(set(aliases)) != 1 or aliases[0] in meta or row.get('gse_id', packet['study_id']) != packet['study_id']:
                raise ValueError('RNA metadata sample identity is ambiguous')
            meta[aliases[0]] = {k: v for k, v in row.items() if _rna_sample_field(k)}
        if set(meta) != set(source): raise ValueError('RNA metadata/primary sample roster differs')
        for sample in packet['samples']:
            sid = sample['sample_id']
            visible = {k: v for k, v in sample.items() if _rna_sample_field(k)}
            if visible != meta[sid] or any(v != source[sid].get(k, '') for k, v in visible.items()):
                raise ValueError('RNA source fields differ from primary metadata: ' + sid)
        domain_unresolved, result['source_witnesses'] = _rna_domain(source)
        result['unresolved'].extend(domain_unresolved)
        raw = _rna_matrix(files['raw_counts'])
        sample_ids = list(raw.columns)
        if set(sample_ids) != set(source): raise ValueError('RNA raw/source sample membership differs')
        if result['unresolved']: return result
        preprocessing = json.loads(files['preprocessing'].read_text())
        if (preprocessing['normalization_profile'] != 'published_mouse_rle' or preprocessing['filter_order'] != 'raw_then_mapping'
                or preprocessing['filter_count_minimum'] != 10 or preprocessing['filter_sample_fraction'] != .2
                or not preprocessing['RLE_applied'] or preprocessing['YuGene_applied']):
            raise ValueError('RNA declared operation order differs from the selected profile')
        stratum = preprocessing['source_stratum']
        if (stratum['sample_ids'] != sample_ids or stratum['reference_ids'] != sample_ids
                or stratum['reference_rule'] != 'Per-gene median across every source sample in this stratum'
                or stratum['species'] != 'Mus musculus' or stratum['molecule'] != 'total RNA'
                or stratum['tissue'] != 'adipose tissue' or stratum['matched_controls_verified'] is not False):
            raise ValueError('RNA declared processing/reference pool differs from the bounded profile')
        preprocessing_sources = required - {'preprocessing', 'transformed_input', 'model', 'feature_coverage'}
        if (not isinstance(preprocessing.get('measurement_sources'), dict)
                or set(preprocessing['measurement_sources']) != preprocessing_sources):
            raise KeyError('Complete preprocessing source bindings are unavailable')
        for label, declaration in preprocessing['measurement_sources'].items():
            if declaration != declarations.get(label): raise ValueError('RNA preprocessing source binding differs: ' + label)
        table = pd.read_csv(files['gene_mapping'], index_col=0).dropna(subset=['Entrez'])
        mapping = {str(k): str(int(v)) for k, v in table['Entrez'].items()}
        registry = files['feature_registry'].read_text().splitlines()
        mapped, scaled, centered, factors = _rna_source_calculation(raw, mapping, registry)
        if (preprocessing['sample_count'] != len(sample_ids) or clocks['sample_count'] != len(sample_ids)
                or preprocessing['mapped_filtered_genes'] != len(mapped)
                or preprocessing['clock_feature_columns'] != len(registry)):
            raise ValueError('RNA declared source/feature dimensions differ from reconstruction')
        expected_outputs = {'mapped_filtered_counts.csv': mapped, 'scaled.csv': scaled,
                            'scaled_diff.csv': centered, 'rle_factors.csv': factors}
        for name, expected in expected_outputs.items():
            path = register('preprocessing_output', {'path': str(files['preprocessing'].parent / name),
                                                    'sha256': preprocessing['output_sha256'][name]})
            compare(name, _rna_matrix(path), expected, exact=name == 'mapped_filtered_counts.csv')
        if (Path(preprocessing['scaled_diff']).resolve() != files['transformed_input'].resolve()
                or declarations['transformed_input']['sha256'] != preprocessing['output_sha256']['scaled_diff.csv']):
            raise ValueError('RNA transformed-input identity differs from the checked preprocessing output')
        model = joblib.load(files['model'])  # Only allowlisted bytes reach deserialization.
        if [type(step).__name__ for _, step in model.steps] != ['SimpleImputer', 'StandardScaler', 'SelectKBest', 'ElasticNet']:
            raise ValueError('RNA fitted estimator sequence differs')
        imputer, scaler, selector, estimator = [step for _, step in model.steps]
        features = [str(value) for value in model.feature_names]
        if features != [str(value) for value in model.feature_names_in_] or len(features) != len(set(features)):
            raise ValueError('RNA fitted feature order is ambiguous')
        if (imputer.strategy != 'median' or imputer.add_indicator or not np.isnan(imputer.missing_values)
                or not np.isfinite(imputer.statistics_).all() or scaler.with_std or not scaler.with_mean):
            raise ValueError('RNA fitted imputation/scaling differs from the allowlisted semantics')
        values = centered[features].to_numpy(dtype=float)
        missing = np.isnan(values)
        if (check['model_sha256'] != declarations['model']['sha256'] or check['required_features'] != len(features)
                or check['observed_features_min'] != int((~missing).sum(axis=1).min())
                or check['score_units'] != 'log10 mortality hazard-ratio score; not years'):
            raise ValueError('RNA model identity, feature summary or native score units differ from reconstruction')
        scores = np.asarray(selector.scores_).copy()
        scores[np.isnan(scores)] = np.finfo(scores.dtype).min
        mask = np.ones(len(scores), dtype=bool) if selector.k == 'all' else np.zeros(len(scores), dtype=bool)
        if selector.k != 'all' and selector.k > 0: mask[np.argsort(scores, kind='mergesort')[-selector.k:]] = True
        selected = np.asarray(features)[mask]
        used = selected[np.asarray(estimator.coef_) != 0]
        expected_coverage = {'feature_ids': features, 'training_retained_ids': features, 'training_dropped_ids': [],
            'selected_ids': selected.tolist(), 'nonzero_coefficient_ids': used.tolist(), 'imputation_strategy': 'median',
            'missing_ids_by_sample': {sid: np.asarray(features)[missing[i]].tolist() for i, sid in enumerate(sample_ids)},
            'imputed_ids_by_sample': {sid: np.asarray(features)[missing[i]].tolist() for i, sid in enumerate(sample_ids)},
            'observed_selected_by_sample': np.isfinite(centered[selected]).sum(axis=1).astype(int).to_dict(),
            'observed_nonzero_by_sample': np.isfinite(centered[used]).sum(axis=1).astype(int).to_dict()}
        if json.loads(files['feature_coverage'].read_text()) != expected_coverage:
            raise ValueError('RNA feature coverage or fitted imputation masks differ')
        wholly_imputed = any(count == 0 for count in expected_coverage['observed_nonzero_by_sample'].values())
        if wholly_imputed:
            result['unresolved'].append('RNA score has no measured nonzero-coefficient model features; wholly imputed measurement is unresolved')
        filled = np.where(missing, imputer.statistics_, values)
        calculated = (filled - scaler.mean_)[:, mask] @ estimator.coef_ + estimator.intercept_
        predictions = _rna_matrix(prediction_path)
        column = clock if clock.endswith('_tAge') else clock + '_tAge'
        if column not in predictions or list(predictions.columns) != clocks.get('clocks'):
            raise ValueError('RNA prediction clock axis differs from the executed receipt/request')
        compare('native fitted-parameter score', predictions[[column]], pd.DataFrame({column: calculated}, index=sample_ids))
        result.update(measured_input_prediction_check='CANNOT VERIFY' if wholly_imputed else 'PASS', raw_to_score_verification='PASS',
            processing_sample_ids=sample_ids, reference_sample_ids=sample_ids, selected_sample_ids=list(selected_ids),
            raw_feature_rows=len(raw), mapped_features=len(mapped), required_features=len(features),
            observed_features_by_sample=dict(zip(sample_ids, (~missing).sum(axis=1).astype(int).tolist())),
            imputed_features_by_sample=dict(zip(sample_ids, missing.sum(axis=1).astype(int).tolist())),
            observed_nonzero_by_sample=expected_coverage['observed_nonzero_by_sample'],
            score_units='log10 mortality hazard-ratio score; not years',
            imputation='Persisted fitted training medians; no cohort refitting or zero padding')
        for evidence in result['evidence_files']:
            if inspect_file(Path(evidence['path'])).sha256 != evidence['sha256']:
                raise ArtifactIdentityError('RNA evidence changed during verification: ' + evidence['kind'])
    except FileNotFoundError as error:
        result['unresolved'].append('Cannot proceed - required data file not found: ' + str(error.filename))
    except KeyError as error:
        result['unresolved'].append('RNA required source or receipt field is unavailable: ' + str(error))
    except OSError as error:
        result['unresolved'].append('RNA required evidence cannot be accessed: ' + str(error))
    except (ValueError, TypeError, IndexError, ArtifactIdentityError) as error:
        result['failures'].append(f'RNA measurement verification failed: {type(error).__name__}: {error}')
    result['status'] = 'FAIL' if result['failures'] else 'CANNOT VERIFY' if result['unresolved'] else 'PASS'
    return result


def verify_clock_measurement(packet, clocks, clock, selected_ids):
    """Verify the declared blood-clock or public mouse RNA numerical profile.

    PASS is limited to a numerical clock-output comparison. Reusing population
    imputation for DNAm and unsupported tissues/platforms/ages require further
    profiles. The RNA route checks a source-bound score recipe, not biological
    transportability. Unsupported domains remain unresolved.
    """
    from ..source_packet import structured_characteristics, normalized_source_label
    result = {'clock': clock, 'status': 'CANNOT VERIFY', 'failures': [], 'unresolved': [],
              'evidence_files': [], 'reference_judgment': False,
              'scope': 'Measured-feature clock output within a declared input domain; no causal or aging-reversal inference'}
    candidates = [r for r in clocks.get('checks', []) if isinstance(r.get('clock'), str)
                  and r['clock'] in {clock, clock.removesuffix('_tAge')}]
    if len(candidates) != 1:
        result['unresolved'].append('Unique executed clock receipt is unavailable')
        return result
    check = candidates[0]
    if check.get('numerical_status') != 'PASS':
        result['unresolved'].append('Clock parameter verification did not pass')
    if clock.removesuffix('_tAge') in MOUSE_RNA_MODELS:
        return _verify_mouse_rna(packet, clocks, clock, selected_ids, check, result)
    if clock not in BLOOD_CLOCK_PROFILES:
        result['unresolved'].append('No source-verified input-domain profile is implemented for this clock')
        return result
    profile = BLOOD_CLOCK_PROFILES[clock]
    result.update(profile='complete-measured-cpgs-450k-whole-blood-v1', profile_source=profile['publication'],
                  declared_age_range_years=list(profile['age_years']))
    if packet is None or not selected_ids:
        result['unresolved'].append('Source sample context/contrast is unavailable')
        return result
    source_samples = {row['sample_id']: row for row in packet['samples']}
    if len(selected_ids) != len(set(selected_ids)) or set(selected_ids) - set(source_samples):
        result.update(status='FAIL', failures=['Measurement comparison has duplicate or unavailable source samples'])
        return result
    for sid in selected_ids:
        sample = source_samples[sid]
        structured = structured_characteristics(sample)
        if sample.get('organism_ch1') != 'Homo sapiens':
            result['unresolved'].append(sid + ': outside the declared human clock domain')
        if normalized_source_label(sample.get('molecule_ch1', '')) != 'genomic dna':
            result['unresolved'].append(sid + ': genomic DNA source is not explicit')
        if sample.get('platform_id') != 'GPL13534':
            result['unresolved'].append(sid + ': outside the declared 450k assay profile')
        tissues = {normalized_source_label(item['value']) for key in ('tissue', 'tissue type')
                   for item in structured.get(key, [])}
        if not tissues:
            tissues = {normalized_source_label(sample.get('source_name_ch1', ''))}
        if not tissues or not tissues <= {'whole blood', 'peripheral whole blood'}:
            result['unresolved'].append(sid + ': whole-blood source is not explicit')
        ages = [item['value'] + (' years' if key in {'age y', 'age yr', 'age yrs', 'age years'}
                and re.fullmatch(r'\d+(?:\.\d+)?', item['value'].strip()) else '')
                for key in ('age', 'age at collection', 'age at harvest', 'age y', 'age yr', 'age yrs', 'age years')
                for item in structured.get(key, [])]
        values = [re.fullmatch(r'(\d+(?:\.\d+)?)\s*(?:years?|yrs?|y)', value.strip(), re.I) for value in ages]
        if not values or any(value is None for value in values):
            result['unresolved'].append(sid + ': age in years is unavailable/ambiguous')
        else:
            numeric = {float(value.group(1)) for value in values}
            if len(numeric) != 1 or not all(profile['age_years'][0] <= age <= profile['age_years'][1] for age in numeric):
                result['unresolved'].append(sid + ': age is outside the declared original-study range')
    declarations = check.get('measurement_sources', {})
    if set(declarations) != {'raw_beta', 'coefficients'}:
        result['unresolved'].append('Raw-beta and coefficient source identities are unavailable')
        return result
    try:
        files = {}
        for label, declaration in declarations.items():
            if not Path(declaration['path']).is_file():
                result['unresolved'].append('Cannot proceed - required data file not found: ' + declaration['path'])
                return result
            identity = inspect_file(Path(declaration['path']))
            if identity.sha256 != declaration['sha256']:
                raise ArtifactIdentityError('Measurement source hash differs: ' + label)
            files[label] = identity.path
            result['evidence_files'].append({'kind': label, 'path': str(identity.path), 'sha256': identity.sha256})
        import numpy as np
        import pandas as pd
        from biolearn.model_gallery import ModelGallery
        from src.pipeline_inputs_2026_09_10 import read_matrix
        raw = read_matrix(files['raw_beta'], beta=True)
        if set(raw.columns) != set(source_samples):
            raise ArtifactIdentityError('Raw beta membership differs from source packet')
        coefficients = pd.read_csv(files['coefficients'], index_col=0).iloc[:, 0]
        reference = ModelGallery().get(clock, imputation_method='none').coefficients.iloc[:, 0]
        if coefficients.index.has_duplicates or not coefficients.index.equals(reference.index) or not np.allclose(
                coefficients.to_numpy(float), reference.to_numpy(float), rtol=1e-12, atol=1e-14):
            raise ArtifactIdentityError('Clock coefficients differ from the pinned model implementation')
        sites = [site for site in coefficients.index if site != 'intercept']
        values = raw.reindex(sites).loc[:, selected_ids]
        observed_by_sample = values.notna().sum(axis=0).astype(int).to_dict()
        result.update(required_features=len(sites), observed_features_by_sample=observed_by_sample)
        if not all(observed_by_sample.values()):
            result['failures'].append('A selected sample has no measured required CpGs')
        elif not all(count == len(sites) for count in observed_by_sample.values()):
            result['unresolved'].append('Imputed required CpGs need a separately justified missing-feature profile')
        else:
            calculated = coefficients.loc[sites].to_numpy(float) @ values.to_numpy(float) + float(coefficients.get('intercept', 0))
            predictions = pd.read_csv(clocks['path']).set_index('sample_id', verify_integrity=True)
            reported = predictions.loc[selected_ids, clock].to_numpy(float)
            if not np.isfinite(reported).all() or not np.allclose(calculated, reported, rtol=1e-10, atol=1e-10):
                result['failures'].append('Clock predictions do not match measured CpGs and published coefficients')
            else:
                result['measured_input_prediction_check'] = 'PASS'
        for evidence in result['evidence_files']:
            if inspect_file(Path(evidence['path'])).sha256 != evidence['sha256']:
                raise ArtifactIdentityError('Measurement source changed during verification')
    except FileNotFoundError as error:
        result['unresolved'].append('Cannot proceed - required data file not found: ' + str(error.filename))
    except (OSError, ValueError, KeyError, ArtifactIdentityError) as error:
        result['failures'].append(f'Measurement source verification failed: {type(error).__name__}: {error}')
    result['status'] = ('FAIL' if result['failures'] else
                        'CANNOT VERIFY' if result['unresolved'] else 'PASS')
    return result


def _timestamp(path: Path) -> str:
    value = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
    return value.isoformat().replace("+00:00", "Z")


def _assert_identity(actual: FileIdentity, expected: FileIdentity) -> None:
    if actual.sha256 != expected.sha256 or actual.byte_count != expected.byte_count:
        raise ArtifactIdentityError(
            f"Artifact identity changed for {actual.path}: expected "
            f"{expected.sha256}/{expected.byte_count}, observed "
            f"{actual.sha256}/{actual.byte_count}"
        )


def _declared_identity(value: dict, label: str) -> tuple[Path, str, int]:
    try:
        return Path(value["path"]), str(value["sha256"]), int(value["size_bytes"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ArtifactIdentityError(f"Invalid declared identity for {label}: {exc}") from exc


def _verify_declared_file(value: dict, label: str) -> FileIdentity:
    path, digest, byte_count = _declared_identity(value, label)
    actual = inspect_file(path)
    if actual.sha256 != digest or actual.byte_count != byte_count:
        raise ArtifactIdentityError(
            f"Declared identity mismatch for {label} at {path}"
        )
    return actual


def adapt_measurement_pilot(
    path: Path,
    release_id: str,
    expected_identity: FileIdentity | None = None,
) -> tuple[AuditRecord, ...]:
    """Read and translate one pilot report without changing source bytes."""

    report_path = Path(path)
    report_identity = inspect_file(report_path)
    if expected_identity is not None:
        _assert_identity(report_identity, expected_identity)
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ArtifactIdentityError(f"Cannot parse pilot report {report_path}: {exc}") from exc

    audit = report.get("audit")
    if not isinstance(audit, dict) or not audit.get("version"):
        raise ArtifactIdentityError(f"Pilot report lacks an audit version: {report_path}")
    source_identity = _verify_declared_file(report["input_source"], "input_source")
    module_identity = _verify_declared_file(
        report["measurement_module"],
        "measurement_module",
    )
    created_at = _timestamp(report_path)
    suffix = report_identity.sha256[:16]
    adapter_id = f"actor:measurement-adapter:{suffix}"
    checker_id = f"actor:measurement-checker:{module_identity.sha256[:16]}"
    report_id = f"artifact:measurement-report:{suffix}"
    source_id = f"artifact:measurement-source:{source_identity.sha256[:16]}"
    module_id = f"artifact:measurement-module:{module_identity.sha256[:16]}"
    predicates = tuple(
        sorted(
            {
                *(str(item.get("name")) for item in audit.get("checks", [])),
                "preflight_status",
                "scientific_status",
                "release_ready",
            }
        )
    )
    records: list[AuditRecord] = [
        ActorRecord(
            record_id=adapter_id,
            release_id=release_id,
            created_at=created_at,
            actor_type=ActorType.PROGRAM,
            role=ActorRole.EVIDENCE_CUSTODIAN,
            organization="ClockBase independent audit",
            model_family=None,
            authority_scope=("artifact_registration",),
            identity_assurance=IdentityAssurance.VERIFIED,
        ),
        ActorRecord(
            record_id=checker_id,
            release_id=release_id,
            created_at=created_at,
            actor_type=ActorType.PROGRAM,
            role=ActorRole.DETERMINISTIC_CHECKER,
            organization="ClockBase measurement preflight",
            model_family=None,
            authority_scope=predicates,
            identity_assurance=IdentityAssurance.VERIFIED,
        ),
        ArtifactRecord(
            record_id=report_id,
            release_id=release_id,
            created_at=created_at,
            uri=str(report_identity.path),
            sha256=report_identity.sha256,
            byte_count=report_identity.byte_count,
            classification=EvidenceClassification.DOCUMENT,
            created_by_actor_id=adapter_id,
            parent_artifact_ids=(source_id, module_id),
        ),
        ArtifactRecord(
            record_id=source_id,
            release_id=release_id,
            created_at=created_at,
            uri=str(source_identity.path),
            sha256=source_identity.sha256,
            byte_count=source_identity.byte_count,
            classification=EvidenceClassification.EMPIRICAL,
            created_by_actor_id=adapter_id,
        ),
        ArtifactRecord(
            record_id=module_id,
            release_id=release_id,
            created_at=created_at,
            uri=str(module_identity.path),
            sha256=module_identity.sha256,
            byte_count=module_identity.byte_count,
            classification=EvidenceClassification.DOCUMENT,
            created_by_actor_id=adapter_id,
        ),
    ]
    for item in audit.get("checks", []):
        name = str(item["name"])
        status = CheckStatus(str(item["status"]))
        records.append(
            CheckRecord(
                record_id=f"check:measurement:{name}:{suffix}",
                release_id=release_id,
                created_at=created_at,
                subject_id=report_id,
                predicate=name,
                status=status,
                checker_actor_id=checker_id,
                evidence_artifact_ids=(report_id, source_id, module_id),
                expected={"message": str(item.get("message", ""))},
                observed={"details": item.get("details", {})},
            )
        )
    records.extend(
        (
            CheckRecord(
                record_id=f"check:measurement:preflight:{suffix}",
                release_id=release_id,
                created_at=created_at,
                subject_id=report_id,
                predicate="preflight_status",
                status=CheckStatus(str(audit["preflight_status"])),
                checker_actor_id=checker_id,
                evidence_artifact_ids=(report_id,),
                expected={"status": "PASS"},
                observed={"raw_status": audit["preflight_status"]},
            ),
            CheckRecord(
                record_id=f"check:measurement:scientific:{suffix}",
                release_id=release_id,
                created_at=created_at,
                subject_id=report_id,
                predicate="scientific_status",
                status=CheckStatus.BLOCKED,
                checker_actor_id=checker_id,
                evidence_artifact_ids=(report_id,),
                expected={"status": "EVALUATED"},
                observed={"raw_status": report.get("scientific_status")},
            ),
            CheckRecord(
                record_id=f"check:measurement:release-ready:{suffix}",
                release_id=release_id,
                created_at=created_at,
                subject_id=report_id,
                predicate="release_ready",
                status=(
                    CheckStatus.NEEDS_REVIEW
                    if report.get("release_ready")
                    else CheckStatus.BLOCKED
                ),
                checker_actor_id=checker_id,
                evidence_artifact_ids=(report_id,),
                expected={"release_ready": True},
                observed={"release_ready": bool(report.get("release_ready"))},
            ),
        )
    )
    return tuple(records)
