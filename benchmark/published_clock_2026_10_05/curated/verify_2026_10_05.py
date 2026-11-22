"""Verify final curated artifacts against cached bytes and an offline rebuild."""
import argparse
import hashlib
import json
import logging
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from benchmark.published_clock_2026_10_05.curated.retrieve_2026_10_05 import ROOT


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(root, repository, rebuilt):
    """Check evidence, cache receipts, atlas identity and copied output equality."""
    manifest = json.loads((root / 'sources.json').read_text())
    cache = Path(manifest['cache_folder'])
    records = [json.loads(line) for line in (root / 'records.jsonl').read_text().splitlines()]
    if digest(Path(manifest['atlas']['path'])) != manifest['atlas']['sha256']:
        raise ValueError('Atlas differs from the assembly source')
    if len({r['record_id'] for r in records}) != len(records):
        raise ValueError('Record IDs are not unique')
    counts = {
        'records': len(records),
        'by_type': dict(Counter(r['ground_truth_type'] for r in records)),
        'atlas_series_records': sum(r['atlas']['series_present'] for r in records),
        'atlas_series_clock_records': sum(r['atlas']['series_clock_match'] for r in records),
        'direct_case_control_series_clock_records': sum(
            r['direct_case_control_statistic'] and r['atlas']['series_clock_match'] for r in records),
    }
    if counts != manifest['counts']:
        raise ValueError('Manifest counts differ from the final records')
    for row in records:
        for quote in row['quote_evidence']:
            if digest(cache / quote['file']) != quote['sha256']:
                raise ValueError(f"Quote evidence changed: {quote['file']}")
    successful = 0
    error_bodies = 0
    for receipt in manifest['downloads']:
        if 'sha256' in receipt:
            if digest(Path(receipt['path'])) != receipt['sha256']:
                raise ValueError(f"Download bytes changed: {receipt['path']}")
            successful += 1
        if 'error_sha256' in receipt:
            if digest(Path(receipt['error_path'])) != receipt['error_sha256']:
                raise ValueError(f"Cached error body changed: {receipt['error_path']}")
            error_bodies += 1
    for entry in manifest['derived_files']:
        if digest(Path(entry['path'])) != entry['sha256']:
            raise ValueError(f"Derived bytes changed: {entry['path']}")
        if not (cache / entry['parent_file']).exists():
            raise ValueError(f"Missing derivation parent: {entry['parent_file']}")
    hashes = {}
    for name in ['records.jsonl', 'sources.json', 'README.md']:
        hashes[name] = digest(root / name)
        if hashes[name] != digest(repository / name):
            raise ValueError(f'Repository copy differs: {name}')
    if digest(rebuilt / 'records.jsonl') != hashes['records.jsonl']:
        raise ValueError('Offline rebuild records differ')
    rebuilt_manifest = json.loads((rebuilt / 'sources.json').read_text())
    for key in ['counts', 'sources', 'downloads', 'derived_files', 'atlas']:
        if rebuilt_manifest[key] != manifest[key]:
            raise ValueError(f'Offline rebuild provenance differs: {key}')
    return {'verified_at': datetime.now(timezone.utc).isoformat(), 'result': 'passed',
            'counts': counts, 'download_receipts': len(manifest['downloads']),
            'successful_response_hashes_verified': successful,
            'cached_error_body_hashes_verified': error_bodies,
            'derived_file_hashes_verified': len(manifest['derived_files']),
            'file_hashes': hashes, 'repository_copy_equal': True,
            'offline_rebuild_records_byte_identical': True,
            'offline_rebuild_provenance_equal_except_output_folder_and_timestamp': True,
            'unique_series': len({r['gse_id'] for r in records}),
            'unique_atlas_series': len({r['gse_id'] for r in records if r['atlas']['series_present']}),
            'direct_matched_series_clock_pairs': len({(r['gse_id'], r['clock_clockbase']) for r in records
                                                     if r['direct_case_control_statistic']
                                                     and r['atlas']['series_clock_match']})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--repository', type=Path, default=Path(__file__).parent)
    parser.add_argument('--rebuilt', type=Path, default=ROOT / 'offline-rebuild-2026-10-05')
    parser.add_argument('--write-report', action='store_true')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    try:
        report = verify(args.root, args.repository, args.rebuilt)
        payload = json.dumps(report, indent=2) + '\n'
        if args.write_report:
            paths = [args.root / 'validation_2026_10_05.json', args.repository / 'validation_2026_10_05.json']
            if any(path.exists() for path in paths):
                raise FileExistsError('Preserve existing validation reports; omit --write-report to verify again')
            for path in paths:
                with path.open('x') as handle:
                    handle.write(payload)
        print(payload, end='')
    except (OSError, ValueError, KeyError) as exc:
        logging.error('Curated verification failed: %s', exc)
        raise


if __name__ == '__main__':
    main()
