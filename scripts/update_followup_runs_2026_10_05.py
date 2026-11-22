#!/usr/bin/env python3
"""Reconcile RNA80 controller completion with its saved follow-up evidence."""
import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from audit_system import paths


EVIDENCE_RELATIVE = Path(
    'paper/revision_2026-09-07/evidence/'
    'continuous-audit-pipeline-2026-09-20'
)


def read_json(path):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f'Could not read evidence JSON {path}: {error}') from error


def sha256(path):
    digest = hashlib.sha256()
    try:
        with path.open('rb') as source:
            for block in iter(lambda: source.read(1024 * 1024), b''):
                digest.update(block)
    except OSError as error:
        raise ValueError(f'Could not hash evidence file {path}: {error}') from error
    return digest.hexdigest()


def relative(path, root):
    return path.relative_to(root).as_posix()


def source_record(path, evidence_root):
    try:
        modified = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()
    except OSError as error:
        raise ValueError(f'Could not stat evidence file {path}: {error}') from error
    return {'path': path.relative_to(evidence_root).as_posix(),
            'sha256': sha256(path), 'filesystem_mtime_utc': modified}


def read_study_status(path):
    try:
        with path.open(newline='') as source:
            return {row['study_id']: row for row in csv.DictReader(source)}
    except (OSError, KeyError, csv.Error) as error:
        raise ValueError(f'Could not read study-status CSV {path}: {error}') from error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-root', type=Path,
                        default=paths.REPO / EVIDENCE_RELATIVE,
                        help='Evidence directory; defaults to the revision evidence tree.')
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='New scratch directory for JSON and TSV output.')
    args = parser.parse_args()

    evidence = args.evidence_root.resolve()
    controller = evidence / 'update-controller'
    recovery = controller / 'remaining-rna-recovery-2026-09-23'
    if args.output_dir.exists():
        raise ValueError(f'Output directory already exists; choose a new scratch path: {args.output_dir}')

    first_launch_path = controller / 'downstream-rna-preparation-001' / 'report-launch-001.json'
    first_readback_path = controller / 'downstream-rna-preparation-001' / 'completed-readback-002' / 'summary.json'
    if not first_readback_path.exists():
        first_readback_path = controller / 'downstream-rna-preparation-001' / 'completed-readback-002' / 'studies.csv'
        rows = read_study_status(first_readback_path)
        initial_reports = sum(row['pipeline_status'] == 'COMPLETED' for row in rows.values())
        initial_audits = sum(row['hook_status'] == 'COMPLETED' for row in rows.values())
        initial_records = sum(int(row['source_samples']) for row in rows.values())
        initial_studies = len(rows)
        first_model = None
        first_designs = None
        first_report_records = sum(int(row['source_samples']) for row in rows.values()
                                   if row['pipeline_status'] == 'COMPLETED')
        first_audit_records = sum(int(row['source_samples']) for row in rows.values()
                                  if row['hook_status'] == 'COMPLETED')
    else:
        readback = read_json(first_readback_path)
        initial_reports = readback['completed_stages']['08_report']
        initial_audits = readback['comprehensive_status_counts']['COMPLETED']
        initial_records = readback['source_samples']
        initial_studies = readback['studies']
        first_model = ', '.join(sorted(readback.get('requested_models', {}))) or None
        first_designs = readback['completed_stages']['03_blinded_design']
        first_rows_path = controller / 'downstream-rna-preparation-001' / 'completed-readback-002' / 'studies.csv'
        first_rows = list(csv.DictReader(first_rows_path.open(newline='')))
        first_report_records = sum(int(row['source_samples']) for row in first_rows
                                   if row['pipeline_status'] == 'COMPLETED')
        first_audit_records = sum(int(row['source_samples']) for row in first_rows
                                  if row['hook_status'] == 'COMPLETED')
    first_launch = read_json(first_launch_path) if first_launch_path.exists() else {}
    first_rows_path = controller / 'downstream-rna-preparation-001' / 'completed-readback-002' / 'studies.csv'
    first_rows = list(csv.DictReader(first_rows_path.open(newline='')))
    initial_report_studies = sorted(row['study_id'] for row in first_rows
                                    if row['pipeline_status'] == 'COMPLETED')
    initial_audit_studies = sorted(row['study_id'] for row in first_rows
                                   if row['hook_status'] == 'COMPLETED')

    format_launch_path = controller / 'format-recovery-development-002' / 'report-launch-001.json'
    format_launch = read_json(format_launch_path)
    format_plan_path = controller / 'ingestion-plan-015.json'
    format_plan = read_json(format_plan_path)
    format_studies = ['GSE203635', 'GSE227798', 'GSE240194']
    format_terminals = []
    format_registry = controller / 'registry-014' / 'studies' / 'human_rna'
    for study in format_studies:
        versions = sorted((format_registry / study).glob('*/pipeline-report.json'))
        if len(versions) != 1:
            raise ValueError(f'Expected one registry-014 report terminal for {study}, found {len(versions)}')
        report_terminal = read_json(versions[0])
        audit_path = versions[0].parent / 'comprehensive-audit.json'
        audit_terminal = read_json(audit_path)
        if report_terminal.get('status') != 'COMPLETED' or audit_terminal.get('status') != 'COMPLETED':
            raise ValueError(f'Registry-014 terminals are not complete for {study}')
        format_terminals.append((study, versions[0], audit_path))

    five_rows_path = recovery / 'completion-readback-001' / 'run-002' / 'studies.json'
    five_receipt_path = recovery / 'completion-readback-001' / 'run-002' / 'receipt.json'
    five_rows = read_json(five_rows_path)
    five_studies = sorted({item['study_id'] for item in five_rows
                           if item.get('comprehensive_summary', {}).get('execution_status') == 'COMPLETED'})

    six_protocol_path = recovery / 'producer-output-completion-001' / 'completion-readback-001' / 'PROTOCOL.md'
    six_report_path = recovery / 'producer-output-completion-001' / 'completion-readback-001' / 'REPORT.md'
    five_launch_path = recovery / 'execution-launch-001' / 'launch.json'
    five_freeze_path = recovery / 'execution-launch-001' / 'execution-freeze.json'
    six_launch_path = recovery / 'producer-output-completion-001' / 'execution-launch-001' / 'launch.json'
    six_freeze_path = recovery / 'producer-output-completion-001' / 'execution-launch-001' / 'execution-freeze-002.json'
    five_launch = read_json(five_launch_path)
    five_freeze = read_json(five_freeze_path)
    five_plan_path = recovery / 'cohort-preparation-001' / 'prepared-001' / 'ingestion-plan.json'
    five_plan = read_json(five_plan_path)
    six_launch = read_json(six_launch_path)
    six_freeze = read_json(six_freeze_path)
    six_readback_path = recovery / 'producer-output-completion-001' / 'completion-readback-001' / 'run-001' / 'studies.json'
    six_readback = read_json(six_readback_path)
    six_studies = sorted(item['study_id'] for item in six_readback
                         if item.get('comprehensive_summary', {}).get('execution_status') == 'COMPLETED')

    audit_summary_path = recovery / 'backlog-after-source-followup-001' / 'summary.json'
    audit_summary = read_json(audit_summary_path)
    audit_studies = sorted(audit_summary['completed_audit_only_followup_records'])
    audit_launch_dir = recovery / 'source-stage-followup-implementation-001' / 'execution-launch-001' / 'run-001'
    audit_processes = [read_json(audit_launch_dir / study / 'process.json') for study in audit_studies]
    audit_freeze_path = recovery / 'source-stage-followup-implementation-001' / 'execution-launch-001' / 'execution-freeze.json'
    audit_freeze = read_json(audit_freeze_path)

    previous_status_path = recovery / 'backlog-after-source-followup-001' / 'study-status.csv'
    producer_status_path = recovery / 'backlog-after-remaining-producer-001' / 'study-status.csv'
    previous_status = read_study_status(previous_status_path)
    producer_status = read_study_status(producer_status_path)
    producer_studies = sorted(study for study, row in producer_status.items()
                              if row['producer_report_completed'].lower() == 'true'
                              and previous_status[study]['producer_report_completed'].lower() != 'true')
    producer_audits = sorted(study for study, row in producer_status.items()
                            if row['comprehensive_audit_completed'].lower() == 'true'
                            and previous_status[study]['comprehensive_audit_completed'].lower() != 'true')
    remaining_launch_dir = recovery / 'remaining-producer-readiness-2026-09-23-001' / 'execution-launch-001' / 'run-001'
    remaining_processes = [read_json(remaining_launch_dir / shard / 'process.json')
                           for shard in ('shard-000', 'shard-001')]
    remaining_freeze_path = recovery / 'remaining-producer-readiness-2026-09-23-001' / 'execution-launch-001' / 'execution-freeze.json'
    remaining_freeze = read_json(remaining_freeze_path)

    gpt_base = recovery / 'gpt6-development-2026-09-24'
    two_summary_path = gpt_base / 'incremental-census-001' / 'run-002' / 'summary.json'
    two_summary = read_json(two_summary_path)
    two_freeze_path = gpt_base / 'allocation-freeze-001.json'
    two_freeze = read_json(two_freeze_path)
    two_start_path = gpt_base / 'execution-evidence-001' / 'execution-001' / 'GSE290668' / 'start.json'
    two_start = read_json(two_start_path)
    two_studies = sorted(two_summary['new_completed_producer_accessions'])
    continuation_summary_path = gpt_base / 'continuation-v3-001' / 'incremental-census-001' / 'run-001' / 'summary.json'
    continuation_summary = read_json(continuation_summary_path)
    continuation_base = gpt_base / 'continuation-v3-001'
    continuation_freeze_path = continuation_base / 'protocol-freeze-001.json'
    continuation_freeze = read_json(continuation_freeze_path)
    continuation_start_path = continuation_base / 'execution-evidence-001' / 'execution-001' / 'start.json'
    continuation_start = read_json(continuation_start_path)
    continuation_studies = sorted(continuation_summary['new_completed_producer_accessions'])
    final_summary_path = gpt_base / 'coding-contract-followup-001' / 'incremental-census-001' / 'run-002' / 'summary.json'
    final_summary = read_json(final_summary_path)
    final_freeze_path = gpt_base / 'coding-contract-followup-001' / 'report-freeze-001' / 'report-freeze-001' / 'freeze.json'
    final_freeze = read_json(final_freeze_path)
    final_launch_path = gpt_base / 'coding-contract-followup-001' / 'completed-report-readback-001' / 'execution-launch-001' / 'process.json'
    final_launch = read_json(final_launch_path)
    final_studies = sorted(set(final_summary['new_completed_report_accessions']) |
                           set(final_summary['new_completed_audit_accessions']))

    # GSE127814 was complete before the first batch launch. It is imported into
    # the roster count but is not a follow-up run.
    pre_batch_import_path = controller / 'full-run-execution-007.json'
    pre_batch_import = read_json(pre_batch_import_path)
    prior_imports = {'studies': ['GSE127814'], 'reports_added': 1, 'audits_added': 1,
                     'source': source_record(pre_batch_import_path, evidence),
                     'completed_at_utc': pre_batch_import['completed_at_utc'],
                     'classification': 'pre-batch completion; excluded from follow-up run count'}
    report_studies = set(initial_report_studies) | set(prior_imports['studies'])
    audit_studies_set = set(initial_audit_studies) | set(prior_imports['studies'])

    # These cumulative values are derived from the initial run, the prior
    # completion import, and each measured study-set increment below.
    records = [
        {'group': 'initial_controller_run', 'directory': first_launch_path.parent,
         'launch_date': first_launch.get('started_at_utc'), 'model': first_model,
         'report_studies': initial_report_studies, 'audit_studies': initial_audit_studies,
         'process_count': 1, 'reports_added': 0, 'audits_added': 0,
         'cumulative_reports': initial_reports, 'cumulative_audits': initial_audits,
         'evidence': [first_launch_path, first_readback_path]},
        {'group': 'format-recovery-development-002', 'directory': format_launch_path.parent,
         'launch_date': format_launch.get('recorded_at_utc'),
         'readback_date': None,
         'model': format_plan['pipeline']['model'], 'report_studies': format_studies,
         'audit_studies': format_studies, 'process_count': 1,
         'reports_added': len(format_studies), 'audits_added': len(format_studies),
         'evidence': [format_launch_path, format_plan_path] +
                     [path for _, report, audit in format_terminals for path in (report, audit)]},
        {'group': 'five-study-recovery', 'directory': five_launch_path.parent,
         'launch_date': five_launch.get('started_at_utc'), 'readback_date': None,
         'model': five_plan.get('pipeline', {}).get('model'), 'report_studies': five_studies,
         'audit_studies': five_studies, 'process_count': 1,
         'reports_added': len(five_studies), 'audits_added': len(five_studies),
         'evidence': [five_launch_path, five_freeze_path, five_plan_path, five_rows_path, five_receipt_path]},
        {'group': 'six-case-producer-followup',
         'directory': (recovery / 'producer-output-completion-001' / 'execution-launch-001'),
         'launch_date': six_launch.get('started_at_utc'), 'readback_date': None,
         'model': six_freeze.get('model'), 'report_studies': six_studies,
         'audit_studies': six_studies, 'process_count': 1,
         'reports_added': len(six_studies), 'audits_added': len(six_studies),
         'evidence': [six_launch_path, six_freeze_path, six_protocol_path, six_report_path, six_readback_path]},
        {'group': 'audit-only-source-followup', 'directory': audit_launch_dir,
         'launch_date': min(process['started_at_utc'] for process in audit_processes),
         'readback_date': audit_summary.get('created_at_utc'), 'model': audit_freeze.get('model'),
         'report_studies': [], 'audit_studies': audit_studies, 'process_count': len(audit_processes),
         'reports_added': 0, 'audits_added': len(audit_studies),
         'evidence': [audit_freeze_path, audit_summary_path] +
                     [audit_launch_dir / study / 'process.json' for study in audit_studies]},
        {'group': 'remaining-producer-followup', 'directory': remaining_launch_dir,
         'launch_date': min(process['started_at_utc'] for process in remaining_processes),
         'readback_date': read_json(producer_status_path.parent / 'summary.json').get('created_at_utc'),
         'model': remaining_freeze.get('model'), 'report_studies': producer_studies,
         'audit_studies': producer_audits, 'process_count': len(remaining_processes),
         'reports_added': len(producer_studies), 'audits_added': len(producer_audits),
         'evidence': [remaining_freeze_path] +
                     [remaining_launch_dir / shard / 'process.json' for shard in ('shard-000', 'shard-001')] +
                     [previous_status_path, producer_status_path,
                      producer_status_path.parent / 'summary.json']},
        {'group': 'gpt6-development-two-study', 'directory': two_start_path.parent.parent,
         'launch_date': two_start.get('timestamp_utc'),
         'readback_date': two_summary.get('recorded_at_utc'),
         'model': two_start.get('model'), 'report_studies': two_studies,
         'audit_studies': two_studies, 'process_count': 1,
         'reports_added': len(two_studies), 'audits_added': len(two_studies),
         'evidence': [two_freeze_path, two_start_path, two_summary_path]},
        {'group': 'continuation-v3-thirteen-study',
         'directory': continuation_start_path.parent,
         'launch_date': continuation_start.get('at_utc'),
         'readback_date': continuation_summary.get('recorded_at_utc'),
         'model': continuation_start.get('model'), 'report_studies': continuation_studies,
         'audit_studies': continuation_studies, 'process_count': 1,
         'reports_added': len(continuation_studies), 'audits_added': len(continuation_studies),
         'evidence': [continuation_freeze_path, continuation_start_path, continuation_summary_path]},
        {'group': 'coding-contract-three-study', 'directory': final_launch_path.parent,
         'launch_date': final_launch.get('started_at_utc'), 'readback_date': final_summary.get('at_utc'),
         'model': final_launch.get('model'),
         'report_studies': sorted(final_summary['new_completed_report_accessions']),
         'audit_studies': sorted(final_summary['new_completed_audit_accessions']), 'process_count': 1,
         'reports_added': len(final_summary['new_completed_report_accessions']),
         'audits_added': len(final_summary['new_completed_audit_accessions']),
         'evidence': [final_launch_path, final_freeze_path, final_summary_path]},
    ]

    cumulative_reports = len(report_studies)
    cumulative_audits = len(audit_studies_set)
    for row in records[1:]:
        if set(row['report_studies']) & report_studies or set(row['audit_studies']) & audit_studies_set:
            raise ValueError(f"A follow-up completion overlaps a prior roster: {row['group']}")
        report_studies.update(row['report_studies'])
        audit_studies_set.update(row['audit_studies'])
        if row['reports_added'] != len(row['report_studies']) or row['audits_added'] != len(row['audit_studies']):
            raise ValueError(f"Study-set deltas do not match reported counts: {row['group']}")
        cumulative_reports = len(report_studies)
        cumulative_audits = len(audit_studies_set)
        row['cumulative_reports'] = cumulative_reports
        row['cumulative_audits'] = cumulative_audits
    if initial_reports != 26 or initial_audits != 24:
        raise ValueError(f'Initial controller readback was {initial_reports}/{initial_audits}, expected 26/24')
    if len(producer_studies) != 25 or len(producer_audits) != 24:
        raise ValueError(f'Remaining-producer snapshot delta is {len(producer_studies)}/{len(producer_audits)}, expected 25/24')
    if len(two_studies) != 1 or len(continuation_studies) != 10 or len(final_studies) != 2:
        raise ValueError('A later census accession delta did not match its checked completion count')
    if final_summary['counts']['producer_report_source_records'] != 875 or final_summary['counts']['comprehensive_audit_source_records'] != 866:
        raise ValueError('Final source-record totals do not match the requested 875/866 census')
    if cumulative_reports != 78 or cumulative_audits != 77:
        raise ValueError(f'Computed cumulative totals are {cumulative_reports}/{cumulative_audits}; expected final census 78/77')
    roster_samples = {row['study_id']: int(row['source_samples']) for row in first_rows}
    if (report_studies | audit_studies_set) - roster_samples.keys():
        raise ValueError('A follow-up completion is outside the original 80-study roster')
    report_records = sum(roster_samples[study] for study in report_studies)
    audit_records = sum(roster_samples[study] for study in audit_studies_set)
    if (report_records != final_summary['counts']['producer_report_source_records']
            or audit_records != final_summary['counts']['comprehensive_audit_source_records']):
        raise ValueError('Study-set source-record totals differ from the final census')

    run_rows = records[1:]
    result = {
        'schema': 'clockbase-update-followup-runs-v1',
        'evidence_root': str(evidence),
        'initial_controller_run': {
            'directory': records[0]['directory'].relative_to(evidence).as_posix(),
            'launch_date': records[0]['launch_date'], 'model': records[0]['model'],
            'studies': len(first_rows), 'source_records': initial_records,
            'designs': first_designs,
            'reports': initial_reports, 'report_records': first_report_records,
            'audits': initial_audits, 'audit_records': first_audit_records,
            'evidence': [source_record(path, evidence) for path in
                         (records[0]['evidence'] + [first_rows_path] if first_readback_path.suffix == '.json'
                          else records[0]['evidence'])],
        },
        'prior_imports': prior_imports,
        'followup_execution_group_count': len(run_rows),
        'report_producing_followup_group_count': sum(row['reports_added'] > 0 for row in run_rows),
        'audit_only_followup_group_count': sum(row['reports_added'] == 0 and row['audits_added'] > 0
                                               for row in run_rows),
        'all_group_launch_dates_recovered': all(row['launch_date'] is not None for row in run_rows),
        'group_count_definition': ('Execution groups: seven report-producing follow-up batches and one audit-only '
                                   'source-follow-up group. The latter has two study processes; the remaining-producer '
                                   'group has two shards. Group count is not a count of individual processes.'),
        'followup_runs': [
            {'group': row['group'],
             'directory': row['directory'].relative_to(evidence).as_posix(),
             'launch_date': row['launch_date'], 'model': row['model'],
             'readback_date': row.get('readback_date'),
             'studies_completed_in_run': sorted(set(row['report_studies']) | set(row['audit_studies'])),
             'reports_completed_in_run': row['report_studies'],
             'audits_completed_in_run': row['audit_studies'],
             'process_count': row['process_count'],
             'reports_added': row['reports_added'], 'audits_added': row['audits_added'],
             'cumulative_reports': row['cumulative_reports'],
             'cumulative_audits': row['cumulative_audits'],
             'evidence': [source_record(path, evidence) for path in row['evidence']]}
            for row in run_rows
        ],
        'final_census': {
            'studies': final_summary['counts']['study_count'],
            'source_records': final_summary['counts']['source_records'],
            'reports': final_summary['counts']['producer_report_studies'],
            'report_records': final_summary['counts']['producer_report_source_records'],
            'audits': final_summary['counts']['comprehensive_audit_studies'],
            'audit_records': final_summary['counts']['comprehensive_audit_source_records'],
            'at_utc': final_summary['at_utc'],
            'evidence': source_record(final_summary_path, evidence),
        },
    }

    output = args.output_dir
    try:
        output.mkdir(parents=True, exist_ok=False)
        json_path = output / 'update_followup_runs.json'
        tsv_path = output / 'update_followup_runs.tsv'
        json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
        fields = ['group', 'directory', 'launch_date', 'model', 'studies_completed_in_run',
                  'reports_completed_in_run', 'audits_completed_in_run', 'process_count',
                  'reports_added', 'audits_added', 'cumulative_reports', 'cumulative_audits']
        with tsv_path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields)
            writer.writeheader()
            for row in result['followup_runs']:
                writer.writerow({field: (','.join(row[field])
                                         if field in {'studies_completed_in_run', 'reports_completed_in_run',
                                                      'audits_completed_in_run'} else row[field])
                                 for field in fields})
    except OSError as error:
        raise ValueError(f'Could not write Part C outputs under {output}: {error}') from error
    print(json.dumps({'json': str(json_path), 'tsv': str(tsv_path),
                      'followup_execution_group_count': result['followup_execution_group_count'],
                      'report_producing_followup_group_count': result['report_producing_followup_group_count'],
                      'audit_only_followup_group_count': result['audit_only_followup_group_count'],
                      'final_census': result['final_census']}, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise SystemExit(f'update_followup_runs_2026_10_05: {error}') from error
