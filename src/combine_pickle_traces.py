"""Combine request-bound trace records; preserve unresolved records in a ledger.

Saved status='success' is not a study-identity check or scientific validation.
The original input archives remain unchanged. No archived code is executed.
"""
import sys
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_system import paths

import os, glob, pickle, datetime
import argparse
import copy
import hashlib
import json
import re

try:
    from .detailed_trace_inspection import SourceOnlyUnpickler, verify_trace_request
except ImportError:  # Direct `python src/combine_pickle_traces.py` entry point.
    from detailed_trace_inspection import SourceOnlyUnpickler, verify_trace_request


def trace_binding(study_id, record):
    if (not isinstance(study_id, str) or not re.fullmatch(r'GSE[0-9]+', study_id)
            or not isinstance(record, dict) or record.get('study_id') != study_id):
        return {'status': 'FAIL', 'issues': ['Archive key and declared study must be the same canonical GSE']}
    result = record.get('result')
    messages = result[0] if isinstance(result, tuple) and result else result
    return verify_trace_request(messages, study_id)

def get_trace_priority(trace):
    """Get priority of a trace based on status (lower number = higher priority)"""
    if isinstance(trace, dict):
        status = trace.get('status', 'unknown')
        if status == 'success':
            return 1
        elif status == 'error':
            return 2
        elif status == 'skipped':
            return 3
    return 4  # unknown/other

def combine_pickle_files(input_dir=None, audit_records=None, max_bytes=2 * 1024**3):
    """Combine all pickle files from result_prompt3 directory, keeping only most informative traces"""
    base_dir = paths.DATA
    result_dir = "result_prompt3"
    dir_path = os.fspath(input_dir) if input_dir is not None else os.path.join(base_dir, result_dir)
    records = audit_records if audit_records is not None else []
    
    if not os.path.exists(dir_path):
        print(f"Directory not found: {dir_path}")
        return {}, 0, 0
    
    all_data = {}
    file_count = 0
    total_traces = 0
    
    # Find all pickle files
    pkl_files = sorted(glob.glob(os.path.join(dir_path, "*.pkl")))
    print(f"Found {len(pkl_files)} pickle files in {result_dir}")
    
    for pkl_file in pkl_files:
        try:
            if not 0 < os.path.getsize(pkl_file) <= max_bytes:
                raise ValueError('Trace archive is empty or exceeds the file byte limit')
            source_hash = hashlib.sha256(Path(pkl_file).read_bytes()).hexdigest()
            with open(pkl_file, 'rb') as f:
                data = SourceOnlyUnpickler(f).load()
                if f.read(1):
                    raise ValueError('Unparsed trailing archive content')
            if not isinstance(data, dict):
                raise ValueError('Trace archive must be a dictionary')
            if hashlib.sha256(Path(pkl_file).read_bytes()).hexdigest() != source_hash:
                raise ValueError('Trace archive changed during loading')
            
            filename = os.path.basename(pkl_file)
            print(f"Loaded: {filename}")
            file_count += 1
            
            # Process data
            if isinstance(data, dict):
                for study_id, study_data in data.items():
                    total_traces += 1
                    binding = trace_binding(study_id, study_data)
                    records.append({'source_file': filename, 'source_sha256': source_hash,
                                    'study_id': study_id if isinstance(study_id, str) else None,
                                    'request_binding': binding,
                                    'disposition': 'ELIGIBLE_FOR_MERGE' if binding['status'] == 'PASS' else 'WITHHELD',
                                    'scientific_status': 'CANNOT VERIFY'})
                    if binding['status'] != 'PASS':
                        continue
                    
                    if study_id not in all_data:
                        # First time seeing this study
                        all_data[study_id] = {
                            'best_trace': copy.deepcopy(study_data),
                            'source': {
                                'file': filename,
                                'dir': result_dir
                            }
                        }
                    else:
                        # Compare with existing trace and keep the more informative one
                        current_priority = get_trace_priority(all_data[study_id]['best_trace'])
                        new_priority = get_trace_priority(study_data)
                        
                        if new_priority < current_priority:
                            # New trace is more informative
                            all_data[study_id] = {
                                'best_trace': copy.deepcopy(study_data),
                                'source': {
                                    'file': filename,
                                    'dir': result_dir
                                }
                            }
                    
        except Exception as e:
            records.append({'source_file': os.path.basename(pkl_file), 'status': 'FAIL',
                            'error': type(e).__name__ + ': ' + str(e), 'disposition': 'WITHHELD',
                            'scientific_status': 'CANNOT VERIFY'})
            print(f"Error loading {pkl_file}: {e}")
    
    # Convert to simpler format for saving
    simplified_data = {}
    for study_id, data in all_data.items():
        simplified_data[study_id] = data['best_trace']
    
    return simplified_data, file_count, total_traces

# Run the combination
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--max-file-bytes', type=int, default=2 * 1024**3)
    args = parser.parse_args()
    if args.max_file_bytes <= 0:
        parser.error('max-file-bytes must be positive')
    if args.output_dir.exists():
        parser.error('Output directory already exists; originals are never overwritten')
    args.output_dir.mkdir(parents=True)
    print("Combining agent traces from result_prompt3 directory (keeping only most informative)...")
    audit_records = []
    combined_data, total_files, total_traces = combine_pickle_files(args.input_dir, audit_records, args.max_file_bytes)
    with (args.output_dir / 'trace-binding-ledger.json').open('x') as handle:
        json.dump({'records': audit_records, 'scope': 'Request binding only; no scientific acceptance',
                   'tie_break': 'Sorted input filenames, then status priority; no outcome-based selection',
                   'input_dir': str(args.input_dir.resolve())}, handle, indent=2)
    
    if total_files == 0:
        print("No files processed. Exiting.")
        exit(1)
    
    # Save combined data
    output_dir = os.fspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"combined_traces_prompt3_deduped_{timestamp}.pkl")
    
    with open(output_file, 'xb') as f:
        pickle.dump(combined_data, f)
    
    # Count by status
    status_counts = {}
    for study_id, trace in combined_data.items():
        if isinstance(trace, dict):
            status = trace.get('status', 'unknown')
            status_counts[status] = status_counts.get(status, 0) + 1
    
    print(f"\nSummary:")
    print(f"Total pickle files processed: {total_files}")
    print(f"Total traces found: {total_traces}")
    print(f"Unique studies after deduplication: {len(combined_data)}")
    print(f"Reduction: {total_traces - len(combined_data)} duplicate traces removed")
    
    print(f"\nStatus distribution:")
    for status, count in sorted(status_counts.items()):
        print(f"  {status}: {count}")
    
    print(f"\nCombined data saved to: {output_file}")
    
    # Create a summary report
    summary_file = os.path.join(output_dir, f"combination_summary_prompt3_deduped_{timestamp}.txt")
    result_dir_path = os.fspath(args.input_dir)
    with open(summary_file, 'x') as f:
        f.write(f"Agent Traces Combination Summary (result_prompt3 - deduplicated)\n")
        f.write(f"=" * 60 + "\n")
        f.write(f"Generated at: {datetime.datetime.now()}\n\n")
        f.write(f"Directory processed: {result_dir_path}\n")
        f.write(f"Total files processed: {total_files}\n")
        f.write(f"Total traces found: {total_traces}\n")
        f.write(f"Unique studies after deduplication: {len(combined_data)}\n")
        f.write(f"Duplicate traces removed: {total_traces - len(combined_data)}\n\n")
        
        f.write("Status distribution of kept traces:\n")
        for status, count in sorted(status_counts.items()):
            f.write(f"  {status}: {count}\n")
    
    print(f"Summary report saved to: {summary_file}")
