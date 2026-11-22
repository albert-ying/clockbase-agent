import sys
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_system import paths

import os
import pickle
import datetime
from collections import Counter

def load_combined_traces(filepath):
    """Load the combined traces pickle file"""
    with open(filepath, 'rb') as f:
        data = pickle.load(f)
    return data

def analyze_traces(combined_data):
    """Analyze the combined agent traces"""
    
    # Initialize analysis results
    analysis = {
        'total_studies': len(combined_data),
        'studies_with_multiple_traces': 0,
        'status_distribution': Counter(),
        'execution_times': [],
        'error_messages': [],
        'successful_studies': [],
        'failed_studies': [],
        'skipped_studies': [],
        'trace_counts': Counter()
    }
    
    # Analyze each study
    for study_id, study_data in combined_data.items():
        traces = study_data['traces']
        sources = study_data['sources']
        
        # Count traces per study
        trace_count = len(traces)
        analysis['trace_counts'][trace_count] += 1
        
        if trace_count > 1:
            analysis['studies_with_multiple_traces'] += 1
        
        # Analyze each trace
        for trace in traces:
            if isinstance(trace, dict):
                # Get status
                status = trace.get('status', 'unknown')
                analysis['status_distribution'][status] += 1
                
                # Collect execution times
                if 'execution_time' in trace:
                    analysis['execution_times'].append(trace['execution_time'])
                
                # Categorize by status
                if status == 'success':
                    analysis['successful_studies'].append(study_id)
                elif status == 'error':
                    analysis['failed_studies'].append(study_id)
                    if 'message' in trace:
                        analysis['error_messages'].append({
                            'study_id': study_id,
                            'error': trace['message']
                        })
                elif status == 'skipped':
                    analysis['skipped_studies'].append(study_id)
    
    # Remove duplicates from study lists
    analysis['successful_studies'] = list(set(analysis['successful_studies']))
    analysis['failed_studies'] = list(set(analysis['failed_studies']))
    analysis['skipped_studies'] = list(set(analysis['skipped_studies']))
    
    return analysis

def generate_report(analysis, output_dir):
    """Generate a detailed analysis report"""
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = os.path.join(output_dir, f"detailed_analysis_report_{timestamp}.txt")
    
    with open(report_file, 'w') as f:
        f.write("Detailed Agent Traces Analysis Report\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Generated at: {datetime.datetime.now()}\n\n")
        
        # Basic statistics
        f.write("Basic Statistics:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Total unique studies: {analysis['total_studies']}\n")
        f.write(f"Studies with multiple traces: {analysis['studies_with_multiple_traces']}\n\n")
        
        # Status distribution
        f.write("Status Distribution:\n")
        f.write("-" * 30 + "\n")
        for status, count in analysis['status_distribution'].most_common():
            f.write(f"  {status}: {count}\n")
        f.write("\n")
        
        # Trace count distribution
        f.write("Trace Count Distribution:\n")
        f.write("-" * 30 + "\n")
        for trace_count, study_count in sorted(analysis['trace_counts'].items()):
            f.write(f"  {trace_count} trace(s): {study_count} studies\n")
        f.write("\n")
        
        # Execution time statistics
        if analysis['execution_times']:
            f.write("Execution Time Statistics:\n")
            f.write("-" * 30 + "\n")
            exec_times = analysis['execution_times']
            f.write(f"  Total executions: {len(exec_times)}\n")
            f.write(f"  Average time: {sum(exec_times)/len(exec_times):.2f} seconds\n")
            f.write(f"  Min time: {min(exec_times):.2f} seconds\n")
            f.write(f"  Max time: {max(exec_times):.2f} seconds\n")
            f.write(f"  Total time: {sum(exec_times):.2f} seconds ({sum(exec_times)/3600:.2f} hours)\n\n")
        
        # Study counts by status
        f.write("Study Counts by Status:\n")
        f.write("-" * 30 + "\n")
        f.write(f"  Successful studies: {len(analysis['successful_studies'])}\n")
        f.write(f"  Failed studies: {len(analysis['failed_studies'])}\n")
        f.write(f"  Skipped studies: {len(analysis['skipped_studies'])}\n\n")
        
        # Sample of error messages
        if analysis['error_messages']:
            f.write("Sample Error Messages (first 10):\n")
            f.write("-" * 30 + "\n")
            for i, error_info in enumerate(analysis['error_messages'][:10]):
                f.write(f"  {error_info['study_id']}: {error_info['error'][:100]}...\n")
            f.write("\n")
    
    print(f"Detailed report saved to: {report_file}")
    
    # Also save lists of studies by status
    for status_type in ['successful', 'failed', 'skipped']:
        studies = analysis[f'{status_type}_studies']
        if studies:
            status_file = os.path.join(output_dir, f"{status_type}_studies_{timestamp}.txt")
            with open(status_file, 'w') as f:
                f.write(f"{status_type.capitalize()} Studies\n")
                f.write("=" * 30 + "\n")
                for study_id in sorted(studies):
                    f.write(f"{study_id}\n")
            print(f"Saved {status_type} studies list to: {status_file}")

if __name__ == "__main__":
    # Find the most recent combined traces file
    combined_dir = paths.DATA / "combined_results"
    
    import glob
    pkl_files = glob.glob(os.path.join(combined_dir, "combined_traces_*.pkl"))
    
    if not pkl_files:
        print("No combined traces files found!")
        exit(1)
    
    # Use the most recent file
    latest_file = max(pkl_files, key=os.path.getctime)
    print(f"Analyzing: {latest_file}")
    
    # Load and analyze
    combined_data = load_combined_traces(latest_file)
    analysis = analyze_traces(combined_data)
    
    # Generate report
    generate_report(analysis, combined_dir)
    
    print("\nAnalysis complete!") 