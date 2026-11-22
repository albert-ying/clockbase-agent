import json
import os
import pandas as pd
from pathlib import Path

def extract_scores_from_json(json_path):
    """Extract all scores from a single JSON file."""
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        scores = {}
        
        # Extract all question scores
        for key, value in data.items():
            if isinstance(value, dict) and 'points' in value:
                scores[key] = value['points']
            elif key == 'q0_plausibility' and 'penalty_or_bonus_applied' in value:
                scores['q0_plausibility'] = value['penalty_or_bonus_applied']
        
        # Add final score and intervention
        if 'final_score' in data:
            scores['final_score'] = data['final_score']
        if 'q6_intervention' in data:
            scores['intervention'] = data['q6_intervention']
        
        return scores
    except Exception as e:
        print(f"Error processing {json_path}: {e}")
        return None

def main():
    # Path to the output directory
    output_dir = Path('data/output')
    
    # Collect all scores
    all_scores = []
    
    # Iterate through all subdirectories
    for subdir in sorted(output_dir.iterdir()):
        if subdir.is_dir():
            json_file = subdir / 'identified_relevance.json'
            
            if json_file.exists():
                scores = extract_scores_from_json(json_file)
                if scores:
                    # Add dataset name
                    scores['dataset'] = subdir.name
                    all_scores.append(scores)
                    print(f"Processed: {subdir.name}")
            else:
                print(f"No identified_relevance.json found in {subdir.name}")
    
    # Create DataFrame
    df = pd.DataFrame(all_scores)
    
    # Reorder columns - put dataset first, then scores in order, intervention, final score last
    score_columns = [col for col in df.columns if col.startswith('q') and col != 'dataset']
    score_columns.sort()  # Sort question columns
    
    column_order = ['dataset'] + score_columns
    if 'intervention' in df.columns:
        column_order.append('intervention')
    if 'final_score' in df.columns:
        column_order.append('final_score')
    
    # Remove any columns not in our desired order
    column_order = [col for col in column_order if col in df.columns]
    df = df[column_order]
    
    # Sort by dataset name
    df = df.sort_values('dataset')
    
    # Save to CSV
    output_file = 'data/scores_summary.csv'
    df.to_csv(output_file, index=False)
    print(f"\nSummary saved to {output_file}")
    
    # Display basic statistics
    print(f"\nTotal datasets processed: {len(df)}")
    print(f"\nScore columns found: {[col for col in df.columns if col.startswith('q')]}")
    
    if 'final_score' in df.columns:
        print(f"\nFinal score statistics:")
        print(f"  Mean: {df['final_score'].mean():.2f}")
        print(f"  Median: {df['final_score'].median():.2f}")
        print(f"  Min: {df['final_score'].min()}")
        print(f"  Max: {df['final_score'].max()}")

if __name__ == "__main__":
    main()