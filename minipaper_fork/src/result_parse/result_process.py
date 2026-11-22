import pandas as pd
from pathlib import Path
import json

def rank_result(output_dir: str, main_file: str) -> pd.DataFrame:
    relevance_file = 'identified_relevance.json'
    sample_file = 'sample.json'
    od = Path(output_dir)
    data = []
    df = pd.read_csv(main_file)
    for case in od.glob('*'):
        try:
            case_id = case.name
            relevance_path = case / relevance_file
            sample_path = case / sample_file
            with open(relevance_path, 'r') as f:
                relevance = json.load(f)
                final_score = relevance['final_score']
            with open(sample_path, 'r') as f:
                sample = json.load(f)
                gse_id = sample['gse_id']
                comparison_name = sample['comparison_name']
                pvalue = sample['pvalue']
            # Use boolean indexing instead of query() to avoid quote escaping issues
            target = df[(df['gse_id'] == gse_id) & (df['comparison_name'] == comparison_name)]
            if len(target) > 0:
                target = target.iloc[0].to_dict()
                fdr = target['fdr']
            else:
                print(f'{case_id} does not exist')
                continue
            data.append({'case_id': case_id, 'final_score': final_score, 'gse_id': gse_id, 'comparison_name': comparison_name, 'pvalue': pvalue, 'fdr': fdr})
        except Exception as e:
            print(f'Error parsing {case_id}: {e}')
            continue
    df = pd.DataFrame(data)
    df = df.sort_values(by='final_score', ascending=False)
    return df


def filter_fdr(df: pd.DataFrame, fdr_threshold: float = 0.05) -> pd.DataFrame:
    df = df[df['fdr'] < fdr_threshold]
    return df


if __name__ == "__main__":
    # 1. rank result
    # output_dir = '/home/kido/avinasi/clockbase/data/output'
    # main_file = '/home/kido/avinasi/clockbase/data/main.csv'
    # df = rank_result(output_dir, main_file)
    # df.to_csv('/home/kido/avinasi/clockbase/data/output_df.csv', index=False)

    # 2. filter pvalue
    df = pd.read_csv('/home/kido/avinasi/clockbase/data/output_df.csv')
    df = filter_fdr(df)
    df.to_csv('/home/kido/avinasi/clockbase/data/output_df_fdr.csv', index=False)