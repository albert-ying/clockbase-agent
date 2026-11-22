import json
import pandas as pd
from pathlib import Path


main_csv = '/root/ClockBase-Paper-Relevance/src/data/main.csv'
report_dir = '/root/ClockBase-Paper-Relevance/data/output'
sample_file_name = 'sample.json'


cases = Path(report_dir).glob('*')

df = pd.read_csv(main_csv)

cases_1000 = []

for case in cases:
    case_id = case.name
    sample_file = case / sample_file_name
    with open(str(sample_file), 'r') as f:
        sample = json.load(f)
    gse_id = sample['gse_id']
    comparison_name = sample['comparison_name']
    
    # Use boolean indexing instead of query() to avoid quote escaping issues
    data = df[(df['gse_id'] == gse_id) & (df['comparison_name'] == comparison_name)]
    
    if len(data) > 0:
        data = data.iloc[0].to_dict()
        data['case_id'] = case_id
        cases_1000.append(data)
    else:
        print(f'{case_id} does not exist')

cases_1000 = pd.DataFrame(cases_1000)

# Move case_id column to the first position
cols = cases_1000.columns.tolist()
cols.remove('case_id')
cases_1000 = cases_1000[['case_id'] + cols]

cases_1000.to_csv('case_mapping.csv', index=False)