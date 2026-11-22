"""
ClockBase E2E orchestration script

Chains the 6 stages of the full processing pipeline into a single run:
1. Batch processing (pipeline)      → run.main()
2. Scoring & FDR filtering          → result_process.rank_result() + filter_fdr()
3. Case selection & copying          → run.select_cases()
4. Report generation                 → run.process_top()
5. Reference cleaning                → clean.main()
6. HTML conversion                   → convert.markdown_to_html()

Usage:
    cd src && python orchestrate.py --df-path ../data/main.csv
    cd src && python orchestrate.py --df-path ../data/main.csv --workdir /tmp/clockbase_run
    cd src && python orchestrate.py --df-path ../data/main.csv --start-stage 3 --top-n 500
"""

import argparse
import asyncio
from pathlib import Path

import pandas as pd

from run import main as run_pipeline, select_cases, process_top
from result_parse.result_process import rank_result, filter_fdr
from result_parse.clean import main as clean_reports
from result_parse.convert import markdown_to_html
from utils.logger import logger


STAGES = {
    1: 'Batch processing (pipeline)',
    2: 'Scoring & FDR filtering (rank + filter)',
    3: 'Select Top-N cases (select_cases)',
    4: 'Generate reports (generate reports)',
    5: 'Clean reference format (clean references)',
    6: 'Markdown → HTML (convert)',
}


def convert_reports_to_html(parsed_dir: Path):
    """Iterate over the parsed directory, converting each parsed_report.md to HTML."""
    for case_dir in parsed_dir.glob('*'):
        md_path = case_dir / 'parsed_report.md'
        html_path = case_dir / 'parsed_report.html'
        if not md_path.exists():
            continue
        content = md_path.read_text(encoding='utf-8')
        try:
            html = markdown_to_html(content, case_dir.name)
            html_path.write_text(html, encoding='utf-8')
        except Exception as e:
            logger.error(f'HTML conversion failed for {case_dir.name}: {e}')


def orchestrate(
    df_path: str,
    workdir: Path,
    num_semaphore: int = 8,
    top_n: int | None = None,
    fdr_threshold: float = 1.0,
    start_stage: int = 1,
):
    workdir.mkdir(parents=True, exist_ok=True)

    output_dir = workdir / 'output'
    report_output_dir = workdir / 'report_output'
    parsed_output_dir = workdir / 'parsed_output'
    scored_csv = workdir / 'output_df.csv'
    filtered_csv = workdir / 'output_df_fdr.csv'

    for stage_num in range(start_stage, 7):
        logger.info(f'=== Stage {stage_num}/6: {STAGES[stage_num]} ===')

        if stage_num == 1:
            asyncio.run(run_pipeline(num_semaphore, df_path=df_path, workdir=workdir))

        elif stage_num == 2:
            df = rank_result(str(output_dir), df_path)
            df.to_csv(str(scored_csv), index=False)
            logger.info(f'Scoring complete, {len(df)} cases total')
            df = filter_fdr(df, fdr_threshold)
            df.to_csv(str(filtered_csv), index=False)
            logger.info(f'{len(df)} cases remaining after FDR filtering (FDR < {fdr_threshold})')

        elif stage_num == 3:
            effective_top_n = top_n or len(pd.read_csv(str(filtered_csv)))
            select_cases(
                csv_path=str(filtered_csv),
                from_dir=str(output_dir),
                to_dir=str(report_output_dir),
                top_n=effective_top_n,
            )
            logger.info(f'Selected {effective_top_n} cases to {report_output_dir}')

        elif stage_num == 4:
            effective_top_n = top_n or len(pd.read_csv(str(filtered_csv)))
            asyncio.run(process_top(
                report_output_dir, str(filtered_csv),
                top_n=effective_top_n, workdir=workdir,
            ))

        elif stage_num == 5:
            parsed_output_dir.mkdir(parents=True, exist_ok=True)
            clean_reports(str(report_output_dir), parsed_output_dir)
            logger.info(f'Reference cleaning complete, output to {parsed_output_dir}')

        elif stage_num == 6:
            convert_reports_to_html(parsed_output_dir)
            logger.info('HTML conversion complete')

    logger.info('=== All stages finished ===')


def parse_args():
    default_workdir = Path(__file__).parent.parent

    parser = argparse.ArgumentParser(
        description='ClockBase E2E Pipeline Orchestrator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Stage descriptions:
  1  Batch processing (pipeline)     Fetch metadata, generate queries, retrieve docs, analyze relevance
  2  Scoring & FDR filtering         Aggregate final_score, filter by FDR threshold
  3  Select Top-N cases              Copy top-ranked cases to report directory
  4  Generate reports                Generate academic reports for selected cases
  5  Clean reference format          Fix reference numbering and formatting in Markdown
  6  Markdown → HTML                Convert to browsable HTML files

Examples:
  python orchestrate.py --df-path ../data/main.csv
  python orchestrate.py --df-path ../data/main.csv --workdir /tmp/clockbase_run
  python orchestrate.py --df-path ../data/main.csv --start-stage 3 --top-n 500
        """,
    )

    parser.add_argument('--df-path', type=str, required=True,
                        help='Path to the main data CSV file')
    parser.add_argument('--workdir', type=str,
                        default=str(default_workdir),
                        help=f'Working directory for all intermediate artifacts (default: {default_workdir})')
    parser.add_argument('--num-semaphore', type=int, default=8,
                        help='Concurrency limit (default: 8)')
    parser.add_argument('--top-n', type=int, default=None,
                        help='Select the top N ranked cases (default: no limit, process all)')
    parser.add_argument('--fdr-threshold', type=float, default=1.0,
                        help='FDR filtering threshold (default: 1.0, i.e. no filtering)')
    parser.add_argument('--start-stage', type=int, default=1,
                        choices=range(1, 7),
                        help='Stage to start execution from, 1-6 (default: 1)')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    orchestrate(
        df_path=args.df_path,
        workdir=Path(args.workdir),
        num_semaphore=args.num_semaphore,
        top_n=args.top_n,
        fdr_threshold=args.fdr_threshold,
        start_stage=args.start_stage,
    )
