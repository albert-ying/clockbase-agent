"""Regenerate the released prompt mirrors without importing model clients."""
import argparse
import ast
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
PROMPTS = REPO / 'prompts'
MINIPAPER = 'minipaper_fork/src/nlp_process/prompt_template.py'
HUMAN = 'human_prompts_for_minipaper.py'
POST_AUDIT = 'prompts/post_audit_templates.py'
SOURCES = {
    'coding_agent_system.txt': ('src/coding_agent.py', 'DEFAULT_SYSTEM_PROMPT'),
    'task_intention_human_rna.txt': ('src/parallel_hs_rna.py', 'TASK_TEMPLATE'),
    'task_intention_mouse_rna.txt': ('src/parallel_mm_rna.py', 'TASK_TEMPLATE'),
    'task_intention_mouse_rna_all_samples_local.txt': ('src/parallel_mm_rna.py', 'TASK_TEMPLATE_ALL_SAMPLES_LOCAL'),
    'task_intention_mouse_rna_all_samples.txt': ('src/parallel_mm_rna.py', 'TASK_TEMPLATE_ALL_SAMPLES'),
    'task_intention_human_dnam.txt': ('src/parallel_hs_dnam.py', 'TASK_TEMPLATE'),
    'rag_query_generator.txt': (MINIPAPER, 'GENERATE_QUERY'),
    'rag_relevance_classifier.txt': (MINIPAPER, 'HAS_RELEVANCE'),
    'biological_interpreter_mouse_rna.txt': (MINIPAPER, 'IDENTIFY_RELEVANCE'),
    'biological_interpreter_human_rna.txt': (HUMAN, 'IDENTIFY_RELEVANCE_HUMAN_RNA'),
    'biological_interpreter_human_dnam.txt': (HUMAN, 'IDENTIFY_RELEVANCE_HUMAN_DNAM'),
    'generate_paper_mouse_rna.txt': (MINIPAPER, 'GENERATE_PAPER_2'),
    'generate_paper_human_rna.txt': (HUMAN, 'GENERATE_PAPER_2_HUMAN_RNA'),
    'generate_paper_human_dnam.txt': (HUMAN, 'GENERATE_PAPER_2_HUMAN_DNAM'),
    'post_audit_opus47.txt': (POST_AUDIT, 'POST_AUDIT_OPUS47'),
    'case_curator_opus47.txt': (POST_AUDIT, 'CASE_CURATOR_OPUS47'),
    # The verifier's model stages (each prompt is sent to two models; only agreement is accepted).
    'verifier_design_reconstruction.txt': ('audit_system/design_reconstruction_2026_09_28.py', 'PROMPT'),
    'verifier_label_adjudication.txt': ('audit_system/label_adjudication_2026_09_27.py', 'PROMPT'),
    'verifier_intervention_classification.txt': ('audit_system/intervention_classification_2026_09_29.py', 'PROMPT'),
    'verifier_relevance_rubric.txt': ('audit_system/intervention_classification_2026_09_29.py', 'RUBRIC'),
    'verifier_design_audit.txt': ('audit_system/design_audit_2026_09_30.py', 'PROMPT'),
}
DNAM_PROVENANCE = (
    'recovered from successful_traces_human_dnam.pkl, GSE211704, result[0][0]; '
    'removed the Human Message wrapper and substituted study_id, study_design, data_dir, RESULTS_DIR '
    '(all 2296 traces agree)'
)
RNA_PROVENANCE = {
    'task_intention_human_rna.txt': 'successful_traces_human_rna.pkl; SONNET37; all 9546 traces',
    'task_intention_mouse_rna.txt': 'successful_traces_mouse.pkl; base; 11385 traces; '
                                    '4c6aa975e::src/parallel_prompt3.py',
    'task_intention_mouse_rna_all_samples_local.txt': 'successful_traces_mouse.pkl; all_samples_local; '
                                                      '3578 traces; 3faf1f266::src/parallel_prompt3.py',
    'task_intention_mouse_rna_all_samples.txt': 'successful_traces_mouse.pkl; all_samples; 417 traces',
}


def source_body(source, symbol):
    """Read a literal constant or a driver's f-string with its named placeholders."""
    tree = ast.parse((REPO / source).read_text(encoding='utf-8'), filename=source)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == symbol
                                                for t in node.targets):
            value = ast.literal_eval(node.value)
            if not isinstance(value, str):
                raise ValueError(f'{source}::{symbol} is not a string')
            return value
        if isinstance(node, ast.FunctionDef) and node.name == symbol:
            returns = [n for n in ast.walk(node) if isinstance(n, ast.Return)]
            if len(returns) != 1 or not isinstance(returns[0].value, ast.JoinedStr):
                raise ValueError(f'{source}::{symbol} must return one f-string')
            parts = []
            for part in returns[0].value.values:
                if isinstance(part, ast.Constant) and isinstance(part.value, str):
                    parts.append(part.value)
                elif (isinstance(part, ast.FormattedValue) and isinstance(part.value, ast.Name)
                      and part.conversion == -1 and part.format_spec is None):
                    parts.append('{' + part.value.id + '}')
                else:
                    raise ValueError(f'Unsupported placeholder in {source}::{symbol}')
            return ''.join(parts)
    raise ValueError(f'Prompt source not found: {source}::{symbol}')


def released_prompt(name):
    """One provenance line followed by the source body, including its whitespace."""
    source, symbol = SOURCES[name]
    provenance = f'# Source: {source}::{symbol}'
    if name == 'task_intention_human_dnam.txt':
        provenance += '; ' + DNAM_PROVENANCE
    if name in RNA_PROVENANCE:
        provenance += '; recovered from ' + RNA_PROVENANCE[name]
    return provenance + '\n' + source_body(source, symbol)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Report differences without writing files')
    args = parser.parse_args(argv)
    differences = []
    try:
        extras = {path.name for path in PROMPTS.glob('*.txt')} - set(SOURCES)
        if extras:
            raise ValueError('Prompt files without sources: ' + ', '.join(sorted(extras)))
        for name in SOURCES:
            path = PROMPTS / name
            expected = released_prompt(name)
            if args.check:
                if not path.exists() or path.read_bytes() != expected.encode('utf-8'):
                    differences.append(name)
            else:
                path.write_bytes(expected.encode('utf-8'))
    except (OSError, ValueError, SyntaxError) as exc:
        print(f'Cannot build prompt mirrors: {exc}', file=sys.stderr)
        return 1
    if differences:
        print('Prompt differences: ' + ', '.join(differences), file=sys.stderr)
        return 1
    print(f'{len(SOURCES)} prompt mirrors: ' + ('no differences' if args.check else 'regenerated'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
