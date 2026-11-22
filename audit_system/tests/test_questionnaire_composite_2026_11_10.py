"""Verify the questionnaire record against the saved comparison source."""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from audit_system.questionnaire_composite_2026_11_10 import (
    CAP, DIRECTIONS, KEY, P_FLOOR, SCOPE, questionnaire_table)


@pytest.fixture
def source():
    verifier = Path(os.environ.get('CLOCKBASE_VERIFIER_ROOT') or
                    Path(__file__).resolve().parents[3] / 'data/clockbase-verifier-v3-2026-09-27')
    path = verifier / 'final-2026-10-10/comparison-ranking/ranked_comparisons_all.csv'
    if not path.exists():
        pytest.fail(f'Cannot proceed - required data file not found: {path}')
    return pd.read_csv(path, low_memory=False)


def test_formula_order_and_scope(source):
    out, m, M = questionnaire_table(source)
    scored = source.condition_category.isin(SCOPE) & source.final_score.notna()
    assert m == source.loc[scored, 'final_score'].min()
    assert M == source.loc[scored, 'final_score'].max()
    pd.testing.assert_frame_equal(out[KEY], source[KEY])
    pd.testing.assert_series_equal(out.questionnaire_score, source.final_score.where(scored), check_names=False)
    expected = np.minimum(-np.log10(np.maximum(source.loc[scored, 'v3_p'], P_FLOOR)), CAP)
    expected *= (source.loc[scored, 'final_score'] - m) / (M - m)
    np.testing.assert_allclose(out.loc[scored, 'questionnaire_composite'], expected)
    assert out.loc[~scored, ['questionnaire_score', 'questionnaire_composite', 'questionnaire_composite_rank']].isna().all().all()
    assert not any('judgement' in c or 'threshold' in c or 'relevance' in c for c in out)
    for direction in DIRECTIONS:
        ordered = out[scored & out.direction_v3.eq(direction)].sort_values(
            ['questionnaire_composite', 'statistical_rank'], ascending=[False, True], kind='stable')
        assert ordered.questionnaire_composite_rank.tolist() == list(range(1, len(ordered) + 1))
    # The source frame and every statistical position remain unchanged.
    pd.testing.assert_series_equal(out.statistical_rank, source.groupby('direction_v3').cumcount() + 1, check_names=False)


def test_zero_range_and_invalid_inputs_fail(source):
    same = source.copy()
    same['final_score'] = source.final_score.dropna().iloc[0]
    with pytest.raises(ValueError, match='range'):
        questionnaire_table(same)
    invalid = source.copy()
    scored = invalid.condition_category.isin(SCOPE) & invalid.final_score.notna()
    invalid.loc[invalid.index[scored][0], 'v3_p'] = np.nan
    with pytest.raises(ValueError, match='P value'):
        questionnaire_table(invalid)
    with pytest.raises(ValueError, match='duplicate'):
        questionnaire_table(pd.concat([source, source]))
