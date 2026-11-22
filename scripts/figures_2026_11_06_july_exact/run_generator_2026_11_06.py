#!/usr/bin/env python3
"""Run a generator with read-only recording before its unchanged PDF save."""
from pathlib import Path
import json
import os
import runpy
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = ROOT / 'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/july_exact'
sys.path.insert(0, str(WORK / 'dependencies'))
os.environ['MPLBACKEND'] = 'Agg'
import matplotlib
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pandas
import numpy
import scipy
import seaborn
import aquarel

sys.path.insert(0, 'pipeline_assets/figure_style_2026_09_10')
from nature_figure_utils import figure_preflight

original_savefig = Figure.savefig


def recorded_savefig(fig, fname, *args, **kwargs):
    if str(fname).endswith('.pdf'):
        dest = WORK / 'logs' / (Path(fname).stem + '_figure_record.json')
        fig.canvas.draw()
        record = {
            'generator':sys.argv[1], 'figsize_inches':list(fig.get_size_inches()),
            'versions':{'matplotlib':matplotlib.__version__, 'pandas':pandas.__version__,
                'numpy':numpy.__version__, 'scipy':scipy.__version__, 'seaborn':seaborn.__version__,
                'aquarel':'0.0.7'},
            'axes':[{'rect_fraction':list(ax.get_position().bounds),
                'xlim':list(ax.get_xlim()), 'ylim':list(ax.get_ylim()),
                'xlabel':ax.get_xlabel(), 'ylabel':ax.get_ylabel(),
                'ytick_labels':[{'text':t.get_text(), 'weight':t.get_fontweight()} for t in ax.get_yticklabels()],
                'bars':[{'width':p.get_width(), 'color':p.get_facecolor(), 'hatch':p.get_hatch()} for p in ax.patches if hasattr(p, 'get_width')],
                'text':[t.get_text() for t in ax.texts],
                'xticks':list(ax.get_xticks()), 'yticks':list(ax.get_yticks())} for ax in fig.axes],
            'figure_text':[t.get_text() for t in fig.texts],
        }
        try:
            # Findings are recorded, not used to restyle the reference.
            report = figure_preflight(fig, strict=False, verbose=False)
            record['preflight'] = report
        except Exception as exc:
            record['preflight_error'] = f'{type(exc).__name__}: {exc}'
        dest.write_text(json.dumps(record, indent=2,
            default=lambda x:x.item() if hasattr(x,'item') else str(x)) + '\n')
    return original_savefig(fig, fname, *args, **kwargs)


Figure.savefig = recorded_savefig
target = sys.argv[1]
try:
    if target in ['pathway','fda']:
        sys.argv = [str(HERE/'build_july_exact_2026_11_06.py'), target]
        runpy.run_path(str(HERE/'build_july_exact_2026_11_06.py'), run_name='__main__')
    else:
        runpy.run_path(str(HERE/target), run_name='__main__')
except Exception as exc:
    print(f'Generator failed: {target}: {exc}', file=sys.stderr)
    raise
