#!/usr/bin/env python3
"""Draw the eight frozen conditional estimates and replace only ED6 in a new PDF."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import DictionaryObject, NameObject

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STYLE = ROOT / 'pipeline_assets/figure_style_2026_09_10/nature_figure_utils.py'
STYLE_SHA = '605d8efaaefe1ec76f991708d7ccb657250e50e318b21c633b9ed5a100cdb21c'
sys.path.insert(0, str(STYLE.parent))
from nature_figure_utils import (retro_style, figsize_nature, trim_axes,
    add_panel_labels, save_figure, figure_preflight, set_domain, BLUE, RUST)


def sha(path):
    with Path(path).open('rb') as h:
        return hashlib.file_digest(h, 'sha256').hexdigest()


def write(path, value):
    with path.open('x') as h:
        json.dump(value,h,indent=2,allow_nan=False)
        h.write('\n')


def figure(output):
    run = HERE / 'numerical-run-001'
    receipt = json.loads((run/'receipt.json').read_text())
    if receipt['status']!='PASS' or receipt['scientific_status']!='CANNOT VERIFY':
        raise ValueError('Expected completed conditional calculation')
    for name, digest in receipt['output_sha256'].items():
        if sha(run/name)!=digest:
            raise ValueError('Calculation artifact changed: '+name)
    if sha(STYLE)!=STYLE_SHA:
        raise ValueError('Bundled figure helper changed')
    with (run/'conditional_models.csv').open() as h:
        rows=list(csv.DictReader(h))
    plan=json.loads((HERE/'plan.json').read_text())
    keys={(r['clock'],r['formula']) for r in rows}
    if len(rows)!=8 or keys!={(c,f) for c in plan['clocks'] for f in plan['models']}:
        raise ValueError('Expected all four models for both clocks')
    retro_style('soft')
    fig, ax = plt.subplots(1,2,figsize=figsize_nature('double',aspect=.38),layout=None)
    fig.subplots_adjust(left=.16,right=.98,bottom=.24,top=.83,wspace=.62)
    labels=['Condition only','With batch','With pool','With batch + pool']
    for axis,clock,title,color in zip(ax,plan['clocks'],['Rodents clock','Mouse clock'],[BLUE['400'],RUST['400']]):
        subset=[next(r for r in rows if r['clock']==clock and r['formula']==f) for f in plan['models']]
        values=np.array([float(r['condition_coefficient']) for r in subset])
        low=np.array([float(r['ci_lower']) for r in subset])
        high=np.array([float(r['ci_upper']) for r in subset])
        if not np.isfinite(np.r_[values,low,high]).all():
            raise ValueError('Nonfinite confidence interval')
        axis.axvline(0,color='#999999',lw=.6,ls='--',zorder=0)
        axis.errorbar(values,np.arange(4),xerr=[values-low,high-values],fmt='o',
                      color=color,ms=4,elinewidth=.85,capsize=0)
        axis.set_yticks(np.arange(4),labels,fontsize=7)
        axis.set_title(title+r' ($n$ = 15 libraries)',fontsize=8,loc='left')
        axis.set_xlabel('Lesion-class coefficient\n(native clock units)',fontsize=8)
        axis.set_xticks([-.4,-.2,0,.2])
        axis.set_xlim(-.4,.25)
        axis.set_ylim(3.4,-.4)
    trim_axes(fig)
    for axis in ax:
        set_domain(axis,xlim=(-.4,.25),ylim=(3.4,-.4))
    axes={'a':ax[0],'b':ax[1]}
    add_panel_labels(fig,axes)
    footer='Saved predictions; conditional 95% intervals; donor independence and clock applicability unresolved'
    preflight=figure_preflight(fig,axes)
    write(output/'preflight.json',preflight)
    if not preflight['passed']:
        raise ValueError('Figure preflight failed; inspect preflight.json')
    save_figure(fig,output/'ED6_covariate_sensitivity.pdf',axes_dict=axes,
                footer=footer,review_png=output/'ED6_covariate_sensitivity.png',dpi=300)
    fig.savefig(output/'ED6_covariate_sensitivity.svg',bbox_inches='tight')
    plt.close(fig)
    return {'calculation_receipt_sha256':sha(run/'receipt.json'),
            'calculation_table_sha256':sha(run/'conditional_models.csv'),
            'figure_style_sha256':sha(STYLE),'estimate_count':len(rows),
            'scope':'Saved-prediction conditional OLS estimates and 95% t intervals; no independent animal or clock validation'}


def assembly(output, previous_path, new_path):
    if new_path.exists():
        raise FileExistsError('Refusing to replace an existing assembly')
    old=PdfReader(previous_path)
    replacement=PdfReader(output/'ED6_covariate_sensitivity.pdf').pages[0]
    target=11
    if len(old.pages)!=26 or ' '.join(old.pages[target].extract_text().split())!='Extended Data Figure 6':
        raise ValueError('Assembly identity or ED6 title mismatch')
    writer=PdfWriter()
    original_contents={i:p.get_contents().get_data() for i,p in enumerate(old.pages)}
    for i,page in enumerate(old.pages):
        if i==target:
            resources=DictionaryObject(dict(page['/Resources']))
            xobjects=DictionaryObject(dict(resources['/XObject']))
            resources[NameObject('/XObject')]=xobjects
            page[NameObject('/Resources')]=resources
            names={name for name in xobjects if xobjects[name].get_object().get('/Subtype')=='/Image'}
            content=page.get_contents()
            removed=[v[0] for v,op in content.operations if op==b'Do' and v[0] in names]
            if len(removed)!=1:
                raise ValueError('Expected one historical ED6 raster')
            content.operations=[(v,op) for v,op in content.operations if not(op==b'Do' and v[0] in names)]
            page.replace_contents(content)
            for name in removed:
                del xobjects[name]
            width,height=float(page.mediabox.width),float(page.mediabox.height)
            sw,sh=float(replacement.mediabox.width),float(replacement.mediabox.height)
            scale=min(1.,(width-72)/sw,(height-100)/sh)
            if scale<.97:
                raise ValueError('Replacement would shrink below readable size')
            dx,dy=(width-sw*scale)/2,height-80-sh*scale
            page.merge_transformed_page(replacement,Transformation().scale(scale).translate(dx,dy))
        writer.add_page(page)
    with new_path.open('xb') as h:
        writer.write(h)
    new=PdfReader(new_path)
    checks={f'page_{i+1}_content_preserved':new.pages[i].get_contents().get_data()==original_contents[i]
            for i in range(26) if i!=target}
    checks['page_count_preserved']=len(new.pages)==26
    checks['ed6_title_preserved']='Extended Data Figure 6' in new.pages[target].extract_text()
    checks['scope_in_figure']='clock applicability unresolved' in new.pages[target].extract_text()
    if not all(checks.values()):
        raise ValueError('Assembly preservation checks failed')
    subprocess.run(['pdftoppm','-f','12','-l','12','-singlefile','-png','-r','130',str(new_path),str(output/'ED6_assembled_page')],check=True,capture_output=True)
    # Compare every untouched rendered page as well as its PDF content.
    def render_bytes(path, index):
        return subprocess.run(['pdftoppm','-f',str(index+1),'-l',str(index+1),
            '-singlefile','-r','72',str(path)],check=True,capture_output=True).stdout
    unchanged=[i+1 for i in range(26) if i!=target
               and render_bytes(previous_path,i)==render_bytes(new_path,i)]
    if len(unchanged)!=25:
        raise ValueError('An unrelated page render changed')
    return {'previous_sha256':sha(previous_path),'assembly_sha256':sha(new_path),
            'assembly_path':str(new_path),'replaced_page':12,'unchanged_rendered_pages':unchanged,
            'checks':checks,'scale':scale,'translation':[dx,dy]}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True)
    p.add_argument('--assembly',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():
        p.error('Use a new output directory')
    args.output.mkdir(parents=True)
    result={'started_utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),
            'command':sys.argv,'matplotlib':matplotlib.__version__,'python':sys.version}
    try:
        result.update(figure(args.output))
        result.update(assembly(args.output,args.previous,args.assembly))
        result['status']='PASS'
    except Exception as error:
        result.update(status='FAIL',error=str(error),error_type=type(error).__name__)
    result['output_sha256']={p.name:sha(p) for p in args.output.iterdir() if p.is_file()}
    write(args.output/'figure_receipt.json',result)
    print(json.dumps({'status':result['status'],'error':result.get('error')}))
    if result['status']!='PASS':
        raise SystemExit(1)


if __name__=='__main__':
    main()
