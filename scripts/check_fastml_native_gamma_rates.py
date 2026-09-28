#!/usr/bin/env python3
"""Compare native versus SciPy category rates at all completed fitted alphas."""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha,write_json


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--readback',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    receipt=json.loads((args.readback/'receipt.json').read_text());assert not receipt['unresolved']
    records={}
    for key,entry in receipt['readbacks'].items():
        assert sha(entry['path'])==entry['sha256'];r=json.loads(Path(entry['path']).read_text())
        if r['status']=='no_coded_characters':continue
        assert r['status']=='integrity_and_independent_replay_complete_not_fit_qualification';records[key]=r
    assert len(records)==306
    library=Path('data/software_audits/fastml-3.11/source/FastML.v3.11/libs/phylogeny')
    archive=library/'libEvolTree.a';assert sha(archive)=='5052c19e72fc2d30353afb397180c9c604aa2407b48e1c0644296c04bbcf4372'
    for variant in ['precision-only','precision-cache-refresh']:
        assert sha(Path('data/software_audits/fastml-precision-cache-20260928-v1')/variant/'libs/phylogeny/libEvolTree.a')==sha(archive)
    args.output.mkdir(parents=True,exist_ok=False);binary=args.output/'export_rates'
    command=['g++','-O3','-DLOG','-I',str(library),'scripts/export_fastml_gamma_rates.cpp',str(archive),'-o',str(binary)]
    subprocess.run(command,check=True,capture_output=True)
    text=''.join(key+' '+repr(r['fitted_parameters']['alpha'])+'\n' for key,r in records.items())
    (args.output/'requests.tsv').write_text(text)
    output=subprocess.check_output([str(binary)],input=text,text=True);(args.output/'native_rates.tsv').write_text(output)
    rows=[];seen=set()
    for line in output.splitlines():
        key,*fields=line.split('\t');assert key in records and key not in seen;seen.add(key)
        native=np.array(fields,dtype=float);reference=np.array(records[key]['gamma_category_rates'])
        assert native.shape==(4,) and np.isfinite(native).all() and (native>0).all()
        rows.append(dict(id=key,alpha=records[key]['fitted_parameters']['alpha'],maximum_rate_difference=float(np.max(abs(native-reference))),
                         native_mean_rate=float(native.mean()),scipy_mean_rate=float(reference.mean()),native_rates=native.tolist(),scipy_rates=reference.tolist()))
    assert seen==set(records)
    (args.output/'rate_comparison.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    result=dict(status='full_fitted_alpha_native_rate_comparison_complete',fits=len(rows),
        maximum_absolute_category_rate_difference=max(r['maximum_rate_difference'] for r in rows),
        maximum_absolute_native_mean_rate_error=max(abs(r['native_mean_rate']-1) for r in rows),
        fits_rate_difference_over_1e_minus8=sum(r['maximum_rate_difference']>1e-8 for r in rows),
        build_command=command,compiler=subprocess.check_output(['g++','--version'],text=True),
        pins={str(p):sha(p) for p in [archive,Path('scripts/export_fastml_gamma_rates.cpp'),Path(__file__),args.readback/'receipt.json']},
        artifacts={p.name:sha(p) for p in args.output.iterdir()},
        scope='Uses the identical frozen FastML static library shared by both variants. Measures category discretization differences at every reported fitted alpha. Does not itself show how much of the likelihood or marginal discrepancy these rates explain; native-rate replay remains required.')
    write_json(args.output/'receipt.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
