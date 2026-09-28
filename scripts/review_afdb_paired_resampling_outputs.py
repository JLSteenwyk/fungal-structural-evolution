#!/usr/bin/env python3
"""Read back completed interval tables and census native fit warnings."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time
import numpy as np
import pandas as pd


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def category(line):
    if re.fullmatch(r'WARNING: \d+ near-zero internal branches \(<[^)]+\) should be treated with caution',line):
        return 'near_zero_internal_branches'
    if re.fullmatch(r'WARNING: States\(s\) .+ rarely appear in alignment and may cause numerical problems',line):
        return 'rare_alignment_states'
    return 'other_warning_requires_review'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    started=time.perf_counter();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    for path,h in plan['pins'].items():assert sha(path)==h,path
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show','fungal-recovered-afdb-paired-resampling-20260926.service','-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),state
    control=Path(plan['control']);audit=Path(plan['audit'])
    terminal=json.loads((control/'receipt.json').read_text());ar=json.loads((audit/'receipt.json').read_text())
    assert terminal['plan_sha256']==ph and terminal['audit_receipt_sha256']==sha(audit/'receipt.json')
    assert terminal['status']=='complete_collection_paired_resampling_and_audit'
    assert ar['status']=='complete_paired_resampling_audit'
    for name,h in ar['artifacts'].items():assert sha(audit/name)==h,name
    arguments=plan['stages'][0]['arguments'];root=Path(arguments[arguments.index('--output')+1])
    assert sha(root/'receipt.json')==ar['resampling_receipt_sha256']
    rr=json.loads((root/'receipt.json').read_text());cfg=json.loads((root/'config.json').read_text())
    assert rr['config_sha256']==sha(root/'config.json')
    batches=pd.read_csv(audit/'batch_summary.tsv',sep='\t');table=pd.read_csv(audit/'conditional_intervals.tsv',sep='\t')
    assert len(batches)==ar['batches']==375 and batches.marker.nunique()==125
    assert not batches.duplicated(['marker','block_length']).any()
    assert set(batches.block_length)=={1,10,30} and (batches.draws==200).all()
    assert (batches.valid_draws+batches.unestimable_draws==batches.draws).all()
    assert int(batches.draws.sum())==ar['attempted_paired_draws']==75000
    assert int(batches.unestimable_draws.sum())==ar['unestimable_draws']
    assert len(table)==ar['branch_interval_rows'] and not table.duplicated(['marker','block_length','split_taxa','fit']).any()
    assert set(table.fit)=={'aa','3di'}
    merged=table.merge(batches,on=['marker','block_length'],suffixes=('','_batch'),validate='many_to_one')
    assert len(merged)==len(table) and (merged.valid_draws==merged.valid_draws_batch).all()
    assert (merged.planned_draws==merged.draws).all()
    assert table.interval_status.eq('conditional_percentile_sensitivity').all()
    numeric=['point_estimate','percentile_2_5','median','percentile_97_5','fraction_le_1e_minus5','sampling_sd','paired_sampling_covariance']
    assert np.isfinite(table[numeric].to_numpy()).all()
    assert (table.percentile_2_5<=table['median']).all() and (table['median']<=table.percentile_97_5).all()
    assert (table[numeric[:-1]]>=0).all().all() and (table.fraction_le_1e_minus5<=1).all()
    for key,frame in table.groupby(['marker','block_length']):
        row=batches.set_index(['marker','block_length']).loc[key]
        assert len(frame)==2*row.branches
        paired=frame.pivot(index='split_taxa',columns='fit',values='paired_sampling_covariance')
        np.testing.assert_array_equal(paired.aa,paired['3di'])
    summaries=[]
    for (block,fit),frame in table.groupby(['block_length','fit']):
        summaries.append(dict(block_length=int(block),fit=fit,rows=len(frame),
            rows_with_at_least_half_draws_near_zero=int((frame.fraction_le_1e_minus5>=.5).sum()),
            median_conditional_interval_width=float((frame.percentile_97_5-frame.percentile_2_5).median())))
    args.output.mkdir(parents=True,exist_ok=False)
    counts=Counter();by_fit=Counter();examples={};invalid=Counter();warned=fits=receipts=0
    with gzip.open(args.output/'fit_warning_inventory.jsonl.gz','wt') as handle:
        for batch in rr['results']:
            folder=root/batch['marker']/('block-'+str(batch['block_length']))
            assert sha(folder/'receipt.json')==batch['receipt_sha256']
            br=json.loads((folder/'receipt.json').read_text())
            assert set(br['receipt_hashes'])=={f'replicate-{i:04d}/receipt.json' for i in range(cfg['replicates'])}
            for name,h in br['receipt_hashes'].items():
                path=folder/name;assert sha(path)==h
                r=json.loads(path.read_text());receipts+=1
                if r['status']=='unestimable_draw':
                    assert not r['fits'];invalid[r['reason']]+=1;continue
                assert r['status']=='completed_draw'
                for fit in ['aa','3di']:
                    warnings=set();hashes={}
                    for suffix in ['.iqtree','.log']:
                        p=path.parent/(fit+suffix);actual=sha(p);assert actual==r['artifacts'][p.name]
                        hashes[p.name]=actual
                        warnings.update(line.strip() for line in p.read_text().splitlines() if 'WARNING:' in line)
                    kinds={category(line) for line in warnings}
                    for kind in kinds:
                        counts[kind]+=1;by_fit[fit+'|'+kind]+=1
                        examples.setdefault(kind,next(line for line in sorted(warnings) if category(line)==kind))
                    fits+=1;warned+=bool(warnings)
                    handle.write(json.dumps(dict(marker=batch['marker'],block_length=batch['block_length'],replicate=path.parent.name,fit=fit,warnings=sorted(warnings),categories=sorted(kinds),source_receipt_sha256=h,source_artifacts=hashes))+'\n')
            print('Reviewed native warnings',receipts,'/ 75000',flush=True)
    assert fits==ar['validated_fits'] and warned==ar['fits_with_warnings']
    assert receipts==75000 and sum(invalid.values())==ar['unestimable_draws']
    result=dict(status='completed_interval_readback_and_native_warning_census',source_plan_sha256=ph,
        source_audit_receipt_sha256=sha(audit/'receipt.json'),fits=fits,fits_with_warnings=warned,
        unestimable_draw_reasons=dict(invalid),warning_categories=dict(counts),warnings_by_fit=dict(by_fit),warning_examples=examples,
        interval_summaries=summaries,minimum_valid_draws=int(batches.valid_draws.min()),
        batches_with_unestimable_draws=batches.loc[batches.unestimable_draws.gt(0)].to_dict('records'),
        wall_seconds=time.perf_counter()-started,script_sha256=sha(__file__),
        artifacts={'fit_warning_inventory.jsonl.gz':sha(args.output/'fit_warning_inventory.jsonl.gz')},
        scope='Complete interval-table consistency checks and receipt-bound native warning census. Does not rerun phylogenetic fits or independent quantile reconstruction. Percentile intervals condition on fitted sequence topology and model; near-zero and rare-state warnings remain substantive limitations.')
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
