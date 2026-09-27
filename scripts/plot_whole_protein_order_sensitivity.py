#!/usr/bin/env python3
"""Describe primary/reference order sensitivity with explicit conditional denominators."""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from screen_duplication_domain_alignment_coverage import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        assert sha(args.plan)==ph
        for path,digest in plan['pins'].items():assert sha(path)==digest,path
    verify();dep=plan['producer']
    while True:
        try:
            proc=psutil.Process(dep['pid'])
            if proc.create_time()!=dep['created'] or proc.status()==psutil.STATUS_ZOMBIE:break
            assert proc.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        print('waiting_for_full_primary_mapping_readback',dep['pid'],flush=True);time.sleep(30)
    terminal=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert terminal==dict(ActiveState='inactive',Result='success',ExecMainStatus='0'),terminal
    verify();datasets={};summaries=[];quantiles=[];hashes={}
    for label,spec in plan['sources'].items():
        root=Path(spec['root']);receipt=json.loads((root/'receipt.json').read_text());audit=json.loads(Path(spec['audit']).read_text())
        assert audit['status']==f'passed_full_{label}_order_sensitivity_readback'
        assert audit['source_receipt_sha256']==sha(root/'receipt.json')
        assert audit['pair_mask_rows']==receipt['pair_mask_rows'] and audit['mapping_counts']==receipt['mapping_counts']
        for name,h in receipt['artifacts'].items():assert sha(root/name)==h
        for path in [root/'receipt.json',Path(spec['audit'])]:hashes[str(path)]=sha(path)
        data=pd.read_csv(root/'pair_mask_sensitivity.tsv',sep='\t',float_precision='round_trip')
        assert len(data)==audit['pair_mask_rows'] and not data.duplicated(['pair_key','mask']).any()
        assert data.in_both_mask_cohort.isin([True,False]).all()
        common=data[data.in_both_mask_cohort].copy()
        assert set(common[common['mask'].eq('full')].pair_key)==set(common[common['mask'].eq('plddt70')].pair_key)
        assert len(common)==2*receipt['common_mask_pairs']
        datasets[label]=common
        for cohort,sub in [('full_all',data[data['mask'].eq('full')]),('full_common',common[common['mask'].eq('full')]),('plddt70_common',common[common['mask'].eq('plddt70')])]:
            counts=sub.mapping_status.value_counts().to_dict();assert counts==receipt['mapping_counts'][cohort]
            identical=counts.get('identical',0);same=counts.get('different_same_count',0);different=counts.get('different_count',0)
            assert identical+same+different==len(sub)
            summaries.append(dict(group=label,cohort=cohort,pairs=len(sub),identical=identical,different_same_count=same,different_count=different,changed_pairs=same+different,changed_fraction=(same+different)/len(sub)))
        q=pd.read_csv(root/'quantiles.tsv',sep='\t',float_precision='round_trip');q.insert(0,'group',label);quantiles.append(q)
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    summary=pd.DataFrame(summaries);summary.to_csv(out/'mapping_summary.tsv',sep='\t',index=False)
    joined=pd.concat(quantiles,ignore_index=True);joined.to_csv(out/'all_metric_quantiles.tsv',sep='\t',index=False)
    for name,expected in [('mapping_summary.tsv',summary),('all_metric_quantiles.tsv',joined)]:
        pd.testing.assert_frame_equal(pd.read_csv(out/name,sep='\t',float_precision='round_trip'),expected)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(1,3,figsize=(13,4.4),gridspec_kw={'width_ratios':[1.2,1,1]})
    labels=[];same=[];diff=[];notes=[]
    for group in ['primary','reference']:
        for cohort in ['full_common','plddt70_common']:
            row=summary[(summary.group==group)&(summary.cohort==cohort)].iloc[0]
            labels.append(group.title()+'\n'+('Full' if cohort=='full_common' else 'pLDDT ≥70'))
            same.append(100*row.different_same_count/row.pairs);diff.append(100*row.different_count/row.pairs)
            notes.append(f'{row.changed_pairs:,}/{row.pairs:,}')
    xx=np.arange(4);axes[0].bar(xx,same,color='#0072B2',label='Same count, different mapping');axes[0].bar(xx,diff,bottom=same,color='#E69F00',label='Different count')
    for x,y,z,txt in zip(xx,same,diff,notes):axes[0].text(x,y+z+.009,txt,ha='center',va='bottom',fontsize=7.5)
    axes[0].set_xticks(xx,labels);axes[0].set_ylim(0,max(np.array(same)+diff)*1.35)
    axes[0].set_ylabel('Pairs with changed residue mapping (%)');axes[0].set_title('A  All common-mask pairs',loc='left');axes[0].legend(loc='upper left',fontsize=7,frameon=False)
    plot_counts={}
    for ax,mask,panel in zip(axes[1:],['full','plddt70'],['B','C']):
        for label,color in [('primary','#0072B2'),('reference','#D55E00')]:
            sub=datasets[label];sub=sub[sub['mask'].eq(mask)&sub.mapping_status.ne('identical')]
            values=np.sort(sub.rmsd_recomputed_absolute_difference.to_numpy());assert len(values)>0 and np.isfinite(values).all()
            # Every changed-mapping observation is retained, including zero differences.
            ax.step(np.r_[values[0],values],np.r_[0,np.arange(1,len(values)+1)/len(values)],where='post',color=color,label=f'{label.title()} (n={len(values):,})')
            plot_counts[label+':'+mask]=len(values)
        ax.set_xscale('symlog',linthresh=1e-6);ax.set_xlim(left=0);ax.set_ylim(0,1.02)
        ax.set_xlabel('Absolute input-order RMSD difference (Å)');ax.set_ylabel('Cumulative fraction of changed mappings')
        ax.set_title(panel+'  Changed mappings: '+('full' if mask=='full' else 'pLDDT ≥70'),loc='left');ax.legend(frameon=False,fontsize=8)
    fig.suptitle('Whole-protein alignment order sensitivity',fontsize=13,y=.99)
    fig.text(.5,.015,'Pairs are held constant across masks within each group. B–C show changed mappings only. Descriptive checks; no independent-event or biological-effect claim.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.065,1,.94))
    for suffix in ['png','pdf','svg']:fig.savefig(out/('order_sensitivity.'+suffix),dpi=180)
    plt.close(fig);verify()
    for path,h in hashes.items():assert sha(path)==h
    result=dict(status='complete_audited_source_order_sensitivity_figure_pending_visual_review',plan_sha256=ph,source_hashes=hashes,summary_rows=len(summary),metric_quantile_rows=len(joined),conditional_ecdf_counts=plot_counts,prerequisite_terminal_state=terminal,artifacts={p.name:sha(p) for p in out.iterdir()},scope='Descriptive primary/reference comparison. Same pairs across masks within each group; groups are not matched to each other. Panel A retains all common-mask pairs; B/C explicitly condition on changed mappings. Repeated/phylogenetically dependent pairs are not independent observations. No confidence interval, test or biological effect estimated.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
