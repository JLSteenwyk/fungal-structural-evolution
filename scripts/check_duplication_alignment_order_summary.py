#!/usr/bin/env python3
"""Known endpoint reversal, failed-order retention, and complete-grid fixture."""
import csv,json,tempfile
from pathlib import Path
from summarize_duplication_alignment_orders import summarize_pair,summarize,sha

row=dict(aligned_length='10',rmsd_recomputed='1.2',sequence_identity_exact='.5',joint_plddt70_fraction='.8',tm_left_native='.6',tm_right_native='.8',coverage_left='.5',coverage_right='1')
reverse=dict(row,tm_left_native='.7',tm_right_native='.55',coverage_left='1',coverage_right='.5')
x=summarize_pair('x','full',{0:'aligned',1:'aligned'},{0:row,1:reverse})
assert x['order1_tm_a']==.55 and x['order1_tm_b']==.7
assert abs(x['tm_a_order_absolute_difference']-.05)<1e-14 and x['coverage_a_order_absolute_difference']==0
x=summarize_pair('x','full',{0:'aligned',1:'timeout'},{0:row})
assert x['order1_tm_a']=='' and x['tm_a_order_absolute_difference']=='' and x['order_summary_status']=='one_order_aligned'
with tempfile.TemporaryDirectory() as td:
    root=Path(td);prod=root/'producer';audit=root/'audit';prod.mkdir();audit.mkdir();pair='a'*64
    def write(p,data):p.write_text(json.dumps(data)+'\n')
    def table(p,rows):
        with p.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    dispositions=[dict(path=f'pairs/aa/{pair}-{mask}-{order}.json',sha256='fixture',status='aligned' if mask=='full' else 'input_unavailable') for mask in ['full','plddt70'] for order in [0,1]]
    table(prod/'checkpoint_manifest.tsv',dispositions)
    table(audit/'numeric_readback.tsv',[dict(pair_key=pair,mask='full',order=o,**r) for o,r in [(0,row),(1,reverse)]])
    counts={'full:aligned':2,'plddt70:input_unavailable':2}
    receipt=dict(status='complete_duplication_alignment_dispositions_pending_readback',distinct_model_pairs=1,directed_dispositions=4,counts=counts,artifacts={'checkpoint_manifest.tsv':sha(prod/'checkpoint_manifest.tsv')})
    write(prod/'receipt.json',receipt)
    ar=dict(status='passed_full_duplication_alignment_mapping_rmsd_identity_readback',producer_receipt_sha256=sha(prod/'receipt.json'),mode='primary',directed_dispositions=4,numerically_checked_alignments=2,counts=counts,artifacts={'numeric_readback.tsv':sha(audit/'numeric_readback.tsv')})
    write(audit/'receipt.json',ar)
    result=summarize(prod,audit,root/'summary');assert result['pair_mask_rows']==2
    table(prod/'checkpoint_manifest.tsv',dispositions[:-1]);receipt['artifacts']['checkpoint_manifest.tsv']=sha(prod/'checkpoint_manifest.tsv');write(prod/'receipt.json',receipt);ar['producer_receipt_sha256']=sha(prod/'receipt.json');write(audit/'receipt.json',ar)
    try:summarize(prod,audit,root/'bad')
    except ValueError:pass
    else:raise AssertionError('Missing disposition accepted')
print('Endpoint reversal, explicit failed order, complete grid and rehashed missing-disposition rejection passed.')
