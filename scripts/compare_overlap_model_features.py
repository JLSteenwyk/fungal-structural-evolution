#!/usr/bin/env python3
"""Compare full native features for every distinct complete-sequence-matched overlap pair."""
import argparse,csv,json,hashlib
from collections import Counter
from pathlib import Path
import numpy as np
from compare_predictor_alphabets import compare_states
from prepare_paired_phylogenetic_inputs import write_table
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['context','context-readback','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();rp=a.context/'receipt.json';r=json.loads(rp.read_text());audit=json.loads(a.context_readback.read_text())
    if r['status']!='complete_paired_source_full_sequence_context' or audit['status']!='passed_paired_source_model_context_readback' or audit['producer_receipt_sha256']!=sha(rp):raise ValueError('Unverified model context')
    table=a.context/'source_model_context.tsv'
    if sha(table)!=r['artifacts'][table.name]:raise ValueError('Changed context table')
    pins={str(x):sha(x) for x in [rp,a.context_readback,table,Path('scripts/compare_predictor_alphabets.py')]}
    pairs={}
    for row in csv.DictReader(table.open(),delimiter='\t'):
        if row['full_sequence_status']!='identical_complete_encoded_sequence':raise ValueError('Different full sequence')
        key=tuple(row[k] for k in ['reference_model_id','reference_model_version','local_model_id','local_model_version'])
        descriptor={k:v for k,v in row.items() if k.endswith(('_encoding_path','_encoding_sha256','_sequence_sha256'))}
        if key in pairs and pairs[key]!=descriptor:raise ValueError('Inconsistent repeated model pair')
        pairs[key]=descriptor
    if len(pairs)!=r['distinct_model_pairs']:raise ValueError('Pair scope mismatch')
    rows=[];confusion=Counter()
    for key,desc in sorted(pairs.items()):
        features=[]
        for prefix in ['reference','local']:
            path=Path(desc[prefix+'_encoding_path']);digest=desc[prefix+'_encoding_sha256']
            if sha(path)!=digest:raise ValueError('Changed feature encoding')
            pins[str(path)]=digest
            with np.load(path,allow_pickle=False) as f:data={k:f[k].copy() for k in f.files}
            data['sequence']=str(data['sequence']);data['states_array']=np.asarray(list(str(data['states'])))
            if hashlib.sha256(data['sequence'].encode()).hexdigest()!=desc[prefix+'_sequence_sha256']:raise ValueError('Changed full sequence')
            features.append(data)
        for cutoff in [0,70,90]:
            for pae in [None,5,10,15]:
                row,mask=compare_states(*features,cutoff,pae)
                rows.append(dict(reference_model_id=key[0],reference_version=key[1],local_model_id=key[2],local_version=key[3],sequence_sha256=desc['reference_sequence_sha256'],**row))
                if cutoff==70 and pae==10 and row['status']=='compared':confusion.update(zip(features[0]['states_array'][mask],features[1]['states_array'][mask]))
    summary=[]
    for cutoff in [0,70,90]:
        for pae in ['unfiltered',5,10,15]:
            subset=[r for r in rows if r['plddt_cutoff']==cutoff and r['pae_cutoff']==pae];eligible=[r for r in subset if r['status']=='compared']
            record=dict(plddt_cutoff=cutoff,pae_cutoff=pae,model_pairs=len(subset),compared=len(eligible),coverage_excluded=len(subset)-len(eligible))
            for name in ['retained_residues','state_mismatches','partner_changes','same_partner_residues','same_partner_state_mismatches','changed_partner_residues','changed_partner_state_mismatches']:record[name]=sum(x[name] for x in eligible)
            record['pooled_state_mismatch_fraction']=record['state_mismatches']/record['retained_residues'] if record['retained_residues'] else ''
            record['median_state_mismatch_fraction']=float(np.median([x['state_mismatch_fraction'] for x in eligible])) if eligible else ''
            for context in ['same_partner','changed_partner']:
                record[context+'_pooled_mismatch_fraction']=record[context+'_state_mismatches']/record[context+'_residues'] if record[context+'_residues'] else ''
            summary.append(record)
    a.output.mkdir(parents=True,exist_ok=False)
    write_table(a.output/'model_pair_comparisons.tsv',rows);write_table(a.output/'threshold_summary.tsv',summary)
    alphabet='ACDEFGHIKLMNPQRSTVWY';write_table(a.output/'state_confusion_plddt70_pae10.tsv',[dict(reference_state=x,local_state=y,residues=confusion[x,y]) for x in alphabet for y in alphabet])
    for path,h in pins.items():
        if sha(path)!=h:raise ValueError('Source changed during comparison')
    result=dict(status='complete_overlap_model_feature_comparison',model_pairs=len(pairs),comparison_rows=len(rows),source_sha256=pins,script_sha256=sha(__file__),artifacts={x.name:sha(x) for x in a.output.glob('*.tsv')},scope='All distinct complete-sequence-matched model pairs in qualified cohort overlap. Full encoded protein features, not only marker projection. Twelve joint confidence alternatives; summary requires >=50 and >=half full protein retained. Native partner changes associated with state differences, not established causes. Dependent models/residues, changing threshold cohorts and acquisition bias preclude population/error/causal calibration.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
