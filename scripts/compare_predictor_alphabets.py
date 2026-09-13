#!/usr/bin/env python3
"""Measure same-sequence predictor effects on audited native structural states."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha
from prepare_paired_phylogenetic_inputs import write_table


def compare_states(a, b, plddt, pae):
    if a['sequence'] != b['sequence']:
        raise ValueError('Predictor sequences differ')
    n = len(a['sequence'])
    for d in (a,b):
        if any(np.asarray(d[k]).shape != (n,) for k in ['states_array','valid','feature_min_plddt','feature_max_pae','partner_residue_1based']):
            raise ValueError('Feature array dimensions differ')
        valid = d['valid']
        if valid.dtype != bool or not np.isfinite(d['feature_min_plddt'][valid]).all() or not np.isfinite(d['feature_max_pae'][valid]).all():
            raise ValueError('Invalid feature confidence')
    mask = a['valid'] & b['valid'] & (a['feature_min_plddt']>=plddt) & (b['feature_min_plddt']>=plddt)
    if pae is not None:
        mask &= (a['feature_max_pae']<=pae) & (b['feature_max_pae']<=pae)
    mismatch = a['states_array'] != b['states_array']
    partner_changed = a['partner_residue_1based'] != b['partner_residue_1based']
    retained = int(mask.sum())
    row = {'plddt_cutoff':plddt,'pae_cutoff':pae if pae is not None else 'unfiltered',
           'protein_length':n,'retained_residues':retained,'retained_fraction':retained/n,
           'status':'compared' if retained>=50 and retained>=n*.5 else 'insufficient_coverage',
           'state_mismatches':int((mask&mismatch).sum()),
           'state_mismatch_fraction':float(mismatch[mask].mean()) if retained else '',
           'partner_changes':int((mask&partner_changed).sum()),
           'partner_change_fraction':float(partner_changed[mask].mean()) if retained else ''}
    for label, context in [('same_partner',~partner_changed),('changed_partner',partner_changed)]:
        subset=mask&context;count=int(subset.sum())
        row[label+'_residues']=count
        row[label+'_state_mismatches']=int((subset&mismatch).sum())
        row[label+'_state_mismatch_fraction']=float(mismatch[subset].mean()) if count else ''
    return row, mask


def load_encoding(row):
    p=ROOT/row['encoding_path']
    if sha(p)!=row['encoding_sha256']:raise ValueError('Changed native encoding')
    with np.load(p,allow_pickle=False) as source:d={k:source[k].copy() for k in source.files}
    d['sequence']=str(d['sequence']);d['states_array']=np.array(list(str(d['states'])))
    if hashlib.sha256(d['sequence'].encode()).hexdigest()!=row['sequence_sha256']:raise ValueError('Encoding sequence differs')
    return d


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for n in ['inputs','reference','local','output']:parser.add_argument('--'+n,type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use new immutable output')
    receipts={k:checked_receipt(getattr(args,k)) for k in ['inputs','reference','local']}
    for k in ['reference','local']:
        if receipts[k]['status']!='complete_native_3di_feature_audit':raise ValueError('Complete qualified native encodings required')
    if receipts['reference']['mapping_receipt_sha256']!=receipts['inputs']['source_mapping_receipt_sha256']:
        raise ValueError('Different reference mapping')
    # Bind native provenance to the exact complete control inference configuration.
    native=ROOT/receipts['local']['native_path']
    if sha(native/'config.json')!=receipts['local']['native_config_sha256']:
        raise ValueError('Changed local native configuration')
    provenance=native/'model_provenance.json'
    if sha(provenance)!=receipts['local']['native_artifacts']['model_provenance.json']:
        raise ValueError('Changed local native model provenance')
    local_models=json.loads(provenance.read_text())
    for model in local_models:
        pred=ROOT/model['prediction_receipt_path'];config=pred.parent/'config.json'
        if model['provider']!='local' or model['tool']!='ESMFold v1' or sha(pred)!=model['prediction_receipt_sha256']:
            raise ValueError('Wrong local prediction source')
        if sha(config)!=model['prediction_config_sha256'] or json.loads(config.read_text())['input_receipt_sha256']!=sha(args.inputs/'receipt.json'):
            raise ValueError('Local prediction configuration is not bound to controls')
    indexes={}
    for label in ['reference','local']:
        with (getattr(args,label)/'model_summary.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        indexes[label]={(r['model_id'],str(r['version'])):r for r in rows}
        if len(indexes[label])!=len(rows):raise ValueError('Duplicate model identities')
    if set(indexes['local'])!={(m['model_id'],str(m['version'])) for m in local_models}:raise ValueError('Local provenance model universe differs')
    local_by_seq={r['sequence_sha256']:r for r in indexes['local'].values()}
    if len(local_by_seq)!=len(indexes['local']):raise ValueError('Multiple local models per sequence')
    models=json.loads((args.inputs/'reference_models.json').read_text())
    seq={r.id:str(r.seq) for r in SeqIO.parse(args.inputs/'candidates.faa','fasta')}
    expected={m['sequence_sha256'] for m in models}
    if len(expected)!=len(models) or expected!=set(local_by_seq) or {'S'+s for s in expected}!=set(seq):raise ValueError('Different control sequence universes')
    rows=[];confusion=Counter()
    for m in models:
        af=indexes['reference'][(m['model_id'],str(m['version']))];esm=local_by_seq[m['sequence_sha256']]
        a,b=load_encoding(af),load_encoding(esm)
        if a['sequence']!=seq['S'+m['sequence_sha256']] or b['sequence']!=a['sequence']:raise ValueError('Control sequence identity differs')
        for plddt in (0,70,90):
            for pae in (None,5,10,15):
                r,mask=compare_states(a,b,plddt,pae)
                rows.append({'sequence_sha256':m['sequence_sha256'],'reference_model_id':af['model_id'],'local_model_id':esm['model_id'],**r})
                if plddt==70 and pae==10 and r['status']=='compared':
                    confusion.update(zip(a['states_array'][mask],b['states_array'][mask]))
    summary=[]
    for plddt in (0,70,90):
        for pae in ('unfiltered',5,10,15):
            subset=[r for r in rows if r['plddt_cutoff']==plddt and r['pae_cutoff']==pae]
            eligible=[r for r in subset if r['status']=='compared']
            summary.append({'plddt_cutoff':plddt,'pae_cutoff':pae,'proteins':len(subset),'compared':len(eligible),
              'coverage_excluded':len(subset)-len(eligible),
              'median_state_mismatch_fraction':float(np.median([r['state_mismatch_fraction'] for r in eligible])) if eligible else '',
              'median_partner_change_fraction':float(np.median([r['partner_change_fraction'] for r in eligible])) if eligible else ''})
    args.output.mkdir(parents=True)
    write_table(args.output/'protein_comparisons.tsv',rows)
    write_table(args.output/'threshold_summary.tsv',summary)
    # Fixed 20-state alphabet grid includes zero cells; these letters are not amino acids.
    alphabet='ACDEFGHIKLMNPQRSTVWY'
    write_table(args.output/'state_confusion_plddt70_pae10.tsv',[{'reference_state':a,'local_state':b,'residues':confusion[a,b]} for a in alphabet for b in alphabet])
    result={'status':'complete_same_sequence_predictor_alphabet_comparison','proteins':len(models),'comparison_rows':len(rows),
        'source_receipts':{k:{'path':str(getattr(args,k)),'sha256':sha(getattr(args,k)/'receipt.json')} for k in receipts},
        'script_sha256':sha(Path(__file__)),
        'interpretation':'Same-sequence predictor disagreement, not evolutionary substitutions or experimental accuracy. Joint masks require valid native features and six-residue confidence in both models. PAE already represents maximum directional context error. Compared summaries require >=50 and >=half the full protein; cohort changes across thresholds. Partner changes can accompany state changes but are not established causes. No independence, significance or population calibration claim.',
        'artifacts':{p.name:sha(p) for p in args.output.iterdir() if p.is_file()}}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
