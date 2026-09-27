#!/usr/bin/env python3
"""Compare exact gene-pair/Pfam identities across guides and describe sequence/structure directions."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha

KEY=['family','gene_a','gene_b','pfam_accession','screen']
COPY=['gene_node','taxon_id','chosen_reference_gene','reference_choices','complete','eligible_contrast_min','eligible_contrast_max','duplicate_pair_sequence_distance','sequence_signed_difference','sequence_normalized_contrast','sequence_direction','sequence_covariate_status']
STRUCT={'a_farther_from_reference':1,'b_farther_from_reference':-1}
SEQ={'a_longer':1,'b_longer':-1}


def relation(row,margin):
    if row is None:return 'missing_guide'
    structural=STRUCT.get(row[margin]);sequence=SEQ.get(row['sequence_direction'])
    if structural is None:return 'structural_direction_unresolved'
    if sequence is None or row['sequence_covariate_status']!='resolved_pair_distance':return 'sequence_direction_unresolved'
    return 'concordant' if structural==sequence else 'discordant'


def compare(a,b,margin):
    presence='both_guides' if a and b else ('mafft_only' if a else 'profile_only')
    if a and b:
        sa,sb=STRUCT.get(a[margin]),STRUCT.get(b[margin])
        structural=('same_stable_direction' if sa==sb else 'opposite_stable_directions') if sa is not None and sb is not None else 'unresolved_in_one_or_both'
        qa,qb=SEQ.get(a['sequence_direction']),SEQ.get(b['sequence_direction'])
        sequence=('same_resolved_direction' if qa==qb else 'opposite_resolved_directions') if qa is not None and qb is not None and a['sequence_covariate_status']==b['sequence_covariate_status']=='resolved_pair_distance' else 'unresolved_in_one_or_both'
    else:structural=sequence='missing_guide'
    ra,rb=relation(a,margin),relation(b,margin)
    joint=ra+'_in_both_guides' if ra==rb and ra in ['concordant','discordant'] else ('guide_sensitive_relationship' if {ra,rb}=={'concordant','discordant'} else 'unresolved_or_missing')
    return dict(presence=presence,structural_guide_agreement=structural,sequence_guide_agreement=sequence,mafft_sequence_structure_relation=ra,profile_sequence_structure_relation=rb,joint_sequence_structure_relation=joint)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();p=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for path,h in p['pins'].items():
            if sha(path)!=h:raise ValueError('Changed input: '+path)
    verify();root=Path(p['source']);r=json.loads((root/'receipt.json').read_text());a=json.loads(Path(p['readback']).read_text())
    if a['status']!='passed_full_domain_event_robustness_readback' or a['producer_receipt_sha256']!=sha(root/'receipt.json'):raise ValueError('Unbound event readback')
    path=root/'event_domain_robustness.tsv'
    if sha(path)!=r['artifacts'][path.name]:raise ValueError('Changed event-domain table')
    groups={};source_rows=0
    for row in csv.DictReader(path.open(),delimiter='\t'):
        key=tuple(row[k] for k in KEY);guide=row['guide']
        if guide not in ['mafft','profile'] or row['gene_a']>=row['gene_b']:raise ValueError('Unexpected guide or gene orientation')
        if guide in groups.setdefault(key,{}):raise ValueError('Ambiguous guide-specific event for exact gene pair')
        groups[key][guide]=row;source_rows+=1
    out=Path(p['output']);out.mkdir(parents=True,exist_ok=False);counts=Counter();rows=0;presence_counts=Counter()
    extra=list(compare(None,{'sequence_direction':'unresolved','sequence_covariate_status':'unresolved','m':'incomplete'},'m'))
    fields=KEY+['margin']+[g+'_'+f for g in ['mafft','profile'] for f in COPY+['structural_direction']]+extra
    with (out/'guide_comparison.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for key,guides in sorted(groups.items()):
            ma,pr=guides.get('mafft'),guides.get('profile')
            if ma and pr and ma['taxon_id']!=pr['taxon_id']:raise ValueError('Taxon identity differs')
            for margin in p['margins']:
                x=dict(zip(KEY,key));x['margin']=margin;x.update(compare(ma,pr,margin))
                for guide,row in [('mafft',ma),('profile',pr)]:
                    for field in COPY:x[guide+'_'+field]=row[field] if row else ''
                    x[guide+'_structural_direction']=row[margin] if row else ''
                w.writerow(x);rows+=1;presence_counts[x['screen']+'|'+margin+'|'+x['presence']]+=1
                for field in extra[1:]:counts['|'.join([x['screen'],margin,field,x[field]])]+=1
    if source_rows!=r['event_domain_screen_rows'] or rows!=len(groups)*len(p['margins']):raise ValueError('Incomplete guide comparison')
    verify()
    result=dict(status='complete_domain_guide_direction_comparison_pending_readback',plan_sha256=ph,source_rows=source_rows,gene_pair_domain_screen_union=len(groups),comparison_rows=rows,presence_counts=dict(presence_counts),classification_counts=dict(counts),artifacts={'guide_comparison.tsv':sha(out/'guide_comparison.tsv')},scope='Exact family/lexical-gene-pair/Pfam/screen outer join; guide-specific nodes retained, never used as cross-guide identifiers. All margins retained. Whole-protein sequence direction compared descriptively with domain reference-similarity direction; unresolved cases explicit. Shared guides are sensitivity alternatives, not independent replication; no association significance, rate or causal inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['presence_counts','classification_counts']},indent=2))


if __name__=='__main__':main()
