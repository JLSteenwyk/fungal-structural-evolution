#!/usr/bin/env python3
"""Compare whole-protein and domain contrasts at matched margins and exact reference sets."""
import csv
import json
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha

KEY=['family','gene_a','gene_b','screen']
MARGINS={'direction_margin_0':0.,'direction_margin_0_01':0.01,'direction_margin_0_1':0.1}


def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def sign(lo,hi,margin):
    if float(lo)>margin:return 1
    if float(hi)<-margin:return -1
    return 0


def main():
    whole=Path('results/structural_comparisons/whole-protein-cross-guide-sensitivity-20260927-v1')
    domain=Path('results/structural_comparisons/duplication-domain-candidate-sampling-20260927-v1')
    maps=Path('results/structural_comparisons/whole-protein-common-residues-20260927-v1')
    integration=Path('results/structural_comparisons/duplication-domain-triad-integration-20260927-v1')
    sources={}
    for root in [whole,domain,maps,integration]:
        rp=root/'receipt.json';sources[str(rp)]=sha(rp)
    dr=json.loads((domain/'receipt.json').read_text());da=Path('metadata/duplication_domain_candidate_sampling_completed_20260927.json')
    assert json.loads(da.read_text())['producer_receipt_sha256']==sha(domain/'receipt.json');sources[str(da)]=sha(da)
    paths=[whole/'all_pair_comparisons.tsv',domain/'annotated_comparisons.tsv',maps/'event_reference_triads.tsv',maps/'model_triads.tsv',integration/'domain_triads.sqlite']
    for path in paths:
        r=json.loads((path.parent/'receipt.json').read_text());assert sha(path)==r['artifacts'][path.name];sources[str(path)]=sha(path)
    wr={tuple(r[k] for k in KEY):r for r in rows(paths[0])};assert len(wr)==121068
    triads={r['triad_id']:r for r in rows(paths[3])};wrefs=defaultdict(set);drefs=defaultdict(set)
    for row in rows(paths[2]):
        key=tuple(row[k] for k in ['guide','family','gene_a','gene_b']);t=triads[row['triad_id']]
        wrefs[key].add((row['reference_gene'],t['a_model'],t['a_version'],t['b_model'],t['b_version'],t['reference_model'],t['reference_version']))
    with sqlite3.connect('file:'+str(paths[4])+'?mode=ro',uri=True) as conn:
        for row in conn.execute('SELECT DISTINCT guide,family,gene_a,gene_b,reference_gene,model_a,version_a,model_b,version_b,reference_model,reference_version FROM triads'):
            drefs[row[:4]].add(row[4:])
    joined=[];counts=Counter();covered=set()
    for row in rows(paths[1]):
        key=tuple(row[k] for k in KEY);w=wr.get(key);margin=MARGINS[row['margin']];covered.add(key)
        out=dict(row);out['whole_pair_present']=int(w is not None)
        same_refs=[];whole_sign=[];domain_sign=[]
        for guide in ['mafft','profile']:
            rk=(guide,)+key[:3];same=bool(wrefs[rk]) and wrefs[rk]==drefs[rk];same_refs.append(same)
            out[guide+'_exact_reference_models_match']=int(same)
            complete=bool(w and w['all_references_eligible_'+guide] and float(w['all_references_eligible_'+guide])==1)
            ws=sign(w['structural_min_'+guide],w['structural_max_'+guide],margin) if complete else 0
            ds={'a_farther_from_reference':1,'b_farther_from_reference':-1}.get(row[guide+'_structural_direction'],0)
            whole_sign.append(ws);domain_sign.append(ds)
            out[guide+'_whole_direction_at_matched_margin']=ws
        if not w:status='whole_pair_missing'
        elif not all(same_refs):status='reference_or_model_set_mismatch'
        elif not (whole_sign[0] and whole_sign[0]==whole_sign[1]):status='whole_direction_incomplete_or_unresolved'
        elif not (domain_sign[0] and domain_sign[0]==domain_sign[1]):status='domain_direction_incomplete_or_unresolved'
        else:status='same_direction_both_scales' if whole_sign[0]==domain_sign[0] else 'opposite_direction_between_scales'
        out['whole_domain_relationship']=status;joined.append(out);counts[row['screen'],row['margin'],status]+=1
    assert len(joined)==77760
    unmatched=[dict(zip(KEY,key)) for key in sorted(set(wr)-covered)]
    summary=[dict(screen=k[0],margin=k[1],relationship=k[2],domain_comparisons=v) for k,v in sorted(counts.items())]
    assert sum(counts.values())==77760
    # Independent recomputation of the direction rule from the exported source ranges.
    for row in joined:
        w=wr.get(tuple(row[k] for k in KEY));margin=MARGINS[row['margin']]
        for guide in ['mafft','profile']:
            expected=0
            if w and w['all_references_eligible_'+guide] and float(w['all_references_eligible_'+guide])==1:
                values=[float(w['structural_min_'+guide]),float(w['structural_max_'+guide])]
                if all(v>margin for v in values):expected=1
                elif all(v < -margin for v in values):expected=-1
            assert row[guide+'_whole_direction_at_matched_margin']==expected
    out=Path('results/structural_comparisons/whole-domain-contrast-integration-20260927-v1');out.mkdir(exist_ok=False);artifacts={}
    for name,data in [('domain_comparisons.tsv',joined),('whole_pair_screens_without_domain_rows.tsv',unmatched),('summary.tsv',summary)]:
        path=out/name
        with path.open('w') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');writer.writeheader();writer.writerows(data)
        assert rows(path)==[{k:str(v) for k,v in row.items()} for row in data];artifacts[name]=sha(path)
    for path,digest in sources.items():assert sha(path)==digest
    result=dict(status='complete_descriptive_whole_domain_contrast_integration',source_hashes=sources,script_sha256=sha(__file__),domain_comparison_rows=len(joined),whole_pair_screen_rows_without_domain=len(unmatched),summaries=summary,artifacts=artifacts,
        scope='All six screens and three matched descriptive margins. Exact reference-gene and A/B/reference model/version sets checked per guide. Whole coverage uses original protein lengths; domain coverage uses domain lengths. Same/opposite contrast directions are descriptive, not proof of orientation change, coupling or biological significance; no subtraction of RMSDs across scales.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
