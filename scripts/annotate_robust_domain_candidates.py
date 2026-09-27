#!/usr/bin/env python3
"""Annotate every cross-guide robust domain candidate with Pfam and fit ranges."""
import argparse,csv,json,hashlib
from pathlib import Path
from collections import defaultdict


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open() as f:yield from csv.DictReader(f,delimiter='\t')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text())
    for name,digest in p['pins'].items():assert sha(name)==digest,name
    out=Path(p['output']);assert not out.exists()
    pfams={};record={}
    for line in Path(p['pfam']).read_text().splitlines():
        if line=='//':
            assert record['AC'] not in pfams;pfams[record['AC']]=record;record={}
        elif line.startswith('#=GF '):
            _,key,value=line.split(maxsplit=2)
            if key in ['AC','ID','DE','TP','CL']:record[key]=value.strip()
    selected=[r for r in rows(p['candidates']) if r['screen']=='n30_c70' and r['margin']=='direction_margin_0_1' and r['structural_guide_agreement']=='same_stable_direction']
    assert len(selected)==947
    links=defaultdict(set)
    for r in rows(p['links']):links[tuple(r[k] for k in ['guide','family','gene_node','gene_a','gene_b','pfam_accession'])].add(r['triad_key'])
    triads={};required=set()
    for r in selected:
        key=tuple(r[k] for k in ['family','gene_a','gene_b','pfam_accession']);assert key not in triads
        sets=[links[(guide,r['family'],r[guide+'_gene_node'],r['gene_a'],r['gene_b'],r['pfam_accession'])] for guide in ['mafft','profile']]
        assert all(sets);triads[key]=set.union(*sets);required.update(triads[key])
    metrics=['common_residues','coverage_a','coverage_b','coverage_reference','rmsd_ab','rmsd_ar_minus_br','sequence_identity_ab','mean_plddt_a','mean_plddt_b','mean_plddt_reference','joint_plddt70_fraction']
    fits=defaultdict(list);seen=set()
    for r in rows(p['fits']):
        if r['triad_key'] not in required:continue
        key=tuple(r[k] for k in ['triad_key','mask','order_ab','order_ar','order_br','mapping_definition']);assert key not in seen;seen.add(key)
        assert r['n30_c70_pass']=='1' and r['fit_status']=='computed_unique_at_numeric_tolerance'
        fits[r['triad_key']].append({k:float(r[k]) for k in metrics})
    assert set(fits)==required and all(len(v)==32 for v in fits.values())
    output=[];summaries=defaultdict(list)
    for r in selected:
        key=tuple(r[k] for k in ['family','gene_a','gene_b','pfam_accession']);meta=pfams[r['pfam_accession']];pool=[x for triad in triads[key] for x in fits[triad]]
        x=dict(r,pfam_name=meta['ID'],pfam_description=meta['DE'],pfam_type=meta.get('TP',''),pfam_clan=meta.get('CL',''),unique_interval_triads=len(triads[key]),unique_fit_alternatives=len(pool),triad_keys_json=json.dumps(sorted(triads[key]),separators=(',',':')))
        for metric in metrics:
            x['all_alternatives_min_'+metric]=min(t[metric] for t in pool);x['all_alternatives_max_'+metric]=max(t[metric] for t in pool)
        low=x['all_alternatives_min_rmsd_ar_minus_br'];high=x['all_alternatives_max_rmsd_ar_minus_br']
        assert abs(low-min(float(r[g+'_eligible_contrast_min']) for g in ['mafft','profile']))<1e-12
        assert abs(high-max(float(r[g+'_eligible_contrast_max']) for g in ['mafft','profile']))<1e-12
        assert low>.1 or high<-.1
        x['minimum_absolute_contrast_angstrom']=min(abs(low),abs(high))
        output.append(x);summaries[r['study_role'],r['candidate_class'],r['pfam_accession']].append(x)
    groups=[]
    for (role,category,pfam),xs in sorted(summaries.items()):
        m=pfams[pfam];groups.append(dict(study_role=role,candidate_class=category,pfam_accession=pfam,pfam_name=m['ID'],pfam_description=m['DE'],event_domain_combinations=len(xs),distinct_gene_pairs=len({(x['family'],x['gene_a'],x['gene_b']) for x in xs}),families=len({x['family'] for x in xs}),taxa=len({x['taxon'] for x in xs})))
    out.mkdir(parents=True)
    for name,rs in [('candidate_annotations.tsv',sorted(output,key=lambda x:(x['family'],x['gene_a'],x['gene_b'],x['pfam_accession']))),('pfam_summary.tsv',groups)]:
        with (out/name).open('w') as f:
            w=csv.DictWriter(f,list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
    receipt=dict(status='complete_robust_domain_candidate_annotations_pending_readback',plan_sha256=sha(a.plan),candidates=len(output),unique_interval_triads=len(required),unique_fit_alternatives=len(seen),pfam_summary_rows=len(groups),pfam_release='38.2',artifacts={f.name:sha(f) for f in out.iterdir()},scope='All 947 stable cross-guide event/domain candidates at n30/c70 and descriptive 0.1 A margin; no subset chosen by function or quality. Minima/maxima cover every tied-reference/annotation interval triad and all 32 fitting alternatives, deduplicated across guides. Pfam labels describe families, not experimentally established functions. Confidence ranges are descriptive, not uncertainty calibration or evidence of selection.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
