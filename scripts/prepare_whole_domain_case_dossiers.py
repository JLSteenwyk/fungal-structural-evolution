#!/usr/bin/env python3
"""Prepare every opposing-scale case at the declared descriptive 0.1-A margin."""
import csv,json
from collections import defaultdict
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha

KEY=['family','gene_a','gene_b','pfam_accession']
METRICS=['common_residues','coverage_a','coverage_b','coverage_reference','rmsd_ab','rmsd_ar_minus_br','sequence_identity_ab','mean_plddt_a','mean_plddt_b','mean_plddt_reference','joint_plddt70_fraction']


def rows(path):
    with Path(path).open() as f:yield from csv.DictReader(f,delimiter='\t')


def main():
    roots={k:Path('results/structural_comparisons')/v for k,v in dict(integration='whole-domain-contrast-integration-20260927-v1',whole_maps='whole-protein-common-residues-20260927-v1',domain_maps='duplication-domain-common-residues-20260927-v1',whole_fits='whole-protein-common-core-fits-20260927-v1',domain_fits='duplication-domain-common-core-fits-20260927-v1').items()}
    sources={}
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name
        assert sha(p)==r['artifacts'][name];sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return p
    all_rows=list(rows(checked(roots['integration'],'domain_comparisons.tsv')))
    selected={tuple(r[k] for k in KEY) for r in all_rows if r['margin']=='direction_margin_0_1' and r['whole_domain_relationship']=='opposite_direction_between_scales'}
    assert len(selected)==13
    sensitivity=[r for r in all_rows if tuple(r[k] for k in KEY) in selected];assert len(sensitivity)==13*18
    cases={key:[r for r in sensitivity if tuple(r[k] for k in KEY)==key] for key in selected}
    annotation_plan=Path('metadata/duplication_domain_candidate_annotation_plan_20260927.json');plan=json.loads(annotation_plan.read_text())
    pfam=Path(plan['pfam']);assert sha(pfam)==plan['pins'][str(pfam)];sources[str(pfam)]=sha(pfam)
    annotations={};record={}
    for line in pfam.read_text().splitlines():
        if line=='//':annotations[record['AC']]=record;record={}
        elif line.startswith('#=GF '):
            _,key,value=line.split(maxsplit=2)
            if key in ['AC','ID','DE','TP','CL']:record[key]=value.strip()
    links={};required={};case_triads={}
    for scale,filename,tk in [('whole','event_reference_triads.tsv','triad_id'),('domain','event_domain_links.tsv','triad_key')]:
        source=list(rows(checked(roots[scale+'_maps'],filename))); index=defaultdict(set)
        for r in source:index[tuple(r[k] for k in (KEY[:3] if scale=='whole' else KEY))].add(r[tk])
        case_triads[scale]={k:index[k[:3] if scale=='whole' else k] for k in selected}
        assert all(case_triads[scale].values());required[scale]=set.union(*case_triads[scale].values())
        links[scale]=[r for r in source if r[tk] in required[scale]]
    fits={};exported={}
    for scale,tk in [('whole','triad_id'),('domain','triad_key')]:
        grouped=defaultdict(list)
        for r in rows(checked(roots[scale+'_fits'],'common_residue_fits.tsv')):
            if r[tk] in required[scale]:
                assert r['n30_c50_pass']=='1' and r['fit_status']=='computed_unique_at_numeric_tolerance'
                grouped[r[tk]].append(r)
        assert set(grouped)==required[scale] and all(len(v)==32 for v in grouped.values())
        fits[scale]=grouped;exported[scale]=[r for key in sorted(grouped) for r in grouped[key]]
    dossiers=[]
    for key in sorted(selected):
        observations=cases[key];base=observations[0];meta=annotations[key[-1]]
        r=dict(zip(KEY,key),species_name=base['species_name'],study_role=base['study_role'],lineage_group=base['lineage_group'],pfam_name=meta['ID'],pfam_description=meta['DE'],pfam_type=meta.get('TP',''),pfam_clan=meta.get('CL',''))
        r['opposing_screens_at_margin_0_1']=';'.join(sorted(x['screen'] for x in observations if x['margin']=='direction_margin_0_1' and x['whole_domain_relationship']=='opposite_direction_between_scales'))
        for scale in ['whole','domain']:
            pool=[row for tk in sorted(case_triads[scale][key]) for row in fits[scale][tk]]
            r[scale+'_unique_triads']=len(case_triads[scale][key]);r[scale+'_fit_alternatives']=len(pool)
            for metric in METRICS:
                values=[float(row[metric]) for row in pool];r[scale+'_min_'+metric]=min(values);r[scale+'_max_'+metric]=max(values)
        dossiers.append(r)
    out=Path('results/structural_comparisons/whole-domain-case-dossiers-20260927-v1');out.mkdir(exist_ok=False);artifacts={}
    tables={'case_dossiers.tsv':dossiers,'all_sensitivity_settings.tsv':sensitivity}
    for scale in ['whole','domain']:tables[scale+'_reference_links.tsv']=links[scale];tables[scale+'_fit_alternatives.tsv']=exported[scale]
    for name,data in tables.items():
        p=out/name
        with p.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
        assert list(rows(p))==[{k:str(v) for k,v in r.items()} for r in data];artifacts[name]=sha(p)
    for path,digest in sources.items():assert sha(path)==digest
    result=dict(status='complete_descriptive_whole_domain_case_dossiers',script_sha256=sha(__file__),source_hashes=sources,cases=len(dossiers),sensitivity_rows=len(sensitivity),whole_unique_triads=len(required['whole']),domain_unique_triads=len(required['domain']),whole_fit_rows=len(exported['whole']),domain_fit_rows=len(exported['domain']),artifacts=artifacts,scope='All opposing-scale cases at 0.1 A in any of six coverage screens; all 18 settings retained. Pfam labels are annotation hypotheses, not experimentally established function. Extrema retain all linked references and fit alternatives; no uncertainty calibration, statistical selection, domain-orientation mechanism or significance claim.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
