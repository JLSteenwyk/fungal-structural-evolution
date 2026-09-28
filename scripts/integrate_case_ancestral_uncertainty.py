#!/usr/bin/env python3
"""Attach audited ancestral diagnostics to every existing structural case."""
import csv,hashlib,json
from pathlib import Path


def main():
    pins={}
    def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    def table(root,name):
        root=Path(root);rp=root/'receipt.json';receipt=json.loads(rp.read_text());p=root/name
        assert sha(p)==receipt['artifacts'][name];pins[str(rp)]=sha(rp);pins[str(p)]=sha(p)
        with p.open() as h:return list(csv.DictReader(h,delimiter='\t'))
    cases=table('results/structural_comparisons/whole-domain-case-dossiers-20260927-v1','case_dossiers.tsv')
    optimization=table('results/ancestral/whole-optimization-probability-comparison-full-20260927-v1','comparison_summary.tsv')
    coverage=table('results/ancestral/ancestral-context-coverage-20260927-v1','coverage_summary.tsv')
    context=table('results/ancestral/ancestral-context-sensitivity-figure-20260927-v1','family_summary.tsv')
    cp=Path('results/ancestral/ancestral-context-conflict-localization-20260927-v1');cr=json.loads((cp/'receipt.json').read_text());p=cp/'conflict_dossiers.json';assert sha(p)==cr['artifacts'][p.name];pins[str(cp/'receipt.json')]=sha(cp/'receipt.json');conflicts=json.loads(p.read_text())
    families={r['family'] for r in cases};assert len(cases)==len(families)==13 and families=={r['family'] for r in context}
    output=[]
    for case in cases:
        family=case['family'];r=dict(case)
        for label in ['refined_versus_baseline','alternate_versus_baseline']:
            selected=[x for x in optimization if x['base_job_id'].split('-')[0]==family and x['comparison']==label]
            assert len(selected)==(36 if label=='refined_versus_baseline' else 108)
            prefix='ancestral_'+label
            r.update({prefix+'_node_sites':sum(int(x['columns']) for x in selected),prefix+'_map_changes':sum(int(x['map_disagreements']) for x in selected),prefix+'_maximum_total_variation':max(float(x['maximum_total_variation']) for x in selected)})
        c=next(x for x in context if x['family']==family)
        for k,v in c.items():
            if k!='family':r['ancestral_context_'+k]=v
        cov=[x for x in coverage if x['family']==family]
        r['ancestral_context_no_domain_descendant_comparisons']=sum(int(x['comparisons']) for x in cov if x['coverage_class']!='both_have_known_descendants')
        conf=[x for x in conflicts if x['family']==family]
        r['ancestral_opposing_source_node_coordinate_groups']=len(conf)
        r['ancestral_opposing_exact_coordinate_signatures']=len({x['coordinate_signature_sha256'] for x in conf})
        r['ancestral_readiness']='conditional_amino_acid_diagnostics_only_joint_indel_and_convergence_unresolved'
        output.append(r)
    assert sum(int(r['ancestral_context_comparisons']) for r in output)==269154
    assert sum(int(r['ancestral_context_opposing_with_known_descendants'])+int(r['ancestral_context_opposing_without_domain_descendants']) for r in output)==104
    out=Path('results/structural_comparisons/case-ancestral-uncertainty-integration-20260927-v1');out.mkdir(exist_ok=False)
    with (out/'case_dossiers_with_ancestral_uncertainty.tsv').open('w') as h:
        writer=csv.DictWriter(h,list(output[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(output)
    with (out/'case_dossiers_with_ancestral_uncertainty.tsv').open() as h:back=list(csv.DictReader(h,delimiter='\t'))
    for original,joined in zip(cases,back):assert all(joined[k]==v for k,v in original.items())
    result=dict(status='complete_all_case_ancestral_uncertainty_integration',cases=13,original_case_fields_preserved=len(cases[0]),added_fields=len(output[0])-len(cases[0]),pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All exploratory structural cases and original values retained. Counts are dependent diagnostics, not biological rankings, independent changes, calibrated effect tests or qualified ancestral ensembles. Zero high-confidence conflict does not establish robustness to unresolved indels, roots or model adequacy.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');result.update(completed_receipt_path=str(out/'receipt.json'),completed_receipt_sha256=sha(out/'receipt.json'))
    Path('metadata/case_ancestral_uncertainty_integration_completed_20260927.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['pins','artifacts']}))


if __name__=='__main__':main()
