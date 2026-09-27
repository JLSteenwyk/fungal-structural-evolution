#!/usr/bin/env python3
"""Export all domain-policy annotations for exact ancestral candidate sequences."""
import csv,json,sqlite3
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from prepare_case_ancestral_neighborhoods import sha,read


def main():
    root=Path('results/ancestral/case-sequences-20260927-v1');rp=root/'receipt.json'
    receipt=json.loads(rp.read_text());inventory=root/'sequence_inventory.tsv'
    assert sha(inventory)==receipt['artifacts'][inventory.name]
    cases_path=Path('results/structural_comparisons/whole-domain-case-dossiers-20260927-v1/case_dossiers.tsv')
    cr=cases_path.parent/'receipt.json';assert sha(cases_path)==json.loads(cr.read_text())['artifacts'][cases_path.name]
    cases={r['family']:r for r in read(cases_path)}
    domain=Path('results/domains/full-candidate-architectures-v1');dr=domain/'receipt.json'
    audit=Path('results/domains/full-candidate-architectures-readback-v1/receipt.json')
    proof=json.loads(audit.read_text());assert proof['status']=='passed_complete_independent_candidate_architecture_readback'
    assert proof['production_receipt_sha256']==sha(dr)
    database=domain/'candidate_architectures.sqlite';digest=sha(database)
    assert digest==json.loads(dr.read_text())['artifacts'][database.name]
    pins={str(p):sha(p) for p in [rp,inventory,cases_path,cr,dr,audit]};pins[str(database)]=digest
    connection=sqlite3.connect(database.resolve().as_uri()+'?mode=ro',uri=True);connection.row_factory=sqlite3.Row
    out=Path('results/ancestral/case-domain-context-20260927-v1');out.mkdir(exist_ok=False)
    rows=[];exports=[];policies=['alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore']
    for item in read(inventory):
        gene=item['gene'];taxon,protein=gene.split('_',1)
        records=connection.execute('SELECT * FROM proteins JOIN queries USING(sequence_id) WHERE taxon_id=? AND protein_id=?',(taxon,protein)).fetchall()
        assert len(records)==1;record=dict(records[0]);assert record['sequence_id']=='S'+item['sequence_sha256']
        payload=json.loads(record.pop('candidate_architectures_json'));exports.append(dict(gene=gene,family=item['family'],**record,candidate_architectures=payload))
        assert set(payload['policies'])==set(policies)
        for policy in policies:
            state=payload['policies'][policy];alternative=payload['alternatives'][state['alternative_index']]
            annotations=alternative['annotations'];target=cases[item['family']]['pfam_accession']
            selected=[h for h in annotations if h['pfam_accession']==target]
            for hit in annotations:
                assert 1<=int(hit['envelope_start'])<=int(hit['alignment_start'])<=int(hit['alignment_end'])<=int(hit['envelope_end'])<=int(item['length'])
            rows.append(dict(family=item['family'],gene=gene,length=int(item['length']),policy=policy,target_pfam=target,raw_hits=record['raw_hits'],retained_hits=len(annotations),target_hits=len(selected),target_partial_hmm_hits=sum(Decimal(h['hmm_coverage'])<Decimal('.70') for h in selected),target_max_hmm_coverage=max((h['hmm_coverage'] for h in selected),key=Decimal,default=''),ordered_tokens_json=json.dumps(alternative['ordered_model_tokens'],separators=(',',':')),architecture_signature=alternative['token_signature_sha256'],alignment_overlap_pairs=alternative['alignment_overlap_pairs'],policy_unresolved_overlap_pairs=state['unresolved_overlap_pairs'],policy_primary_rank_ties=state['primary_rank_tie_pairs'],retained_sets_agree=record['retained_sets_agree']))
    connection.close();assert len(exports)==1025 and len(rows)==4100
    (out/'candidate_annotations.json').write_text(json.dumps(exports,indent=2)+'\n')
    grouped=defaultdict(list)
    for row in rows:grouped[row['family'],row['policy']].append(row)
    summaries=[]
    for (family,policy),group in sorted(grouped.items()):
        lookup={r['gene']:r for r in group};case=cases[family];a,b=lookup[case['gene_a']],lookup[case['gene_b']]
        summaries.append(dict(family=family,policy=policy,proteins=len(group),with_target_pfam=sum(r['target_hits']>0 for r in group),without_target_pfam=sum(r['target_hits']==0 for r in group),multiple_target_hits=sum(r['target_hits']>1 for r in group),with_partial_target_hit=sum(r['target_partial_hmm_hits']>0 for r in group),without_any_retained_hit=sum(r['retained_hits']==0 for r in group),unique_ordered_architectures=len({r['architecture_signature'] for r in group}),matching_focal_a=sum(r['architecture_signature']==a['architecture_signature'] for r in group),matching_focal_b=sum(r['architecture_signature']==b['architecture_signature'] for r in group),focal_architectures_equal=int(a['architecture_signature']==b['architecture_signature']),policy_sensitive_proteins=sum(not r['retained_sets_agree'] for r in group)))
    for filename,data in [('protein_policy_context.tsv',rows),('family_policy_summary.tsv',summaries)]:
        with (out/filename).open('w') as f:w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
        assert read(out/filename)==[{k:str(v) for k,v in row.items()} for row in data]
    # Independent aggregation from complete JSON annotations, rather than the
    # derived policy table, checks every family/policy target and token count.
    for summary in summaries:
        data=[r for r in exports if r['family']==summary['family']];policy=summary['policy'];target=cases[summary['family']]['pfam_accession']
        alts=[r['candidate_architectures']['alternatives'][r['candidate_architectures']['policies'][policy]['alternative_index']] for r in data]
        counts=[sum(h['pfam_accession']==target for h in a['annotations']) for a in alts]
        assert sum(n>0 for n in counts)==summary['with_target_pfam'] and sum(n==0 for n in counts)==summary['without_target_pfam']
        assert sum(n>1 for n in counts)==summary['multiple_target_hits']
        assert len({tuple(tuple(t) for t in a['ordered_model_tokens']) for a in alts})==summary['unique_ordered_architectures']
    for p,h in pins.items():assert sha(p)==h
    result=dict(status='complete_ancestral_candidate_domain_context',proteins=1025,protein_policy_rows=4100,family_policy_rows=len(summaries),source_hashes=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Exact sequence-hash join to fully audited Pfam candidate architectures; all four overlap policies and all model types retained. Target absence is no retained annotation, not domain loss. HMM coverage below 0.70 flags partial model alignment, not proven protein truncation. Ordered architecture agreement is not functional equivalence or ASR eligibility.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in summaries:
        if row['policy']=='alignment_evalue':print(row,flush=True)


if __name__=='__main__':main()
