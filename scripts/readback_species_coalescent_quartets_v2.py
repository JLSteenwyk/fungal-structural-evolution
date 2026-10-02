#!/usr/bin/env python3
"""Full independent local and global quartet scoring with numerical readback.

Reuses the immutable local readback and adds an independent exact matching
quartet numerator before accepting each candidate. Native constrained-search
optimality, gene/model uncertainty and biological rooting remain unproven.
"""
import argparse
import json
from pathlib import Path
import time

from coalescent_global_quartet_audit import matching_quartets,species_tripartitions
from coalescent_quartet_audit_v2 import packed,read_tree,resolved_quartet_denominator
from readback_species_coalescent_quartets import audit_case,source_bundle,write,scalar_log
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args()
    plan=json.loads(args.plan.read_text())
    native_plan,cases,bindings=source_bundle(plan,args.plan)
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    write(out/'config.json',dict(plan_sha256=sha(args.plan),source_hashes=bindings))
    bind(bindings,out/'config.json')
    summaries=[]
    for case in cases:
        started=time.monotonic()
        native_root=Path(native_plan['output'])/case['case']
        rp=native_root/'receipt.json'
        producer=json.loads(rp.read_text())
        assert producer['status']=='complete_native_coalescent_species_case_pending_full_independent_scoring'
        assert producer['case']==case and producer['returncode']==0
        for name,digest in producer['artifacts'].items():
            bind(bindings,native_root/name,digest)
        bind(bindings,rp)
        verify({str(native_root/name):digest for name,digest in producer['artifacts'].items()})
        roles=json.loads(Path(case['taxa_roles']).read_text())
        tips=sorted(roles['taxa'])
        genes=Path(case['genes']).read_text().splitlines()
        order=json.loads(Path(case['gene_order']).read_text())
        assert len(genes)==len(order)==125
        species=read_tree(path=native_root/'species.tree')
        partitions=species_tripartitions(species,tips)
        global_rows=[]
        for marker,newick in zip(order,genes):
            gene=read_tree(data=newick)
            score=int(matching_quartets(partitions,*packed(gene,tips)))
            denominator=resolved_quartet_denominator(gene)
            assert 0 <= score <= denominator
            global_rows.append(dict(marker=marker,matching_quartets=score,resolved_quartets=denominator))
        log=(native_root/'native.log').read_text()
        native_score=int(scalar_log(log,r'^Final quartet score is: (\d+)$'))
        assert sum(r['matching_quartets'] for r in global_rows)==native_score
        native_denominator=int(scalar_log(log,r'^Number of quartet trees in the gene trees: (\d+)$'))
        assert sum(r['resolved_quartets'] for r in global_rows)==native_denominator
        local=audit_case(native_plan,plan['native_plan'],case,out,bindings)
        folder=out/case['case']
        write(folder/'global_quartets.json',dict(independent_matching_score=native_score,
            independent_resolved_denominator=native_denominator,markers=global_rows,
            scope='Exact global matching numerator and resolved denominator from independent node-component/LCA counts. No claim of globally unconstrained search optimality or biological model adequacy.'))
        bind(bindings,folder/'global_quartets.json')
        summary=dict(case=case['case'],taxa=len(tips),roles=case['roles'],markers=125,
            internal_branches=local['internal_branches'],branch_gene_states=local['branch_gene_states'],
            zero_evidence_branches=local['zero_evidence_branches'],
            independent_global_matching_score=native_score,independent_global_resolved_denominator=native_denominator,
            native_normalized_score=local['native_normalized_score'],
            global_numerator_independently_verified=True,native_search_optimality_verified=False,
            local_readback_receipt=str(folder/'receipt.json'),local_readback_receipt_sha256=sha(folder/'receipt.json'),
            global_readback=str(folder/'global_quartets.json'),global_readback_sha256=sha(folder/'global_quartets.json'),
            elapsed_seconds=time.monotonic()-started)
        summaries.append(summary)
        print('complete_global_and_local_quartet_case_verified',case['case'],len(summaries),'/',len(cases),flush=True)
    verify(bindings)
    result=dict(status=('passed_full_30_coalescent_global_local_quartet_numeric_readback'
        if plan['mode']=='full_batch' else 'passed_complete_named_native_case_global_local_quartet_preflight'),
        plan_sha256=sha(args.plan),native_plan_sha256=sha(plan['native_plan']),cases=len(cases),
        markers_per_case=125,internal_branches=sum(r['internal_branches'] for r in summaries),
        branch_gene_states=sum(r['branch_gene_states'] for r in summaries),summaries=summaries,
        source_hashes=bindings,global_numerator_independently_verified=True,
        native_search_optimality_verified=False,scientific_eligibility=False,scope=plan['scope'])
    if plan['mode']=='full_batch':
        assert result['cases']==30 and result['internal_branches']==15510
        assert result['branch_gene_states']==1938750
    write(out/'receipt.json',result)
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,float,bool))}),flush=True)


if __name__=='__main__':main()
