#!/usr/bin/env python3
"""Recompute every local quartet with the pinned native effective-N rule.

All 30 full candidate estimates are required by a full-batch plan. A separately
named preflight may verify complete specified cases; it never certifies the
full batch. Global quartet-score numerator and search optimality remain
unverified: only the exact global resolved denominator and score arithmetic
are checked here. NNI alternatives are checked as an unordered pair.
"""
import argparse
from collections import Counter
import gzip
import json
import math
from pathlib import Path
import re
import time

import numpy as np
import psutil

from coalescent_quartet_audit_v2 import (annotated_branches, colored_quartets, close,
                                       packed, read_tree,
                                       resolved_quartet_denominator)
from coalescent_quartet_numeric_v3 import numerical_check
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def write(path, record):
    with Path(path).open('x') as handle:
        handle.write(json.dumps(record,indent=2,allow_nan=False)+'\n')


def scalar_log(log, pattern):
    values = re.findall(pattern,log,re.MULTILINE)
    assert len(values) == 1, (pattern,len(values))
    return values[0]


def source_bundle(plan, plan_path):
    native_plan = json.loads(Path(plan['native_plan']).read_text())
    completion = json.loads(Path(native_plan['input_completion']).read_text())
    assert completion['status'] == 'complete_verified_full_species_coalescent_inputs'
    assert completion['cases'] == 30 and completion['marker_states'] == 3750
    assert completion['exact_process_journals_checked'] == 2
    bindings = dict(plan['pins'])
    bind(bindings,plan_path)
    bind(bindings,plan['native_plan'])
    bind(bindings,native_plan['input_completion'])
    bind(bindings,completion['full_hash_archive'],completion['full_hash_archive_sha256'])
    archive = json.loads(Path(completion['full_hash_archive']).read_text())
    assert len(archive['services']) == 2
    for path,digest in archive['source_hashes'].items(): bind(bindings,path,digest)
    inputs = json.loads(Path(completion['producer_receipt']).read_text())
    cases = inputs['summaries']
    assert len(cases) == 30 and len({c['case'] for c in cases}) == 30
    root = Path(native_plan['output'])
    if plan['mode'] == 'full_batch':
        producer = json.loads((root/'receipt.json').read_text())
        assert producer['status'] == 'complete_full_native_coalescent_species_sensitivities_pending_independent_scoring'
        assert producer['cases'] == len(producer['runs']) == 30
        assert producer['plan_sha256'] == sha(plan['native_plan'])
        assert producer['input_completion_sha256'] == sha(native_plan['input_completion'])
        assert [r['case'] for r in producer['runs']] == cases
        bind(bindings,root/'receipt.json')
        for path,digest in producer['source_hashes'].items(): bind(bindings,path,digest)
    else:
        assert plan['mode'] == 'complete_named_cases_preflight'
        assert len(plan['cases']) == len(set(plan['cases'])) > 0
        assert set(plan['cases']) <= {c['case'] for c in cases}
        cases = [c for c in cases if c['case'] in plan['cases']]
    verify(bindings)
    return native_plan,cases,bindings


def audit_case(native_plan, native_plan_path, case, out, bindings):
    started = time.monotonic()
    native_root = Path(native_plan['output'])/case['case']
    rp = native_root/'receipt.json'
    receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_native_coalescent_species_case_pending_full_independent_scoring'
    assert receipt['case'] == case and receipt['returncode'] == 0
    assert receipt['taxa'] == case['expected_taxa']
    bind(bindings,rp)
    for name,digest in receipt['artifacts'].items(): bind(bindings,native_root/name,digest)
    assert set(receipt['artifacts']) == {'species.tree','native.log','native_launch.json'}
    launch = json.loads((native_root/'native_launch.json').read_text())
    expected_command = [native_plan['java'],'-XX:ActiveProcessorCount=2','-XX:+ExitOnOutOfMemoryError','-Xmx28G',
        '-jar',native_plan['jar'],'-i',str(Path(case['genes']).resolve()),
        '-o',str((native_root/'species.tree').resolve()),'-t','2','-s',str(native_plan['seed'])]
    assert launch['command'] == launch['cmdline'] == receipt['command'] == expected_command
    assert launch['plan_sha256'] == sha(native_plan_path)
    try:
        original = psutil.Process(launch['pid'])
        assert original.create_time() != launch['created'] or original.status() == psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess: pass
    roles = json.loads(Path(case['taxa_roles']).read_text())
    tips = sorted(roles['taxa'])
    assert len(tips) == len(set(tips)) == case['expected_taxa']
    assert set(roles['outgroups']) <= set(tips)
    assert len(roles['outgroups']) == len(set(roles['outgroups'])) == case['roles']['outgroup']
    assert len(tips)-len(roles['outgroups']) == case['roles']['ingroup']
    assert roles['roles'] == case['roles']
    newicks = Path(case['genes']).read_text().splitlines()
    order = json.loads(Path(case['gene_order']).read_text())
    assert len(newicks) == len(order) == case['markers'] == 125
    genes = [read_tree(data=n) for n in newicks]
    arrays = [packed(g,tips) for g in genes]
    assert set.union(*({n.taxon.label for n in g.leaf_node_iter()} for g in genes)) == set(tips)
    global_denominator = sum(resolved_quartet_denominator(g) for g in genes)
    log = (native_root/'native.log').read_text()
    assert 'ASTRAL version 5.7.8' in log and 'ASTRAL finished in ' in log
    assert int(scalar_log(log,r'^Number of gene trees: (\d+)$')) == 125
    assert int(scalar_log(log,r'^Number of taxa: (\d+) \(\d+ species\)$')) == len(tips)
    assert int(scalar_log(log,r'^Number of quartet trees in the gene trees: (\d+)$')) == global_denominator
    score = int(scalar_log(log,r'^Final quartet score is: (\d+)$'))
    normalized_score = float(scalar_log(log,r'^Final normalized quartet score is: ([^\s]+)$'))
    assert 0 <= score <= global_denominator and global_denominator > 0
    close(normalized_score,score/global_denominator)
    branches = annotated_branches(read_tree(path=native_root/'species.tree'),tips)
    folder = out/case['case']; folder.mkdir()
    branch_rows, states = [], 0
    zero_evidence = 0
    with gzip.open(folder/'gene_branch_quartets.jsonl.gz','xt',compresslevel=6) as handle:
        for branch in branches:
            sums = np.zeros(3)
            gene_counts = Counter()
            for marker,gene,packed_gene in zip(order,genes,arrays):
                counts,available = colored_quartets(*packed_gene,branch['colors'])
                assert np.all(counts >= 0) and counts.sum() <= available
                if available:
                    sums += counts/available
                    disposition = 'resolved_evidence' if counts.sum() else 'unresolved_only'
                else: disposition = 'missing_incident_clade'
                gene_counts[disposition] += 1
                handle.write(json.dumps(dict(split=branch['split'],marker=marker,
                    available_quartets=int(available),resolved_counts=[int(c) for c in counts],
                    unresolved_quartets=int(available-counts.sum()),disposition=disposition),
                    separators=(',',':'),allow_nan=False)+'\n')
                states += 1
            available_genes = gene_counts['resolved_evidence'] + gene_counts['unresolved_only']
            checked = numerical_check(branch['native'],branch['length'],sums,available_genes)
            zero_evidence += not checked['q_available']
            values = {k:(v if math.isfinite(v) else None) for k,v in branch['native'].items()}
            branch_rows.append(dict(split=branch['split'],canonical_side=branch['canonical_side'],
                incident_groups=branch['groups'],native=values,native_map_length=branch['length'],
                gene_dispositions=dict(gene_counts),**checked))
    assert states == (len(tips)-3)*125
    write(folder/'branches.json',branch_rows)
    result = dict(case=case['case'],markers=125,taxa=len(tips),roles=case['roles'],
        internal_branches=len(branches),branch_gene_states=states,zero_evidence_branches=zero_evidence,
        independent_global_resolved_denominator=global_denominator,
        native_global_quartet_score=score,native_normalized_score=normalized_score,
        global_numerator_independently_verified=False,
        native_receipt_sha256=sha(rp),elapsed_seconds=time.monotonic()-started,
        artifacts={p.name:sha(p) for p in folder.iterdir() if p.is_file()},
        scope='All local branch/gene quartet counts independently computed with integer LCA-distribution DP; per-gene available-quartet normalization retains missing and unresolved evidence. Native concordant evidence and unordered NNI alternatives, effective genes, fractions, default-prior local posteriors and MAP coalescent lengths checked. Full global resolved quartet denominator independently counted, normalized-score arithmetic checked; native global matching numerator and constrained search optimum not independently established. Roots/model/gene uncertainty and biological interpretation remain unqualified.')
    write(folder/'receipt.json',result)
    for p in folder.iterdir(): bind(bindings,p)
    print('full_case_independent_coalescent_quartets_passed',case['case'],states,flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    native_plan,cases,bindings = source_bundle(plan,args.plan)
    out = Path(plan['output']); out.mkdir(parents=True,exist_ok=False)
    write(out/'config.json',dict(plan_sha256=sha(args.plan),source_hashes=bindings))
    bind(bindings,out/'config.json')
    summaries = []
    for case in cases:
        summaries.append(audit_case(native_plan,plan['native_plan'],case,out,bindings))
    verify(bindings)
    result = dict(status=('passed_full_30_coalescent_local_quartet_and_numeric_readback'
        if plan['mode'] == 'full_batch' else 'passed_complete_named_native_case_coalescent_quartet_preflight'),
        plan_sha256=sha(args.plan),native_plan_sha256=sha(plan['native_plan']),cases=len(cases),
        markers_per_case=125,internal_branches=sum(c['internal_branches'] for c in summaries),
        branch_gene_states=sum(c['branch_gene_states'] for c in summaries),summaries=summaries,
        source_hashes=bindings,scientific_eligibility=False,
        global_numerator_independently_verified=False,scope=plan['scope'])
    if plan['mode'] == 'full_batch':
        assert result['cases'] == 30 and result['internal_branches'] == 15510
        assert result['branch_gene_states'] == 1938750
    write(out/'receipt.json',result)
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,float,bool))}),flush=True)


if __name__ == '__main__': main()
