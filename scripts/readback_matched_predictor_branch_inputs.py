#!/usr/bin/env python3
"""Rebuild every matched control from raw sources; independently prune raw trees."""
import argparse
from collections import Counter, defaultdict
import copy
import gzip
import hashlib
import json
from pathlib import Path

from Bio import Phylo

from audit_selected_taxon_identity_snapshot_v2 import sha
from full_expanded_model_design_sources import digest
from matched_predictor_branch_inputs import load, verify, SCHEMA, SUMMARY_FIELDS, BRANCH_STATUSES
from prepare_matched_predictor_branch_inputs import STATUS as PRODUCER_STATUS

STATUS = 'passed_full_matched_predictor_branch_inputs_independent_readback'


def raw_fasta(path):
    result = {}; name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:]; assert name and name not in result and not any(c.isspace() for c in name)
            result[name] = ''
        else:
            assert name is not None
            result[name] += line.strip()
    assert result and len({len(s) for s in result.values()}) == 1
    return result


def raw_splits(tree, taxa, positions):
    leaves = [n.name for n in tree.get_terminals()]
    assert len(leaves) == len(taxa) and set(leaves) == set(taxa)
    mask = sum(2**positions[t] for t in taxa); splits = set()
    for clade in tree.find_clades():
        if clade is tree.root or clade.is_terminal(): continue
        side = sum(2**positions[t.name] for t in clade.get_terminals())
        other = mask ^ side
        if min(side.bit_count(),other.bit_count()) >= 2:
            splits.add(min([side,other],key=lambda x:(x.bit_count(),x)))
    assert len(splits) == len(taxa)-3
    return splits


def independent_marker(source, marker):
    pair = source['alignments'][marker]
    if pair is None: return [], {}, []
    left = {p:i for i,p in enumerate(pair[0][0])}; right = {p:i for i,p in enumerate(pair[1][0])}
    original = sorted(set(left) & set(right)); data = {}; account = []
    common_taxa = sorted(set(pair[0][1]['aa']) & set(pair[1][1]['aa']))
    for taxon in common_taxa:
        aa = {}; af = {}; esm = {}
        for position in original:
            a = pair[0][1]['aa'][taxon][left[position]]
            b = pair[1][1]['aa'][taxon][right[position]]
            if a != '?' and b != '?':
                assert a == b
                aa[position] = a
                af[position] = pair[0][1]['3di'][taxon][left[position]]
                esm[position] = pair[1][1]['3di'][taxon][right[position]]
                assert af[position] != '?' and esm[position] != '?'
        row = source['contexts'][marker,taxon]
        assert len(aa) == int(row['same_aa_positions']) == int(row['jointly_observed_positions'])
        assert sum(af[p] != esm[p] for p in aa) == int(row['state_mismatches'])
        eligible = len(aa) >= source['required'][marker]
        account.append(dict(marker=marker,taxon=taxon,jointly_observed=len(aa),
                            required_observed=source['required'][marker],eligible=eligible))
        if eligible: data[taxon] = dict(aa=aa,AlphaFold=af,ESMFold=esm)
    return original, data, account


def expected_case(source, vi, marker, data):
    view = source['views'][vi]; taxa = sorted(set(data) & set(view['taxa']))
    columns = sorted(set().union(*(set(data[t]['aa']) for t in taxa))) if taxa else []
    alignments = {label:{t:''.join(data[t][label].get(p,'?') for p in columns) for t in taxa}
                  for label in ['aa','AlphaFold','ESMFold']}
    observed_mask = sum(2**source['positions'][t] for t in taxa)
    mapping = defaultdict(list); counts = Counter()
    for index, branch in enumerate(view['branches']):
        side = int(branch['split_mask_hex'],16) & observed_mask
        complement = observed_mask ^ side
        if len(taxa) < 4: counts[0] += 1
        elif min(side.bit_count(),complement.bit_count()) == 0: counts[1] += 1
        elif min(side.bit_count(),complement.bit_count()) == 1: counts[2] += 1
        else:
            key = min([side,complement],key=lambda x:(x.bit_count(),x))
            mapping[key].append(index)
    for indices in mapping.values(): counts[3 if len(indices) == 1 else 4] += len(indices)
    internal = [dict(split_mask_hex=hex(k),original_branch_indices=v) for k,v in sorted(mapping.items())]
    identity = dict(schema=SCHEMA,marker=marker,columns=columns,taxa=taxa,
                    alignments=alignments,internal_splits=[r['split_mask_hex'] for r in internal])
    ready = len(taxa) >= 4
    checksum = {}
    for label, rows in alignments.items():
        encoded = ''.join('>'+t+'\n'+rows[t]+'\n' for t in taxa).encode()
        checksum[label] = hashlib.sha256(encoded).hexdigest()
    expected = dict(view_index=vi,cohort=view['cohort'],view=view['view'],marker=marker,
        status='ready_for_matched_predictor_branch_fits' if ready else 'insufficient_joint_observation_taxa',
        input_id=digest(identity) if ready else None,taxa=taxa,columns=columns,
        required_observed=source['required'][marker],original_internal_branches=len(view['branches']),
        branch_status_counts={label:counts[i] for i,label in enumerate(BRANCH_STATUSES)},
        internal_branch_mapping=internal,alignment_sha256=checksum,scientific_eligibility=False)
    return expected, alignments, set(mapping)


def check_ready_input(root, key, record, alignments, actual_splits, positions):
    folder = root/'inputs'/key
    assert {p.name for p in folder.iterdir()} == {'aa.faa','AlphaFold.faa','ESMFold.faa','config.json','topology.nwk'}
    for label, expected in alignments.items():
        assert raw_fasta(folder/(label+'.faa')) == expected
        assert sha(folder/(label+'.faa')) == record['alignment_sha256'][label]
    config = json.loads((folder/'config.json').read_text())
    expected = {k:record[k] for k in ['input_id','marker','taxa','columns','required_observed','alignment_sha256']}
    expected.update(internal_splits=[r['split_mask_hex'] for r in record['internal_branch_mapping']],
                    topology_sha256=sha(folder/'topology.nwk'),
                    branch_initialization='all_edges_0.1_expected_state_substitutions_per_site_start_only')
    assert config == expected
    tree = Phylo.read(folder/'topology.nwk','newick')
    assert raw_splits(tree,record['taxa'],positions) == actual_splits
    for clade in tree.find_clades():
        if clade is not tree.root:
            assert clade.branch_length == .1 and clade.confidence is None


def check_grid(source, root, verbose=True):
    axes = json.loads((root/'input_axes.json').read_text())
    assert axes == dict(schema=SCHEMA,views=source['views'],positions=source['positions'],markers=source['markers'],taxa=source['taxa'])
    with gzip.open(root/'comparison_cases.jsonl.gz','rt') as f: records = [json.loads(line) for line in f]
    assert len(records) == len(source['markers'])*len(source['views'])
    index = {(r['marker'],r['view_index']):r for r in records}; assert len(index) == len(records)
    inputs = {}; accounts = []; raw_trees = {}; prunings = {}; counts = Counter(); branches = Counter()
    for marker in source['markers']:
        _, data, account = independent_marker(source,marker); accounts.extend(account)
        for vi, view in enumerate(source['views']):
            expected, alignments, projected_splits = expected_case(source,vi,marker,data)
            record = index[marker,vi]; assert record == expected, (marker,vi)
            counts[record['status']] += 1; branches.update(record['branch_status_counts'])
            key = record['input_id']
            if key is None: continue
            source_tree = view['source_tree']; pruning_key = source_tree, tuple(record['taxa'])
            if pruning_key not in prunings:
                if source_tree not in raw_trees: raw_trees[source_tree] = Phylo.read(source_tree,'newick')
                tree = copy.deepcopy(raw_trees[source_tree]); wanted = set(record['taxa'])
                assert wanted <= {t.name for t in tree.get_terminals()}
                for tip in list(tree.get_terminals()):
                    if tip.name not in wanted: tree.prune(tip)
                prunings[pruning_key] = raw_splits(tree,record['taxa'],source['positions'])
            assert prunings[pruning_key] == projected_splits
            if key not in inputs:
                check_ready_input(root,key,record,alignments,projected_splits,source['positions'])
                inputs[key] = record
        if verbose: print('independent_matched_predictor_marker_verified',marker,'cases',len(counts),'unique_inputs',len(inputs),flush=True)
    assert json.loads((root/'joint_taxon_observation_accounting.json').read_text()) == accounts
    manifest = [dict(input_id=k,marker=r['marker'],taxa=len(r['taxa']),columns=len(r['columns']),
                     internal_splits=len(r['internal_branch_mapping'])) for k,r in sorted(inputs.items())]
    assert json.loads((root/'input_manifest.json').read_text()) == manifest
    assert {p.name for p in (root/'inputs').iterdir()} == set(inputs)
    return dict(tree_views=len(source['views']),marker_slots=len(source['markers']),taxon_entries=len(source['positions']),
        comparison_cases=len(records),case_status_counts=dict(counts),unique_ready_inputs=len(inputs),future_native_roles=7*len(inputs),
        matched_input_taxa=len({t for r in inputs.values() for t in r['taxa']}),matched_input_markers=len({r['marker'] for r in inputs.values()}),
        original_internal_branch_cells=sum(r['original_internal_branches'] for r in records),branch_status_counts=dict(branches),
        unique_ready_internal_splits=sum(len(r['internal_branch_mapping']) for r in inputs.values()),
        maximum_taxa=max((len(r['taxa']) for r in inputs.values()),default=0),maximum_columns=max((len(r['columns']) for r in inputs.values()),default=0))


def run(path):
    path = Path(path); plan = json.loads(path.read_text()); verify(plan['pins']); root = Path(plan['output'])
    rp = root/'receipt.json'; producer = json.loads(rp.read_text())
    assert producer['status'] == PRODUCER_STATUS and producer['plan_sha256'] == sha(path)
    verify(producer['source_hashes'])
    assert all(sha(root/p) == h for p,h in producer['artifacts'].items())
    source, bindings = load(); bindings.update(plan['pins']); bindings[str(path)] = sha(path)
    result = check_grid(source,root)
    assert all(producer[k] == result[k] for k in SUMMARY_FIELDS)
    verify(producer['source_hashes']); verify(bindings)
    assert all(sha(root/p) == h for p,h in producer['artifacts'].items())
    bindings[str(rp)] = sha(rp)
    receipt = dict(status=STATUS,plan_sha256=sha(path),producer_receipt_sha256=sha(rp),**result,
                   source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with (root/'readback.json').open('x') as f:f.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
