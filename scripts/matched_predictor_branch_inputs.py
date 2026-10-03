"""Identical-protein/position predictor controls on the complete candidate tree grid.

This prepares conditional comparisons, not accepted species trees or evolutionary
rates. Computational reuse retains every marker/view and original branch mapping.
"""
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path

import dendropy
import numpy as np

from audit_selected_taxon_identity_snapshot_v2 import rows, sha
from compare_paired_source_states import load as load_alignment
from full_expanded_model_design_sources import digest
from structural_marker_tree_projection import load as load_sources, SOURCES, normalize, mask_for

SCHEMA = 'matched-predictor-branch-inputs-20261003-v1'
SUMMARY_FIELDS = ['tree_views', 'marker_slots', 'taxon_entries', 'comparison_cases',
                  'case_status_counts', 'unique_ready_inputs', 'future_native_roles',
                  'matched_input_taxa', 'matched_input_markers', 'original_internal_branch_cells',
                  'branch_status_counts', 'unique_ready_internal_splits', 'maximum_taxa', 'maximum_columns']
BRANCH_STATUSES = ['insufficient_matched_taxa', 'no_observations_on_one_side',
                   'terminal_projection', 'unique_internal_projection', 'shared_internal_projection']


def verify(bindings):
    for path, expected in bindings.items():
        assert sha(path) == expected, path


def load():
    views, positions, markers, datasets, taxa, bindings = load_sources()
    context_root = Path('results/phylogeny/paired-source-model-context-20260926-v1')
    cp = context_root/'receipt.json'
    proof_path = Path('metadata/paired_source_model_context_completed_20260926.json')
    proof = json.loads(proof_path.read_text()); receipt = json.loads(cp.read_text())
    assert proof['status'] == 'passed_paired_source_model_context_readback'
    assert proof['producer_receipt'] == str(cp) and proof['producer_receipt_sha256'] == sha(cp)
    assert receipt['status'] == 'complete_paired_source_full_sequence_context'
    assert proof['cells'] == receipt['cells'] == 673
    assert proof['complete_sequence_equal_cells'] == 673
    for path, value in receipt['source_sha256'].items():
        assert path not in bindings or bindings[path] == value
        bindings[path] = value
    for name, value in receipt['artifacts'].items():
        bindings[str(context_root/name)] = value
    for path in [cp, proof_path]: bindings[str(path)] = sha(path)
    contexts = {(r['marker'], r['taxon']):r for r in rows(context_root/'source_model_context.tsv')}
    assert len(contexts) == 673
    expected = {(m,t) for m in markers for t in
                datasets['AlphaFold']['usable'][m] & datasets['ESMFold']['usable'][m]}
    assert set(contexts) == expected
    sequences = {}
    for row in contexts.values():
        assert row['full_sequence_status'] == 'identical_complete_encoded_sequence'
        assert row['aa_mismatches'] == '0'
        full = []
        for prefix in ['reference', 'local']:
            path = row[prefix+'_encoding_path']
            assert bindings[path] == row[prefix+'_encoding_sha256']
            assert bindings[row[prefix+'_model_path']] == row[prefix+'_model_sha256']
            if path not in sequences:
                with np.load(path, allow_pickle=False) as a: sequences[path] = str(a['sequence'])
            seq = sequences[path]
            assert len(seq) == int(row[prefix+'_length'])
            assert hashlib.sha256(seq.encode()).hexdigest() == row[prefix+'_sequence_sha256']
            full.append(seq)
        assert full[0] == full[1]
    alignments = {}; required = {}
    for marker in markers:
        overlap = {t for m,t in contexts if m == marker}
        if not overlap:
            alignments[marker] = None; required[marker] = None; continue
        pair = [load_alignment(Path(SOURCES[s][0]), marker) for s in SOURCES]
        coverage = [next(r for r in rows(Path(SOURCES[s][0])/'marker_summary.tsv')
                         if r['marker'] == marker) for s in SOURCES]
        original = [int(r['original_marker_columns']) for r in coverage]
        assert original[0] == original[1]
        required[marker] = max(50, math.ceil(.3*original[0]))
        assert overlap == set(pair[0][1]['aa']) & set(pair[1][1]['aa'])
        alignments[marker] = pair
    verify(bindings)
    return dict(views=views, positions=positions, markers=markers, taxa=taxa,
                contexts=contexts, alignments=alignments, required=required), bindings


def joint_rows(left, right, contexts, required):
    """Build the shared observation mask before selecting any tree membership."""
    pa, da = left; pb, db = right
    ia = {p:i for i,p in enumerate(pa)}; ib = {p:i for i,p in enumerate(pb)}
    columns = sorted(set(pa) & set(pb)); result = {}; accounting = []
    for taxon in sorted(set(da['aa']) & set(db['aa'])):
        assert taxon in contexts
        aa, af, esm = [], [], []
        for p in columns:
            i,j = ia[p], ib[p]; a,b = da['aa'][taxon][i],db['aa'][taxon][j]
            if a == '?' or b == '?':
                aa.append('?'); af.append('?'); esm.append('?')
            else:
                assert a == b, (taxon,p,a,b)
                assert da['3di'][taxon][i] != '?' and db['3di'][taxon][j] != '?'
                aa.append(a); af.append(da['3di'][taxon][i]); esm.append(db['3di'][taxon][j])
        observed = sum(c != '?' for c in aa)
        ctx = contexts[taxon]
        assert observed == int(ctx['jointly_observed_positions']) == int(ctx['same_aa_positions'])
        assert sum(a != b for a,b in zip(af,esm)) == int(ctx['state_mismatches'])
        accounting.append(dict(taxon=taxon, jointly_observed=observed,
                               required_observed=required, eligible=observed >= required))
        if observed >= required:
            result[taxon] = dict(aa=''.join(aa), AlphaFold=''.join(af), ESMFold=''.join(esm))
    return columns, result, accounting


def select(columns, data, membership):
    taxa = sorted(set(data) & set(membership))
    retained = [i for i in range(len(columns)) if any(data[t]['aa'][i] != '?' for t in taxa)]
    selected = {label:{t:''.join(data[t][label][i] for i in retained) for t in taxa}
                for label in ['aa','AlphaFold','ESMFold']}
    return [columns[i] for i in retained], taxa, selected


def branch_mapping(view, taxa, positions):
    retained = mask_for(taxa, positions); groups = defaultdict(list); codes = []
    for i, branch in enumerate(view['branches']):
        side = int(branch['split_mask_hex'],16) & retained
        a,b = side.bit_count(), len(taxa)-side.bit_count()
        if len(taxa) < 4: code = 0
        elif min(a,b) == 0: code = 1
        elif min(a,b) == 1: code = 2
        else:
            code = -1; groups[normalize(side,retained)].append(i)
        codes.append(code)
    for indices in groups.values():
        for i in indices: codes[i] = 3 if len(indices) == 1 else 4
    if len(taxa) >= 4: assert len(groups) == len(taxa)-3
    mapping = [dict(split_mask_hex=hex(mask), original_branch_indices=indices)
               for mask,indices in sorted(groups.items())]
    counts = Counter(codes)
    return mapping, {label:counts[i] for i,label in enumerate(BRANCH_STATUSES)}


def topology(taxa, mapping, positions):
    """Build only the declared compatible unrooted splits with uniform starts."""
    namespace = dendropy.TaxonNamespace(taxa); bits = []
    for item in mapping:
        global_mask = int(item['split_mask_hex'],16)
        bits.append(sum(1 << i for i,t in enumerate(taxa) if global_mask & (1 << positions[t])))
    tree = dendropy.Tree.from_split_bitmasks(bits, namespace, is_rooted=False)
    for node in tree.preorder_node_iter():
        node.label = None
        node.edge.length = None if node is tree.seed_node else .1
    tree.encode_bipartitions()
    assert len([e for e in tree.internal_edges() if e.tail_node is not None]) == len(taxa)-3
    return tree.as_string(schema='newick', suppress_rooting=True, suppress_annotations=True)


def fasta_text(data):
    return ''.join('>'+t+'\n'+data[t]+'\n' for t in sorted(data))


def case(source, view_index, marker, columns, data):
    view = source['views'][view_index]
    selected_columns, taxa, selected = select(columns, data, view['taxa'])
    mapping, counts = branch_mapping(view, taxa, source['positions'])
    status = 'ready_for_matched_predictor_branch_fits' if len(taxa) >= 4 else 'insufficient_joint_observation_taxa'
    identity = dict(schema=SCHEMA, marker=marker, columns=selected_columns, taxa=taxa,
                    alignments=selected, internal_splits=[r['split_mask_hex'] for r in mapping])
    input_id = digest(identity) if len(taxa) >= 4 else None
    record = dict(view_index=view_index, cohort=view['cohort'], view=view['view'], marker=marker,
                  status=status, input_id=input_id, taxa=taxa, columns=selected_columns,
                  required_observed=source['required'][marker], original_internal_branches=len(view['branches']),
                  branch_status_counts=counts, internal_branch_mapping=mapping,
                  alignment_sha256={k:hashlib.sha256(fasta_text(v).encode()).hexdigest() for k,v in selected.items()},
                  scientific_eligibility=False)
    return record, selected


def summary(source, records, inputs):
    counts = Counter(r['status'] for r in records); branches = Counter()
    for record in records: branches.update(record['branch_status_counts'])
    return dict(tree_views=len(source['views']), marker_slots=len(source['markers']), taxon_entries=len(source['positions']),
                comparison_cases=len(records), case_status_counts=dict(counts), unique_ready_inputs=len(inputs),
                future_native_roles=7*len(inputs), matched_input_taxa=len(set().union(*(set(r['taxa']) for r in inputs.values()))) if inputs else 0,
                matched_input_markers=len({r['marker'] for r in inputs.values()}),
                original_internal_branch_cells=sum(r['original_internal_branches'] for r in records),
                branch_status_counts=dict(branches), unique_ready_internal_splits=sum(len(r['internal_branch_mapping']) for r in inputs.values()),
                maximum_taxa=max((len(r['taxa']) for r in inputs.values()),default=0),
                maximum_columns=max((len(r['columns']) for r in inputs.values()),default=0))
