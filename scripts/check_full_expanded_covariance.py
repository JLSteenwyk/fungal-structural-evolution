#!/usr/bin/env python3
"""Check alias/reuse covariance contracts using synthetic prior handoffs.

Prior case closure, tree audits and journals are synthetic contracts. Tests
exercise the software, not a reduced biological experiment or source acceptance.
"""
import argparse
import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import dendropy
import numpy as np
from scipy import sparse
from check_full_matched_coverage_cases import node
from index_full_matching_cases import metadata
from full_matching_case_sources import CASE_FIELDS as ORIGINAL_FIELDS
from full_expanded_covariance_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True); Path(path).write_text(json.dumps(value, indent=2) + '\n')


def table(path, rows, fields):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(rows)


def invoke(command, success=True):
    result = subprocess.run(command, capture_output=True, text=True)
    if success and result.returncode: raise RuntimeError(result.stderr + result.stdout)
    if not success: assert result.returncode != 0, 'False expanded covariance export accepted'


def setup(root):
    graph = root / 'graph'; graph.mkdir(); case_root = root / 'cases'; case_root.mkdir(); kernel = root / 'kernel'; kernel.mkdir()
    targets, backgrounds = {}, {}
    specs = [('A', 'B', 'C', 'D', 'FA', 'profile'), ('A', 'E', 'F', 'G', 'FB', 'profile'),
        ('H', 'H', 'I', 'J', 'FC', 'mafft'), ('H', 'L', 'M', 'N', 'FD', 'mafft'),
        ('A', 'B', 'C', 'D', 'FA', 'mafft'), ('U', 'V', 'W', 'X', 'FE', 'profile')]
    for i, (a, b, c, d, family, guide) in enumerate(specs, 1):
        t = node('T' + str(i), guide, a, b, True); bg = node('B' + str(i), guide, c, d)
        t['family'] = bg['family'] = family
        bg['taxon_a'] = 'taxon-other'; bg['taxon_b'] = 'taxon-mafft'
        if i == 3: bg['taxon_a'] = bg['taxon_b'] = t['taxon_id']
        if i == 4: t['gene_a'] = targets['T3']['gene_b']
        for n in [t, bg]:
            for end in ['a', 'b']: n['sha256_' + end] = hashlib.sha256((n['model_id_' + end] + str(n['version_' + end])).encode()).hexdigest()
        targets[t['node_id']] = t; backgrounds[bg['node_id']] = bg
    for old, new in [('B1', 'B7'), ('B2', 'B8')]:
        bg = copy.deepcopy(backgrounds[old]); bg['node_id'] = new; bg['gene_a'] = new + '-a'; bg['gene_b'] = new + '-b'; backgrounds[new] = bg
    for role, records in [('target', targets), ('background', backgrounds)]:
        (graph / (role + '_nodes.jsonl')).write_text(''.join(json.dumps(r) + '\n' for r in records.values()))
    case_rows = []
    for tid, bid in [(f'T{i}', f'B{i}') for i in range(1, 7)] + [('T1', 'B7'), ('T2', 'B8')]:
        r = metadata(targets[tid], backgrounds[bid]); r.update(selection_records=3, endpoint_order_bits=3, policy_bits=3, scenario_bits=3)
        r.update({k: 0 for k in ORIGINAL_FIELDS if k.endswith('_pass_bits')}); case_rows.append(r)
    case_rows.sort(key=lambda r: (r['target_id'], r['background_id'])); table(case_root / 'case_index.tsv.gz', case_rows, ORIGINAL_FIELDS)
    matching = root / 'matching-plan.json'; write(matching, dict(graph=str(graph)))
    case_plan = root / 'case-plan.json'; write(case_plan, dict(output=str(case_root), matching_plan=str(matching)))
    summary = dict(logical_cases=len(case_rows), physical_cases=len({r['physical_case_id'] for r in case_rows}), selected_records=24)
    rp = case_root / 'receipt.json'; write(rp, dict(status='complete_full_matching_case_index_pending_independent_readback',
        plan_sha256=sha(case_plan), **summary, artifacts={'case_index.tsv.gz': sha(case_root / 'case_index.tsv.gz')}, scientific_eligibility=False))
    bindings = {str(p): sha(p) for p in [rp, case_plan, matching, case_root / 'case_index.tsv.gz', graph / 'target_nodes.jsonl', graph / 'background_nodes.jsonl']}
    archive = root / 'case-archive.json'; write(archive, dict(status='complete_verified_full_matching_logical_case_index_archive', source_hashes=bindings, services=[{}, {}], summary=summary))
    completion = root / 'case-completion.json'; write(completion, dict(status='complete_verified_full_matching_logical_case_index',
        exact_process_journals_checked=2, bound_source_hashes=len(bindings), full_hash_archive=str(archive), full_hash_archive_sha256=sha(archive),
        **summary, producer_receipt=str(rp), producer_receipt_sha256=sha(rp), scientific_eligibility=False))
    taxa = ['taxon-mafft', 'taxon-other', 'taxon-profile']; write(kernel / 'taxa.json', taxa)
    trees = {}
    for label, newick in [('first', '((taxon-mafft:0.4,taxon-profile:0.7):0.3,taxon-other:0.6);'),
        ('second', '((taxon-mafft:0.8,taxon-other:0.2):0.1,taxon-profile:0.5);')]:
        tp = root / (label + '.tree'); tp.write_text(newick)
        ap = root / (label + '-audit.json'); write(ap, dict(status='fixture_tree_passed'))
        trees[label] = dict(tree=str(tp), audit=str(ap), audit_status='fixture_tree_passed')
        tree = dendropy.Tree.get(path=str(tp), schema='newick', preserve_underscores=True); distances = tree.phylogenetic_distance_matrix()
        tips = {leaf.taxon.label: leaf.taxon for leaf in tree.leaf_node_iter()}
        d = np.array([[distances(tips[a], tips[b]) for b in taxa] for a in taxa]); center = np.eye(3) - np.ones((3, 3)) / 3
        np.savez_compressed(kernel / (label + '.npz'), distances=d, centered_kernel=-.5 * center @ d @ center)
    kernel_plan = root / 'kernel-plan.json'; write(kernel_plan, dict(trees=trees, pins={str(p): sha(p) for spec in trees.values() for p in [spec['tree'], spec['audit']]}))
    krp = kernel / 'receipt.json'; write(krp, dict(status='complete_species_distance_kernels_pending_independent_readback', plan_sha256=sha(kernel_plan),
        trees={t: dict(tree_sha256=sha(s['tree']), audit_sha256=sha(s['audit'])) for t, s in trees.items()},
        artifacts={p.name: sha(p) for p in kernel.iterdir()}))
    audit = root / 'kernel-readback.json'; write(audit, dict(status='passed_full_species_distance_kernel_readback', source_receipt_sha256=sha(krp)))
    patterns = set()
    for r in case_rows:
        w = {}
        for t, c in [(r['focal_taxon'], 2), (r['background_taxon_a'], -1), (r['background_taxon_b'], -1)]: w[t] = w.get(t, 0) + c
        patterns.add(tuple(sorted((t, c) for t, c in w.items() if c)))
    pp = root / 'plan.json'; write(pp, dict(case_completion=str(completion), case_plan=str(case_plan), kernel_plan=str(kernel_plan), kernel_root=str(kernel),
        kernel_readback=str(audit), trees=list(trees), expected=dict(**summary, species_columns=3, patterns=len(patterns)),
        output=str(root / 'output'), pins={}, covariance_block_rows=2, resources=dict(minimum_free_disk_gib=0), scope=__doc__))
    return pp


def run():
    with tempfile.TemporaryDirectory(prefix='expanded-covariance-fixture-') as temp:
        root = Path(temp); pp = setup(root); out = root / 'output'
        producer = [sys.executable, 'scripts/prepare_full_expanded_covariance.py', '--plan', str(pp)]
        reader = [sys.executable, 'scripts/readback_full_expanded_covariance.py', '--plan', str(pp), '--output']
        invoke(producer); invoke(reader + [str(root / 'passed.json')]); invoke(producer, False)
        receipt = json.loads((out / 'receipt.json').read_text()); assert receipt['logical_cases'] == 8 and receipt['rank'] == 2 and receipt['zero_patterns'] == 1
        assert receipt['family_components'] == 3 and receipt['families'] == 5 and receipt['entity_occurrences'] == 112
        saved = {n: (out / n).read_bytes() for n in receipt['artifacts']}; saved_receipt = (out / 'receipt.json').read_bytes()
        labels = ['missing_case', 'duplicate_case', 'wrong_case_pattern', 'wrong_taxon_index', 'wrong_reuse',
            'missing_occurrence', 'duplicate_occurrence', 'wrong_signed_loading', 'wrong_unsigned_loading', 'collapsed_gene',
            'ignored_model_version', 'wrong_coordinate_hash', 'wrong_component', 'wrong_component_size',
            'changed_weights', 'changed_design', 'changed_basis', 'changed_projection', 'changed_factor', 'changed_core', 'changed_cholesky', 'changed_census']
        def rows(name):
            opener = gzip.open if name.endswith('.gz') else open
            with opener(out / name, 'rt') as h:
                reader = csv.DictReader(h, delimiter='\t'); return list(reader), reader.fieldnames
        for label in labels:
            for n, content in saved.items(): (out / n).write_bytes(content)
            r = copy.deepcopy(receipt)
            if label in ['missing_case', 'duplicate_case', 'wrong_case_pattern', 'wrong_taxon_index', 'wrong_reuse']:
                name = 'case_covariance_index.tsv.gz'; records, fields = rows(name)
                if label == 'missing_case': records.pop()
                elif label == 'duplicate_case': records.append(dict(records[0]))
                elif label == 'wrong_case_pattern': records[0]['species_pattern_id'] = '0' * 64
                elif label == 'wrong_taxon_index': records[0]['focal_index'] = '99'
                else: records[0]['selection_records'] = '99'
                table(out / name, records, fields)
            elif label in ['missing_occurrence', 'duplicate_occurrence', 'wrong_signed_loading', 'wrong_unsigned_loading', 'collapsed_gene', 'ignored_model_version', 'wrong_coordinate_hash']:
                name = 'entity_incidence.tsv.gz'; records, fields = rows(name)
                if label == 'missing_occurrence': records.pop()
                elif label == 'duplicate_occurrence': records.append(dict(records[0]))
                elif label in ['wrong_signed_loading', 'wrong_unsigned_loading']: records[0]['signed_loading' if 'signed_' in label and 'unsigned_' not in label else 'unsigned_loading'] = '0'
                elif label == 'collapsed_gene': next(r for r in records if r['entity_kind'] == 'gene')['entity_id'] = '0' * 64
                else:
                    row = next(r for r in records if r['entity_kind'] == 'model'); token = json.loads(row['entity_value']); token[1 if label == 'ignored_model_version' else 2] = 999 if label == 'ignored_model_version' else '0' * 64; row['entity_value'] = json.dumps(token)
                table(out / name, records, fields)
            elif label in ['wrong_component', 'wrong_component_size']:
                name = 'family_components.tsv'; records, fields = rows(name); records[0]['family_component' if label == 'wrong_component' else 'component_families'] = 'bad'; table(out / name, records, fields)
            elif label == 'changed_weights':
                name = 'patterns.tsv'; records, fields = rows(name); records[0]['twice_weights'] = '[2, -1]'; table(out / name, records, fields)
            elif label == 'changed_design':
                name = 'species_contrast_design.npz'; array = sparse.load_npz(out / name).toarray(); array[0, 0] += .5; sparse.save_npz(out / name, sparse.csr_matrix(array))
            elif label in ['changed_basis', 'changed_projection']:
                name = 'pattern_basis.npy' if label == 'changed_basis' else 'species_projection.npy'; array = np.load(out / name); array.flat[0] += .1; np.save(out / name, array)
            elif label in ['changed_factor', 'changed_core', 'changed_cholesky']:
                name = 'first.npz'
                with np.load(out / name) as z: arrays = {n: z[n] for n in z.files}
                arrays[{'changed_factor': 'factor', 'changed_core': 'reduced_covariance', 'changed_cholesky': 'cholesky'}[label]].flat[0] += .1; np.savez_compressed(out / name, **arrays)
            else: r['selected_records'] += 1
            r['artifacts'] = {n: sha(out / n) for n in r['artifacts']}; write(out / 'receipt.json', r)
            invoke(reader + [str(root / (label + '.json'))], False)
        for n, content in saved.items(): (out / n).write_bytes(content)
        (out / 'receipt.json').write_bytes(saved_receipt)
        (out / 'receipt.json').unlink(); (out / 'entity_incidence.tsv.gz.tmp').write_text('interrupted')
        invoke(producer); invoke(reader + [str(root / 'replayed.json')]); replay = json.loads((out / 'receipt.json').read_text())
        assert all(receipt[k] == replay[k] for k in SUMMARY_FIELDS)
        return dict(status='passed_full_expanded_covariance_software_contracts', logical_cases=8, patterns=receipt['patterns'], rank=2,
            entity_occurrences=112, zero_pattern_and_same_model_preserved=True, cross_family_model_aliases_preserved=True,
            gene_aliases_and_distinct_genes_preserved=True, completed_restart_refused=True, interrupted_full_replay_passed=True,
            rejected_rehashed_exports=labels, synthetic_prior_case_tree_journal_contracts=True, production_sources_or_journals_tested=False, scope=__doc__)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True); a = p.parse_args(); r = run()
    r['source_hashes'] = {p: sha(p) for p in [__file__, 'scripts/full_expanded_covariance_sources.py', 'scripts/prepare_full_expanded_covariance.py', 'scripts/readback_full_expanded_covariance.py']}
    with a.output.open('x') as handle: handle.write(json.dumps(r, indent=2) + '\n')
    print(json.dumps(r, indent=2))
