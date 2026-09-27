#!/usr/bin/env python3
"""Preserve endpoint taxa and reuse counts for every selected comparison."""
import csv
import gzip
import json
from collections import Counter
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def main():
    graph = Path('results/orthology/background-match-graph-20260927-v1')
    selection = Path('results/orthology/background-control-selection-20260927-v1')
    dependency = Path('results/orthology/selected-control-family-dependencies-20260927-v1')
    kernel = Path('results/phylogeny/species-distance-kernels-20260927-v1')
    out = Path('results/orthology/selected-taxon-inputs-20260927-v1')
    roots = [graph, selection, dependency, kernel]
    pins = {}
    for root in roots:
        receipt = root / 'receipt.json'
        pins[str(receipt)] = sha(receipt)
        for name, digest in json.loads(receipt.read_text())['artifacts'].items():
            path = root / name
            assert sha(path) == digest, path
            pins[str(path)] = digest
    for name in ['background_match_graph_completed_readback_20260927.json',
                 'background_control_selection_completed_readback_20260927.json',
                 'selected_control_family_dependency_readback_20260927.json',
                 'species_distance_kernel_readback_20260927.json']:
        path = Path('metadata') / name
        assert json.loads(path.read_text())['status'].startswith('passed'), path
        pins[str(path)] = sha(path)
    taxa = json.loads((kernel / 'taxa.json').read_text())
    index = {taxon: i for i, taxon in enumerate(taxa)}
    assert len(index) == 526
    components = {(r['guide'], r['family']): r['component_id']
                  for r in csv.DictReader((dependency / 'family_components.tsv').open(), delimiter='\t')
                  if r['entity_kind'] == 'combined'}
    nodes = {}
    for role in ['target', 'background']:
        nodes[role] = {r['node_id']: r for r in map(json.loads, (graph / (role + '_nodes.jsonl')).open())}
    usage = {'target': Counter(), 'background': Counter()}
    records = 0
    focal_present = 0
    with gzip.open(selection / 'selections.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            target = nodes['target'][row['target_id']]
            background = nodes['background'][row['background_id']]
            assert (target['guide'], target['family']) == (background['guide'], background['family'])
            focal_present += target['taxon_id'] in [background['taxon_a'], background['taxon_b']]
            for role in usage:
                usage[role][row[role + '_id']] += 1
            records += 1
    assert records == json.loads((selection / 'receipt.json').read_text())['selected_records']
    out.mkdir(exist_ok=False)
    fields = ['role', 'node_id', 'guide', 'family', 'family_component', 'pair_key',
              'taxon_a', 'taxon_b', 'kernel_index_a', 'kernel_index_b', 'selected_record_uses']
    with (out / 'selected_node_taxa.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for role in usage:
            for node_id, count in sorted(usage[role].items()):
                node = nodes[role][node_id]
                a = node['taxon_id'] if role == 'target' else node['taxon_a']
                b = node['taxon_id'] if role == 'target' else node['taxon_b']
                writer.writerow(dict(role=role, node_id=node_id, guide=node['guide'],
                                     family=node['family'], family_component=components[node['guide'], node['family']],
                                     pair_key=node['pair_key'], taxon_a=a, taxon_b=b,
                                     kernel_index_a=index[a], kernel_index_b=index[b], selected_record_uses=count))
    # Read exported rows back against raw graph fields and independently accumulated selection uses.
    seen = set()
    for row in csv.DictReader((out / 'selected_node_taxa.tsv').open(), delimiter='\t'):
        role, node_id = row['role'], row['node_id']
        assert (role, node_id) not in seen
        seen.add((role, node_id))
        source = nodes[role][node_id]
        for end in ['a', 'b']:
            expected = source['taxon_id'] if role == 'target' else source['taxon_' + end]
            assert row['taxon_' + end] == expected == taxa[int(row['kernel_index_' + end])]
        for key in ['guide', 'family', 'pair_key']:
            assert row[key] == source[key]
        assert row['family_component'] == components[row['guide'], row['family']]
        assert int(row['selected_record_uses']) == usage[role][node_id]
    assert seen == {(role, node) for role in usage for node in usage[role]}
    for path, digest in pins.items():
        assert sha(path) == digest, path
    result = dict(status='complete_selected_taxon_inputs_with_full_source_readback',
                  selected_records=records, selected_nodes=len(seen),
                  selected_records_with_focal_taxon_in_background=focal_present,
                  unique_nodes_by_role={role: len(counts) for role, counts in usage.items()},
                  pins=pins, script_sha256=sha(__file__),
                  artifacts={'selected_node_taxa.tsv': sha(out / 'selected_node_taxa.tsv')},
                  scope='Lossless node-key join to every original selected record; all policies and scenarios remain in the pinned selections table. Endpoint labels retain source order. Reuse counts are not independent observations. No structural filtering, fitted covariance, effect estimate or assumption about additive taxon effects.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['pins', 'artifacts']}, indent=2))


if __name__ == '__main__':
    main()
