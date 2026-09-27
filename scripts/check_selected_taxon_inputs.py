#!/usr/bin/env python3
"""Independent dataframe reconstruction of full selected-node taxon inputs."""
import json
from pathlib import Path
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha


def main():
    out = Path('results/orthology/selected-taxon-inputs-20260927-v1')
    receipt = json.loads((out / 'receipt.json').read_text())
    for path, digest in receipt['pins'].items():
        assert sha(path) == digest, path
    for name, digest in receipt['artifacts'].items():
        assert sha(out / name) == digest
    actual = pd.read_csv(out / 'selected_node_taxa.tsv', sep='\t', dtype=str)
    selection = Path('results/orthology/background-control-selection-20260927-v1')
    selected = pd.read_csv(selection / 'selections.tsv.gz', sep='\t', usecols=['target_id', 'background_id'])
    components = pd.read_csv('results/orthology/selected-control-family-dependencies-20260927-v1/family_components.tsv', sep='\t', dtype=str)
    components = components[components.entity_kind.eq('combined')][['guide', 'family', 'component_id']]
    taxa = json.loads(Path('results/phylogeny/species-distance-kernels-20260927-v1/taxa.json').read_text())
    expected_parts = []
    for role in ['target', 'background']:
        source = pd.read_json('results/orthology/background-match-graph-20260927-v1/' + role + '_nodes.jsonl', lines=True, dtype=False)
        counts = selected[role + '_id'].value_counts().rename_axis('node_id').rename('selected_record_uses').reset_index()
        source = counts.merge(source, on='node_id', validate='one_to_one')
        assert len(source) == len(counts)
        source = source.merge(components, on=['guide', 'family'], validate='many_to_one')
        assert len(source) == len(counts)
        source['family_component'] = source.component_id
        source['role'] = role
        for end in ['a', 'b']:
            if role == 'target':
                source['taxon_' + end] = source.taxon_id
            source['kernel_index_' + end] = source['taxon_' + end].map({t: i for i, t in enumerate(taxa)})
            assert source['kernel_index_' + end].notna().all()
        expected_parts.append(source[list(actual.columns)].astype(str))
    expected = pd.concat(expected_parts, ignore_index=True)
    keys = ['role', 'node_id']
    pd.testing.assert_frame_equal(actual.sort_values(keys).reset_index(drop=True), expected.sort_values(keys).reset_index(drop=True))
    target = actual[actual.role.eq('target')].set_index('node_id')
    background = actual[actual.role.eq('background')].set_index('node_id')
    focal = selected.target_id.map(target.taxon_a)
    present = focal.eq(selected.background_id.map(background.taxon_a)) | focal.eq(selected.background_id.map(background.taxon_b))
    for key in ['guide', 'family']:
        assert selected.target_id.map(target[key]).eq(selected.background_id.map(background[key])).all()
    assert int(present.sum()) == receipt['selected_records_with_focal_taxon_in_background']
    assert len(selected) == receipt['selected_records']
    result = dict(status='passed_full_selected_taxon_input_readback', selected_records=len(selected),
                  selected_nodes=len(actual), exported_cells_checked=int(actual.size),
                  selected_records_with_focal_taxon_in_background=int(present.sum()),
                  producer_receipt_sha256=sha(out / 'receipt.json'), script_sha256=sha(__file__),
                  scope='Independent dataframe joins and grouped reuse counts reconstruct every exported cell; all selected records checked for family/guide consistency and focal/background taxon overlap. No model fitted.')
    Path('metadata/selected_taxon_inputs_readback_20260927.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
