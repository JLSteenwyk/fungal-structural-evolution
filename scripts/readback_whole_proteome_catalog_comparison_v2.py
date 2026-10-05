"""Independently replay both sides of the full catalog-change table using CSV.

Planning: one CPU, up to 4 GiB RAM, under ten minutes; no GPU or new models.
"""
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from ancestral_chain_attempt import write_json
from reference_measurement_union_sources import bind, verify
from readback_whole_proteome_catalog import sha


def rows(path):
    with open(path) as h:
        yield from csv.DictReader(h, delimiter='\t')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();assert not args.receipt.exists()
    plan_path=args.plan
    plan=json.loads(plan_path.read_text())
    adapter_path=Path(plan['comparison_adapter_receipt'])
    transport_path=Path(plan['comparison_original_transport'])
    adapter,transport=[json.loads(path.read_text()) for path in [adapter_path,transport_path]]
    assert adapter['status']=='completed_full_afdb_catalog_comparison_pending_readback'
    assert transport['validation_sha256']==sha(adapter_path)
    assert transport['original_tool_terminal_exit_code']==0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    assert transport['manager_start_records']==transport['manager_completion_records']==1
    pins=dict(adapter['source_hashes'])
    verify(pins)
    for path in [plan_path,adapter_path,transport_path,Path(__file__)]:bind(pins,path)
    plan = json.loads(plan_path.read_text())
    out = Path(plan['output'])
    receipt = json.loads((out / 'receipt.json').read_text())
    assert receipt['status'] == 'complete_validated_catalog_comparison'
    assert receipt['plan_sha256'] == sha(plan_path)
    for name, digest in receipt['artifacts'].items():
        assert sha(out / name) == digest
    for name, digest in receipt['sources'].items():
        assert sha(name) == digest
    fields = ['sequence_sha256', 'model_id', 'version', 'model_path']
    observed = {}
    for side in ('old', 'new'):
        root = Path(plan[side])
        source_receipt = json.loads((root / 'receipt.json').read_text())
        source = root / 'protein_model_links.tsv'
        assert sha(source) == source_receipt['artifacts'][source.name]
        index = {}
        for r in rows(source):
            key = (r['taxon_id'], r['protein_id'])
            assert key not in index
            index[key] = tuple(r[f] for f in fields)
        expected = len(index)
        assert expected == source_receipt['proteins_linked']
        count = 0
        for r in rows(out / 'protein_link_dispositions.tsv'):
            values = tuple(r[f + '_' + side] for f in fields)
            if any(values):
                assert all(values)
                assert index.pop((r['taxon_id'], r['protein_id'])) == values
                count += 1
        assert not index and count == expected
        observed[side] = count
    by_taxon = defaultdict(Counter)
    previous = None
    for r in rows(out / 'protein_link_dispositions.tsv'):
        key = (r['taxon_id'], r['protein_id'])
        assert previous is None or key > previous
        previous = key
        old = tuple(r[f + '_old'] for f in fields)
        new = tuple(r[f + '_new'] for f in fields)
        assert any(old) or any(new)
        if not any(old):
            expected = 'new_catalog_link'
        elif not any(new):
            expected = 'lost_catalog_link'
        else:
            assert old[0] == new[0]
            expected = 'unchanged_model' if old == new else 'changed_selected_model'
        assert expected == r['disposition']
        by_taxon[key[0]][expected] += 1
    totals = sum(by_taxon.values(), Counter())
    assert dict(totals) == receipt['dispositions']
    coverage = {}
    for side in ('old', 'new'):
        path = Path(plan[side]) / 'taxon_coverage.tsv'
        sr = json.loads((path.parent / 'receipt.json').read_text())
        assert sha(path) == sr['artifacts'][path.name]
        coverage[side] = {r['taxon_id']: r for r in rows(path)}
        assert len(coverage[side]) == 526
    seen = set(); unlinked = 0
    for r in rows(out / 'taxon_coverage_change.tsv'):
        taxon = r['taxon_id']; assert taxon not in seen; seen.add(taxon)
        n = int(r['representative_proteins'])
        counts = by_taxon[taxon]
        for disposition in ('new_catalog_link', 'lost_catalog_link', 'changed_selected_model', 'unchanged_model'):
            assert int(r[disposition]) == counts[disposition]
        for side in ('old', 'new'):
            source = coverage[side][taxon]
            for f in ('species_name', 'study_role', 'representative_proteins'):
                assert r[f] == source[f]
            for f in ('proteins_with_model', 'proteins_without_catalog_model'):
                assert r[f + '_' + side] == source[f]
            assert math.isclose(float(r['coverage_fraction_' + side]), int(source['proteins_with_model']) / n, abs_tol=1e-14)
        assert math.isclose(float(r['coverage_fraction_change']), float(r['coverage_fraction_new']) - float(r['coverage_fraction_old']), abs_tol=1e-14)
        assert int(r['unlinked_in_both']) == n - sum(counts.values())
        unlinked += int(r['unlinked_in_both'])
    assert seen == set(coverage['old']) == set(coverage['new'])
    assert unlinked == receipt['unlinked_in_both']
    assert observed['new'] - observed['old'] == receipt['net_link_change']
    result = dict(status='passed_full_catalog_comparison_source_replay',
                  source_receipt_sha256=sha(out / 'receipt.json'), script_sha256=sha(__file__),
                  links_replayed=observed, taxa=len(seen), dispositions=dict(totals),
                  unlinked_in_both=unlinked,
                  scope='Every old and new link field matched to source; sorted unique union, dispositions, all taxon counts and coverage fractions independently checked. No new structural confidence qualification.')
    verify(pins)
    result.update(source_hashes=pins,scientific_eligibility=False,new_predictions=0,gpu=False)
    write_json(args.receipt,result)
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
