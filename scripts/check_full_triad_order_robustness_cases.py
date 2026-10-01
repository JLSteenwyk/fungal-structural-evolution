#!/usr/bin/env python3
"""Synthetic aggregation/readback software checks; not a biological pilot."""
import argparse
import copy
import csv
import gzip
import json
import shutil
import sqlite3
import tempfile
from pathlib import Path
import summarize_full_triad_order_robustness as producer
import readback_full_triad_order_robustness as reader
from full_triad_robustness_sources import DEFINITIONS, NUMERIC_FIELDS
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); assert not args.output.exists()
    screens = [dict(id='n50_c70', minimum_aligned_residues=50, minimum_original_coverage=.7)]
    triads = [dict(triad_id='synthetic-' + str(i), models=[[role + str(i), 1] for role in ['a', 'b', 'r']]) for i in range(4)]
    with tempfile.TemporaryDirectory(prefix='full-triad-order-software-', dir='results') as folder:
        temp = Path(folder); fit = temp / 'fits.tsv.gz'; records = []
        for ix, triad in enumerate(triads):
            for mask in ['full', 'plddt70']:
                for order in range(8):
                    for definition in DEFINITIONS:
                        row = dict(triad_id=triad['triad_id'], mask=mask, mapping_definition=definition,
                                   order_ab=order // 4, order_ar=(order // 2) % 2, order_br=order % 2,
                                   triples_sha256=str(order), source_exclusions='', fit_status='computed_unique_at_numeric_tolerance')
                        for role, model in zip(['a', 'b', 'reference'], triad['models']): row['model_' + role], row['version_' + role] = model
                        row.update({field: str(order + 1) for field in NUMERIC_FIELDS})
                        row['rmsd_ar_minus_br'] = str(order + 1 if ix == 0 else -order - 1 if ix == 3 else order - 3)
                        row.update(n50_c70_pass='1', n50_c70_core_pass='1', n50_c70_three_pair_pass='1', n50_c70_exclusions='')
                        if ix == 1 and mask == 'plddt70':
                            row['fit_status'] = 'source_excluded'; row['source_exclusions'] = 'failed_native_input'
                            row.update({field: '' for field in NUMERIC_FIELDS}); row['common_residues'] = '0'
                            row.update(n50_c70_pass='0', n50_c70_core_pass='0', n50_c70_exclusions='source_excluded;short_core')
                        elif ix == 1 and order == 0:
                            row['fit_status'] = 'fewer_than_three_common_residues'; row.update({field: '' for field in NUMERIC_FIELDS}); row['common_residues'] = '2'
                            row.update(n50_c70_pass='0', n50_c70_core_pass='0', n50_c70_exclusions='short_core')
                        elif ix == 2 and definition == 'cycle_consistent' and order == 7:
                            row['fit_status'] = 'computed_nonunique_at_numeric_tolerance'
                            row.update(n50_c70_pass='0', n50_c70_core_pass='0', n50_c70_exclusions='nonunique_fit')
                        elif ix == 3 and order % 2:
                            row.update(n50_c70_pass='0', n50_c70_three_pair_pass='0', n50_c70_exclusions='inherited_pair_failure')
                        records.append(row)
        with gzip.open(fit, 'wt') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t'); writer.writeheader(); writer.writerows(records)
        def sources(plan, path):
            return fit, triads, temp / 'unused-contexts', {}, {str(path): sha(path), str(fit): sha(fit)}
        producer.load_sources = reader.load_sources = sources
        plan = dict(output=str(temp / 'baseline'), screens=screens, pins={}, resources=dict(minimum_free_disk_gib=0),
                    expected=dict(correspondence_work_triads=4, fit_rows=128, robustness_groups=16), numeric_readback_absolute_tolerance=1e-9,
                    scope='Synthetic full-grid software fixture only; source/proof I/O stubbed; no biological pilot or source qualification.')
        pp = temp / 'baseline-plan.json'; pp.write_text(json.dumps(plan)); producer.run(pp); baseline = reader.run(pp, temp / 'baseline-readback.json')
        with gzip.open(Path(plan['output']) / 'triad_order_robustness.jsonl.gz', 'rt') as handle: groups = [json.loads(line) for line in handle]
        positive, missing, nonunique, negative = groups[0], groups[6], groups[9], groups[12]
        assert positive['metric_ranges']['rmsd_ab'] == dict(available_orders=8, minimum=1., maximum=8., span=7., mean=4.5)
        assert positive['screens']['n50_c70']['order_pass_bits'] == 255 and positive['strict_all_orders_contrast_direction'] == 'positive'
        assert missing['metric_ranges']['rmsd_ab'] == dict(available_orders=0, minimum=None, maximum=None, span=None, mean=None)
        assert missing['computed_orders'] == 0 and missing['strict_all_orders_contrast_direction'] == 'uncomputed_or_excluded_order'
        assert nonunique['strict_all_orders_contrast_direction'] == 'nonunique_order_fit' and nonunique['screens']['n50_c70']['order_pass_bits'] == 127
        assert negative['strict_all_orders_contrast_direction'] == 'negative' and negative['screens']['n50_c70']['order_pass_bits'] == 85
        mutations = {
            'promoted_all_orders': lambda g, r: g[12]['screens']['n50_c70'].update(all_orders_pass=True, order_pass_bits=255),
            'wrong_range': lambda g, r: g[0]['metric_ranges']['rmsd_ab'].update(minimum=0.),
            'wrong_mean': lambda g, r: g[0]['metric_ranges']['rmsd_ab'].update(mean=4.6),
            'missing_distance_as_zero': lambda g, r: g[6]['metric_ranges']['rmsd_ab'].update(minimum=0., maximum=0., span=0., mean=0.),
            'promoted_nonunique': lambda g, r: g[9].update(all_orders_unique=True),
            'invented_direction': lambda g, r: g[8].update(strict_all_orders_contrast_direction='positive'),
            'lost_exclusion': lambda g, r: g[6]['order_source_exclusions'][0].clear(),
            'changed_order_layout': lambda g, r: g[0]['order_bits_layout'].reverse(),
            'changed_model_role': lambda g, r: g[0]['models'].reverse(),
            'duplicated_group': lambda g, r: g.insert(1, copy.deepcopy(g[0])),
            'missing_group': lambda g, r: g.pop(),
            'invented_scientific_eligibility': lambda g, r: r.update(scientific_eligibility=True),
        }
        base_receipt = json.loads((Path(plan['output']) / 'receipt.json').read_text()); rejected = []
        for name, mutate in mutations.items():
            candidate = copy.deepcopy(groups); receipt = copy.deepcopy(base_receipt); mutate(candidate, receipt)
            case = temp / name; case.mkdir(); cp = temp / (name + '-plan.json'); config = {**plan, 'output': str(case)}; cp.write_text(json.dumps(config))
            with gzip.open(case / 'triad_order_robustness.jsonl.gz', 'wt') as handle:
                for group in candidate: handle.write(json.dumps(group) + '\n')
            receipt.update(plan_sha256=sha(cp), source_hashes=sources(config, cp)[-1], artifacts={'triad_order_robustness.jsonl.gz': sha(case / 'triad_order_robustness.jsonl.gz')})
            (case / 'receipt.json').write_text(json.dumps(receipt))
            try: reader.run(cp, case / 'readback.json')
            except (AssertionError, sqlite3.IntegrityError): rejected.append(name)
            else: raise AssertionError('False export accepted: ' + name)
        assert len(rejected) == len(mutations)
    result = dict(status='passed_synthetic_full_triad_all_order_software_checks', synthetic_triads=4, fit_rows=128, robustness_groups=16,
                  baseline_readback_status=baseline['status'], rejected_rehashed_false_exports=rejected,
                  script_hashes={str(p): sha(p) for p in [Path(__file__), Path(producer.__file__), Path(reader.__file__), Path('scripts/full_triad_robustness_sources.py')]},
                  scope='Synthetic aggregation/null/status/order-bit/sign tests and independent SQL rejection of12 rehashed false exports. Source/proof I/O stubbed only here; does not qualify real closed-source hashes, model geometry or scientific inference; not a pilot.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
