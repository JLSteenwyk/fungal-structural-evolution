#!/usr/bin/env python3
"""Exercise native gCF semantics and reject altered exports on synthetic trees."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from readback_species_gene_concordance_v2 import validate_run
from run_ortholog_pair_guide_comparison import sha


REFERENCE = '((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2,(E:0.1,F:0.1):0.2);'
GENES = [REFERENCE,
         '((A:0.1,(C:0.1,D:0.1):0.2):0.2,B:0.1,(E:0.1,F:0.1):0.2);',
         '((A:0.1,(E:0.1,F:0.1):0.2):0.2,B:0.1,(C:0.1,D:0.1):0.2);',
         '((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2,(E:0.1,F:0.1):0.2);',
         '(A:0.1,B:0.1,(C:0.1,D:0.1):0.2);']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    records = []
    for name, genes in [('all_states', GENES), ('zero_decisive', [GENES[-1]]),
                        ('fractional_native_precision', GENES + [GENES[1], GENES[3]])]:
        folder = args.output / name
        folder.mkdir()
        reference = folder / 'reference.tree'
        trees = folder / 'genes.tree'
        reference_text = REFERENCE.replace(':0.2', ':0.1234567890') if name == 'fractional_native_precision' else REFERENCE
        reference.write_text(reference_text + '\n')
        trees.write_text('\n'.join(genes) + '\n')
        command = [args.executable, '-t', str(reference), '--gcf', str(trees), '--cf-verbose',
                   '-T', '2', '--mem', '1G', '--seed', '20261002', '--prefix', str(folder / 'gcf')]
        with (folder / 'stdout.log').open('x') as handle:
            subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
        summaries, cells = validate_run(folder, reference_text, genes, set('ABCDEF'))
        assert len(summaries) == 3 and len(cells) == 3 * len(genes)
        if name == 'all_states':
            by_split = {tuple(row['split_taxa']): row for row in summaries}
            # Integer mask canonicalization chooses AB, CD, and ABCD for the EF edge.
            assert by_split['A', 'B']['decisive'] == 4
            assert by_split['A', 'B']['concordant'] == 1
            assert by_split['A', 'B']['alternative_1'] == by_split['A', 'B']['alternative_2'] == 1
            assert by_split['A', 'B']['residual_discordance'] == 1
            assert by_split['C', 'D']['concordant'] == 3
            assert by_split['A', 'B', 'C', 'D']['concordant'] == 4
            assert sum(row['state'] == 'not_decisive_missing_incident_clade' for row in cells) == 3
        elif name == 'zero_decisive':
            assert all(row['decisive'] == 0 and row['native_statistics']['gCF'] == 'NA' for row in summaries)
        else:
            by_split = {tuple(row['split_taxa']): row for row in summaries}
            assert by_split['A', 'B']['native_statistics']['gCF'] == '16.67'
            assert by_split['A', 'B']['native_statistics']['Length'] == '0.123457'
        records.append(dict(case=name, branches=len(summaries), cells=len(cells), command=command,
                            artifacts={p.name: sha(p) for p in folder.iterdir() if p.is_file()}))
    good = args.output / 'all_states'
    bad_cases = [
        ('cell_concordance', 'gcf.cf.stat_tree', '7\t1\t1\t0\t0', '7\t1\t0\t0\t0'),
        ('missing_clade_as_zero', 'gcf.cf.stat_tree', '7\t5\tNA\tNA\tNA', '7\t5\t0\t0\t0'),
        ('changing_nni_orientation', 'gcf.cf.stat_tree', '7\t2\t0\t1\t0', '7\t2\t0\t0\t1'),
        ('invalid_gene_id', 'gcf.cf.stat_tree', '7\t5\tNA\tNA\tNA', '7\t6\tNA\tNA\tNA'),
        ('duplicate_cell', 'gcf.cf.stat_tree', '7\t5\tNA\tNA\tNA', '7\t1\t1\t0\t0'),
        ('missing_cell', 'gcf.cf.stat_tree', '7\t5\tNA\tNA\tNA\n', ''),
        ('aggregate_concordance', 'gcf.cf.stat', '7\t25\t1', '7\t50\t2'),
        ('aggregate_decisive', 'gcf.cf.stat', '\t1\t4\t\t0.2', '\t1\t5\t\t0.2'),
        ('wrong_percentage', 'gcf.cf.stat', '7\t25\t1', '7\t26\t1'),
        ('incorrect_label', 'gcf.cf.stat', '\t4\t\t0.2', '\t4\twrong\t0.2'),
        ('aggregate_length', 'gcf.cf.stat', '\t0.2\n', '\t0.3\n'),
        ('branch_id_collision', 'gcf.cf.branch', ')8:', ')7:'),
        ('branch_length', 'gcf.cf.branch', '7:0.2000000000', '7:0.3000000000'),
        ('annotated_factor', 'gcf.cf.tree', ')25:', ')26:'),
        ('annotated_taxon', 'gcf.cf.tree', 'A:', 'Z:'),
        ('nexus_cell', 'gcf.cf.tree.nex', 'gC1="1"', 'gC1="0"'),
        ('nexus_factor', 'gcf.cf.tree.nex', 'gCF="25"', 'gCF="26"'),
        ('nexus_count', 'gcf.cf.tree.nex', 'gN="4"', 'gN="5"'),
        ('nexus_combined_factors', 'gcf.cf.tree.nex', '25/25/25/25', '26/25/25/25'),
        ('nexus_combined_counts', 'gcf.cf.tree.nex', '1/1/1/1', '2/1/1/1'),
        ('nexus_length', 'gcf.cf.tree.nex', ':0.2000000000', ':0.3000000000'),
        ('nexus_taxon', 'gcf.cf.tree.nex', 'A:', 'Z:'),
    ]
    rejected = []
    for name, file, old, new in bad_cases:
        with tempfile.TemporaryDirectory(prefix='fungal-gcf-case-') as temporary:
            folder = Path(temporary) / 'native'
            shutil.copytree(good, folder)
            path = folder / file
            text = path.read_text()
            assert old in text, (name, text)
            path.write_text(text.replace(old, new, 1))
            # Record the modified hash to demonstrate that semantic checking is
            # needed even when a receipt has been rehashed after an alteration.
            digest = sha(path)
            try:
                validate_run(folder, REFERENCE, GENES, set('ABCDEF'))
            except (AssertionError, ValueError, KeyError):
                rejected.append(dict(case=name, modified_sha256=digest))
            else:
                raise AssertionError('Accepted false export: ' + name)
    for name, genes in [('gene_order_changed', [GENES[1], GENES[0]] + GENES[2:]),
                        ('taxon_coverage_changed', GENES[:-1] + [REFERENCE])]:
        try:
            validate_run(good, REFERENCE, genes, set('ABCDEF'))
        except AssertionError:
            rejected.append(dict(case=name))
        else:
            raise AssertionError('Accepted wrong gene source: ' + name)
    for name, old, new in [('zero_factor_instead_of_na', '7\tNA\tNA', '7\t0\tNA'),
                           ('zero_count_instead_of_na', '7\tNA\tNA', '7\tNA\t0')]:
        with tempfile.TemporaryDirectory(prefix='fungal-gcf-zero-') as temporary:
            folder = Path(temporary) / 'native'
            shutil.copytree(args.output / 'zero_decisive', folder)
            path = folder / 'gcf.cf.stat'
            text = path.read_text()
            assert old in text
            path.write_text(text.replace(old, new, 1))
            try:
                validate_run(folder, REFERENCE, [GENES[-1]], set('ABCDEF'))
            except AssertionError:
                rejected.append(dict(case=name, modified_sha256=sha(path)))
            else:
                raise AssertionError('Accepted zero substituted for native NA: ' + name)
    result = dict(status='passed_native_synthetic_gcf_all_states_zero_decisive_and_false_export_checks',
                  native_cases=records, false_exports_rejected=len(rejected), rejected=rejected,
                  executable=args.executable, executable_sha256=sha(args.executable),
                  script_sha256=sha(__file__), reader_sha256=sha('scripts/readback_species_gene_concordance_v2.py'),
                  scientific_eligibility=False,
                  scope='Synthetic six-taxon software contracts, not biological inference, a pilot, or production completion. All native cells checked; missing decisiveness stays NA and NNI orientation must be consistent.')
    with (args.output / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], native_cases=len(records), false_exports_rejected=len(rejected))), flush=True)


if __name__ == '__main__':
    main()
