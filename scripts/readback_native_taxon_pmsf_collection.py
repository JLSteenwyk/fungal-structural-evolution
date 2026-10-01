#!/usr/bin/env python3
"""Reconstruct every native sensitivity output with DendroPy and Decimal.

This reader does not import the native runner or its Bio.Phylo audit algorithms.
All profiles, raw tree/report edges, 1,000 raw replicates, complete NEXUS split
weights and saved support rows are checked per run. Role-boundary support is a
diagnostic of an unrooted split, not a root assignment or model adequacy test.
"""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import json
import math
from pathlib import Path
import re
import dendropy
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


SUPPORT_FIELDS = ['tree', 'split_taxa_json', 'branch_length', 'sh_alrt_percent',
                  'reported_ufboot_percent', 'empirical_ufboot_percent']


def demand(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(side, universe):
    left, right = tuple(sorted(side)), tuple(sorted(universe.difference(side)))
    return left if (len(left), left) <= (len(right), right) else right


def matrix_taxa(path, columns):
    lengths = {}
    current = None
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current = line[1:].split()[0]
                demand(current not in lengths, 'Duplicate FASTA tip')
                lengths[current] = 0
            else:
                demand(current is not None and not any(c.isspace() for c in line), 'Malformed FASTA')
                lengths[current] += len(line)
    demand(lengths and all(n == columns for n in lengths.values()), 'Matrix column grid differs')
    return set(lengths)


def parsed_edges(tree, universe, binary=False):
    names = [n.taxon.label for n in tree.leaf_node_iter()]
    demand(len(names) == len(universe) and set(names) == universe, 'Tree tip universe differs')
    descendants, edges = {}, {}
    for node in tree.postorder_node_iter():
        descendants[node] = ({node.taxon.label} if node.is_leaf() else
                             set().union(*(descendants[c] for c in node.child_node_iter())))
        if node is tree.seed_node:
            demand(node.edge.length in (None, 0), 'Unexpected root stem')
            continue
        length = node.edge.length
        demand(length is not None and math.isfinite(length) and length >= 0,
               'Invalid native branch length')
        key = canonical(descendants[node], universe)
        demand(key and key not in edges, 'Repeated or empty native split')
        edges[key] = (length, node.label)
    if binary:
        demand(len(edges) == 2 * len(universe) - 3, 'Nonbinary native ML/bootstrap tree')
    return edges


def tree_file(path, universe, binary=False):
    return parsed_edges(dendropy.Tree.get(path=str(path), schema='newick',
                        rooting='force-unrooted', preserve_underscores=True), universe, binary)


def inspect_run(run, matrix, taxa, columns, audit, bindings, audit_status):
    """Independently reconstruct one entire saved run; no output is trusted alone."""
    run, audit, matrix = map(Path, [run, audit, matrix])
    rp, cp, ap = run/'receipt.json', run/'config.json', audit/'receipt.json'
    r, c, a = [json.loads(p.read_text()) for p in [rp, cp, ap]]
    demand(r['status'] == 'complete_pmsf_execution_pending_full_audit' and r['returncode'] == 0,
           'Successful native execution receipt required')
    demand(r['taxa'] == taxa and r['columns'] == columns and r['config_sha256'] == sha(cp),
           'Native execution identity or dimensions differ')
    demand(a['status'] == audit_status and a['source_receipt_sha256'] == sha(rp),
           'First audit lineage differs')
    for path in [rp, cp, ap, matrix]:
        bind(bindings, path)
    for root, record in [(run, r), (audit, a)]:
        for name, digest in record['artifacts'].items():
            bind(bindings, root/name, digest)
    for path, digest in c['pinned_files'].items():
        bind(bindings, path, digest)
    verify(bindings)
    command = c['command']
    for flag, value in [('-m', 'LG+C20+F+G4'), ('--alrt', '1000'), ('-B', '1000')]:
        demand(command.count(flag) == 1 and command[command.index(flag)+1] == value,
               'Native inference settings differ')
    demand(all(flag in command for flag in ['--tree-freq', '--bnni', '--boot-trees']),
           'Native support/conditioning settings absent')
    demand(Path(command[command.index('-s')+1]).resolve() == matrix.resolve(), 'Wrong native matrix')
    demand(Path(command[command.index('--tree-freq')+1]).resolve() == (run/'input_guide.treefile').resolve(),
           'Wrong frozen guide')
    universe = matrix_taxa(matrix, columns)
    demand(len(universe) == taxa, 'Matrix tip count differs')
    profiles = 0
    max_error = Decimal(0)
    with (run/'pmsf.sitefreq').open() as handle:
        for profiles, line in enumerate(handle, 1):
            fields = line.split()
            demand(len(fields) == 21 and int(fields[0]) == profiles, 'Profile site order/grid differs')
            values = [Decimal(v) for v in fields[1:]]
            demand(all(v.is_finite() and 0 <= v <= 1 for v in values), 'Invalid profile probability')
            error = abs(sum(values) - 1)
            demand(error <= Decimal('0.00001'), 'Profile normalization differs at printed precision')
            max_error = max(error, max_error)
    demand(profiles == columns, 'Profile site count differs')
    views = {kind: tree_file(run/name, universe, kind != 'consensus') for kind, name in
             [('ml', 'pmsf.treefile'), ('consensus', 'pmsf.contree'), ('guide', 'input_guide.treefile')]}
    report, log = [(run/name).read_text() for name in ['pmsf.iqtree', 'pmsf.log']]
    demand(f'Input data: {taxa} sequences with {columns} amino-acid sites' in report and
           'Model of substitution: LG+SSF+F+G4' in report and
           'SH-aLRT support (%) / ultrafast bootstrap support (%)' in report and
           'Computing posterior mean site frequencies' in log, 'Report/model/profile declaration differs')
    reported = [line.strip() for line in report.splitlines() if line.startswith('(') and line.rstrip().endswith(';')]
    demand(len(reported) == 2, 'Report tree count differs')
    for kind, text in zip(['ml', 'consensus'], reported):
        parsed = dendropy.Tree.get(data=text, schema='newick', rooting='force-unrooted', preserve_underscores=True)
        demand(parsed_edges(parsed, universe, kind == 'ml') == views[kind], 'Raw/report tree differs')
    counts, replicates = Counter(), 0
    namespace = dendropy.TaxonNamespace(sorted(universe))
    for parsed in dendropy.Tree.yield_from_files([str(run/'pmsf.ufboot')], schema='newick',
            taxon_namespace=namespace, rooting='force-unrooted', preserve_underscores=True):
        counts.update(parsed_edges(parsed, universe, True).keys())
        replicates += 1
    demand(replicates == 1000, 'Native bootstrap count differs')
    nexus = (run/'pmsf.splits.nex').read_text()
    label_rows = re.findall(r"^\[(\d+)\] '([^']+)'$", nexus, re.M)
    labels = {int(i): name for i, name in label_rows}
    demand(len(label_rows) == taxa and set(labels) == set(range(1, taxa+1)) and
           len(set(labels.values())) == taxa and set(labels.values()) == universe, 'NEXUS tip labels differ')
    match = re.search(r'\bMATRIX\s*\n(.*?)\n\s*;', nexus, re.S)
    demand(match is not None, 'NEXUS split matrix absent')
    seen = set()
    for line in match[1].splitlines():
        fields = line.strip().rstrip(',').split()
        if not fields:
            continue
        weight, members = Decimal(fields[0]), list(map(int, fields[1:]))
        demand(weight.is_finite() and 0 <= weight <= 100 and members and
               len(members) == len(set(members)) and set(members) <= set(labels), 'Invalid NEXUS split')
        key = canonical({labels[i] for i in members}, universe)
        demand(key not in seen and key in counts, 'Repeated or unsupported NEXUS split')
        seen.add(key)
        tolerance = Decimal('0.5') * Decimal(10) ** weight.as_tuple().exponent
        demand(abs(weight - Decimal(100)*counts[key]/replicates) <= tolerance + Decimal('0.00000001'),
               'NEXUS weight differs from all raw bootstrap trees')
    demand(seen == set(counts), 'NEXUS/empirical split universes differ')
    rows = []
    expected_rows = {}
    for kind in ['ml', 'consensus']:
        for key, (length, label) in views[kind].items():
            if len(key) < 2:
                continue
            if kind == 'ml':
                alrt, boot = map(float, (label or '').split('/'))
                demand(math.isfinite(alrt) and 0 <= alrt <= 100, 'Invalid SH-aLRT label')
            else:
                alrt, boot = None, float(label)
            empirical = 100 * counts[key]/replicates
            demand(math.isfinite(boot) and 0 <= boot <= 100 and abs(boot-empirical) <= .500001,
                   'Native tree support differs from empirical frequency')
            row = dict(tree=kind, split_taxa_json=json.dumps(key), branch_length=length,
                       sh_alrt_percent=alrt, reported_ufboot_percent=boot, empirical_ufboot_percent=empirical)
            expected_rows[kind, key] = row
            rows.append(row)
    with (audit/'branch_support.tsv').open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        demand(reader.fieldnames == SUPPORT_FIELDS, 'First-audit support schema differs')
        actual_seen = set()
        for row in reader:
            side = json.loads(row['split_taxa_json'])
            demand(len(side) == len(set(side)) and set(side) <= universe, 'Invalid first-audit split')
            key = row['tree'], canonical(set(side), universe)
            demand(key in expected_rows and key not in actual_seen, 'Repeated or unexpected first-audit support')
            actual_seen.add(key)
            expected = expected_rows[key]
            for field in SUPPORT_FIELDS[2:]:
                demand((row[field] == '' if expected[field] is None else float(row[field]) == expected[field]),
                       'First audit field differs: '+field)
        demand(actual_seen == set(expected_rows), 'Missing first-audit support rows')
    rf = len(set(views['ml']) ^ set(views['consensus']))
    demand(rf == int(re.search(r'Robinson-Foulds distance between ML tree and consensus tree: (\d+)', report)[1]),
           'Report RF differs')
    likelihoods = [float(re.search(pattern, report)[1]) for pattern in
                   [r'Log-likelihood of the tree: ([-\d.]+)', r'Log-likelihood of consensus tree: ([-\d.]+)']]
    demand(all(math.isfinite(v) for v in likelihoods), 'Nonfinite native likelihood')
    total = float(re.search(r'Total tree length \(sum of branch lengths\): ([\d.]+)', report)[1])
    demand(abs(math.fsum(edge[0] for edge in views['ml'].values())-total) < .000051, 'Total branch length differs')
    for field, value in dict(taxa=taxa, sites=columns, bootstrap_trees=1000, empirical_splits=len(counts),
            internal_support_rows=len(rows), ml_consensus_rf=rf, reported_ml_log_likelihood=likelihoods[0],
            reported_consensus_log_likelihood=likelihoods[1]).items():
        demand(a[field] == value, 'First-audit summary differs: '+field)
    verify(bindings)
    summary = dict(taxa=taxa, sites=columns, bootstrap_trees=1000, empirical_splits=len(counts),
                   internal_support_rows=len(rows), ml_consensus_rf=rf,
                   max_decimal_profile_sum_error=str(max_error), source_receipt_sha256=sha(rp),
                   first_audit_receipt_sha256=sha(ap))
    return summary, rows, universe


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = dict(plan['pins'])
    bind(bindings, args.plan)
    np = Path(plan['native_plan'])
    native = json.loads(np.read_text())
    nrp = Path(native['output'])/'receipt.json'
    nr = json.loads(nrp.read_text())
    demand(nr['status'] == 'complete_native_taxon_pmsf_batch_pending_full_independent_collection_readback' and
           nr['plan_sha256'] == sha(np) and len(nr['runs']) == 16, 'Complete native batch required')
    closed = json.loads(Path(native['input_completion']).read_text())
    demand(closed['status'] == 'complete_verified_native_species_taxon_refit_inputs' and
           len(closed['services']) == 2 and closed['summary']['matrices'] == 8, 'Closed full inputs required')
    for r in [nr, closed]:
        for path, digest in r['source_hashes'].items():
            bind(bindings, path, digest)
    for p in [np, nrp, native['input_completion']]:
        bind(bindings, p)
    ip = Path(native['inputs'])/'receipt.json'
    bind(bindings, ip)
    inputs = json.loads(ip.read_text())
    matrices = {(r['alignment'], r['policy']): r for r in inputs['matrices']}
    demand(len(matrices) == 8 and len(native['jobs']) == 16 and
           len({(j['policy'], j['alignment'], j['guide_alignment']) for j in native['jobs']}) == 16,
           'Full sensitivity job grid differs')
    verify(bindings)
    output = Path(plan['output'])
    demand(not output.exists(), 'Use a new immutable collection directory')
    output.mkdir(parents=True)
    summaries, boundary = [], []
    for job, record in zip(native['jobs'], nr['runs']):
        demand(all(record[k] == v for k, v in job.items()), 'Native job order/identity differs')
        run = Path(native['output'])/job['label']
        demand(Path(record['run']) == run and Path(record['audit']) == run/'audit' and
               record['run_receipt_sha256'] == sha(run/'receipt.json') and
               record['audit_receipt_sha256'] == sha(run/'audit'/'receipt.json'), 'Batch run lineage differs')
        config = json.loads((run/'config.json').read_text())
        demand(config['batch_plan_sha256'] == sha(np) and
               all(config[k] == job[k] for k in ['policy', 'alignment', 'guide_alignment']), 'Run config identity differs')
        first_audit = json.loads((run/'audit'/'receipt.json').read_text())
        demand(first_audit['script_sha256'] == native['pins']['scripts/audit_species_taxon_pmsf.py'],
               'First audit implementation differs from pinned native batch')
        spec = matrices[job['alignment'], job['policy']]
        demand(record['taxa'] == spec['taxa'] and record['columns'] == spec['columns'], 'Batch matrix dimensions differ')
        matrix = Path(spec['path'])/'matrix.faa'
        demand(sha(matrix) == spec['matrix_sha256'], 'Prepared matrix changed')
        summary, rows, universe = inspect_run(run, matrix, spec['taxa'], spec['columns'], run/'audit', bindings,
            'passed_taxon_sensitivity_pmsf_profile_tree_and_bootstrap_readback')
        guide = Path(native['output'])/'guides'/(job['guide_alignment']+'-'+job['policy'])
        gp, gr = guide/'receipt.json', json.loads((guide/'receipt.json').read_text())
        guide_spec = matrices[job['guide_alignment'], job['policy']]
        demand(gr['status'] == 'checked_full_taxon_sensitivity_conditioning_guide' and
               gr['plan_sha256'] == sha(np) and gr['matrix_sha256'] == guide_spec['matrix_sha256'] and
               gr['taxa'] == spec['taxa'] and (guide/'guide.treefile').read_bytes() == (run/'input_guide.treefile').read_bytes(),
               'Conditioning guide lineage differs')
        bind(bindings, gp)
        for name, digest in gr['artifacts'].items():
            bind(bindings, guide/name, digest)
        tp = Path(spec['path'])/'taxa.tsv'
        bind(bindings, tp)
        with tp.open() as handle:
            taxa_rows = list(csv.DictReader(handle, delimiter='\t'))
        demand(len(taxa_rows) == len(universe) and {r['taxon_id'] for r in taxa_rows} == universe,
               'Prepared role tip grid differs')
        roles = Counter(r['study_role'] for r in taxa_rows)
        demand(dict(roles) == spec['roles'], 'Sensitivity role counts differ')
        key = canonical({r['taxon_id'] for r in taxa_rows if r['study_role'] == 'outgroup'}, universe)
        for kind in ['ml', 'consensus']:
            found = [row for row in rows if row['tree'] == kind and tuple(json.loads(row['split_taxa_json'])) == key]
            demand(len(found) <= 1, 'Repeated role boundary')
            b = found[0] if found else None
            boundary.append(dict(run=job['label'], tree=kind, ingroup=roles['ingroup'], outgroup=roles['outgroup'],
                boundary_taxa_json=json.dumps(key), present=b is not None,
                sh_alrt_percent=b['sh_alrt_percent'] if b else None,
                empirical_ufboot_percent=b['empirical_ufboot_percent'] if b else None))
        path = output/(job['label']+'.support.tsv')
        with path.open('x') as handle:
            writer = csv.DictWriter(handle, SUPPORT_FIELDS, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda r: (r['tree'], r['split_taxa_json'])))
        summaries.append(dict(**job, **summary))
        print('independently_read_back_full_native_run', job['label'], len(summaries), '/16', flush=True)
    with (output/'role_boundary.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, list(boundary[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(boundary)
    verify(bindings)
    proof = dict(status='passed_full_native_taxon_pmsf_independent_collection_pending_journal_closure',
        plan_sha256=sha(args.plan), native_batch_receipt_sha256=sha(nrp), runs=summaries, run_count=16,
        tree_views=32, raw_bootstrap_trees=16000, site_profiles=sum(r['sites'] for r in summaries),
        internal_support_rows=sum(r['internal_support_rows'] for r in summaries),
        role_boundary_rows=len(boundary), boundary_views_with_role_split=sum(r['present'] for r in boundary),
        dendropy_version=dendropy.__version__, script_sha256=sha(Path(__file__)), source_hashes=bindings,
        artifacts={p.name:sha(p) for p in output.iterdir() if p.is_file()}, scientific_eligibility=False,
        scope=plan['scope'])
    with (output/'receipt.json').open('x') as handle:
        handle.write(json.dumps(proof, indent=2)+'\n')
    print(json.dumps({k:proof[k] for k in ['status','run_count','raw_bootstrap_trees','role_boundary_rows']}), flush=True)


if __name__ == '__main__':
    main()
