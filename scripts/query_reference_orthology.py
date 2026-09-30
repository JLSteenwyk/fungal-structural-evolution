#!/usr/bin/env python3
"""Query native orthology for all nearest-reference ties, preserving every context."""
import argparse
import csv
import hashlib
import json
import sqlite3
import subprocess
from collections import Counter
from pathlib import Path

from run_ortholog_pair_guide_comparison import sha


def json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def verify(bindings):
    checked = {}
    for path, digest in bindings.items():
        resolved = str(Path(path).resolve())
        if resolved not in checked:
            checked[resolved] = sha(path)
        if checked[resolved] != digest:
            raise ValueError('Changed source: ' + path)


def load_sources(plan):
    """Only source I/O and proof lineage; shared with the independent reader."""
    root = Path(plan['source'])
    receipt = json.loads((root / 'receipt.json').read_text())
    proof = json.loads(Path(plan['source_readback']).read_text())
    closure = json.loads(Path(plan['source_completion']).read_text())
    assert receipt['status'] == 'complete_sequence_first_sister_reference_inventory_pending_readback'
    assert proof['status'] == 'passed_full_sequence_first_sister_reference_native_readback'
    assert closure['status'] == 'complete_verified_expanded_sequence_first_references'
    assert receipt['plan_sha256'] == sha(plan['source_plan'])
    assert proof['producer_receipt_sha256'] == sha(root / 'receipt.json')
    source_plan = json.loads(Path(plan['source_plan']).read_text())
    assert Path(source_plan['output']).resolve() == root.resolve()
    bindings = dict(plan['pins'])
    for lineage in [receipt['source_hashes'], proof['source_hashes'], closure['source_hashes']]:
        for path, digest in lineage.items():
            if path in bindings:
                assert bindings[path] == digest
            bindings[path] = digest
    for name, digest in receipt['artifacts'].items():
        bindings[str(root / name)] = digest
    assert closure['source_hashes'][str(root / 'receipt.json')] == proof['producer_receipt_sha256']
    assert closure['source_hashes'][plan['source_readback']] == sha(plan['source_readback'])
    summaries = {r['guide']: r for r in receipt['guides']}
    assert proof['guides'] == closure['summary']['guides']
    reference_plan = json.loads(Path(source_plan['reference_plan']).read_text())
    review_plan = json.loads(Path(reference_plan['review_plan']).read_text())
    assert Path(reference_plan['bridge']).resolve() == Path(plan['bridge']).resolve()
    guide_sources = {r['guide']: r for r in review_plan['guides']}
    assert set(guide_sources) == {'profile', 'mafft'}
    maps = []
    for stream in plan['streams']:
        guide = stream['guide']
        audit = json.loads(Path(stream['receipt']).read_text())
        audit_plan = json.loads(Path(stream['audit_plan']).read_text())
        assert audit['status'] == 'complete_native_ortholog_pair_multiplicity_audit'
        assert audit['plan_sha256'] == sha(stream['audit_plan'])
        assert audit['unique_unordered_pairs'] * 28 == Path(stream['stream']).stat().st_size
        for flag in ['duplicate_directed_incidences', 'pairs_missing_reverse',
                     'pairs_with_repeated_direction', 'pairs_with_unequal_multiplicity']:
            assert audit[flag] == 0
        identity_path = audit_plan['identity_plan']
        assert audit_plan['pins'][identity_path] == sha(identity_path)
        identity_plan = json.loads(Path(identity_path).read_text())
        identity_receipt = Path(identity_plan['output']) / 'receipt.json'
        assert sha(identity_receipt) == audit['identity_receipt_sha256']
        identities = json.loads(identity_receipt.read_text())
        assert identities['status'] == 'complete_grouped_ortholog_protein_family_identity_audit'
        assert identities['plan_sha256'] == sha(identity_path)
        snapshot_path = identity_plan['supplement_receipt']
        assert identity_plan['pins'][snapshot_path] == sha(snapshot_path)
        snapshot = json.loads(Path(snapshot_path).read_text())
        assert snapshot['status'] == 'complete_separate_small_family_ortholog_supplement'
        tree_source = guide_sources[guide]
        tree_proof = json.loads(Path(tree_source['tree_readback']).read_text())
        assert tree_proof['status'] == 'passed_complete_resolved_tree_membership_readback'
        assert tree_proof['resolved_tree_file_sha256'] == sha(tree_source['trees'])
        mapping = []
        for name in ['species_ids', 'sequence_ids']:
            path = stream[name]
            digest = sha(path)
            assert snapshot['input_hashes'][path] == digest
            assert tree_proof['input_hashes'][path] == digest
            bindings[path] = digest
            mapping.append(digest)
        maps.append(mapping)
        for path in [stream['receipt'], stream['audit_plan'], identity_path,
                     str(identity_receipt), snapshot_path, tree_source['tree_readback']]:
            bindings[path] = sha(path)
        bindings[stream['stream']] = audit['sorted_pairs_sha256']
    assert len(maps) == 2 and maps[0] == maps[1], 'Native ordinal universes differ'
    assert {s['guide'] for s in plan['streams']} == {'profile', 'mafft'}
    assert sum(r['candidates'] for r in summaries.values()) == plan['resources']['target_contexts']
    for path, digest in plan['pins'].items():
        assert bindings[path] == digest, 'Overwritten plan binding: ' + path
    verify(bindings)
    return receipt, bindings


def source_rows(plan):
    root = Path(plan['source'])
    for guide in ['profile', 'mafft']:
        with (root / (guide + '_sequence_first_references.tsv')).open() as handle:
            for number, row in enumerate(csv.DictReader(handle, delimiter='\t'), 1):
                yield guide, number, row


def reference_sets(row):
    available = json.loads(row['nearest_reference_genes'])
    sequence = json.loads(row['sequence_first_nearest_genes'])
    for genes in [available, sequence]:
        assert genes == sorted(set(genes))
        assert all(g.split('_', 1)[0] != row['taxon_id'] for g in genes)
    assert int(row['sequence_first_ties']) == len(sequence)
    assert row['sequence_first_chosen_gene'] == (sequence[0] if sequence else '')
    assert row['chosen_reference_gene'] == (available[0] if available else '')
    return [('availability', available), ('sequence_first', sequence)]


def ordinals(plan, needed):
    mapping = plan['streams'][0]
    species = {}
    for line in Path(mapping['species_ids']).read_text().splitlines():
        native, label = line.split(': ', 1)
        assert native not in species
        species[native] = label.rsplit('.', 1)[0]
    assert len(set(species.values())) == len(species)
    ids = {}
    with Path(mapping['sequence_ids']).open() as handle:
        for index, line in enumerate(handle):
            native, protein = line.rstrip('\n').split(': ', 1)
            gene = species[native.split('_')[0]] + '_' + protein
            if gene in needed:
                assert gene not in ids and index < 2**24
                ids[gene] = index
    assert set(ids) == needed
    return ids


def pair_key(ids, left, right):
    a, b = sorted([ids[left], ids[right]])
    assert a < b
    return f'{a:06x}{b:06x}'


def coorthology(a, b):
    return ['neither', 'only_b', 'only_a', 'both'][2 * int(a) + int(b)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    receipt, bindings = load_sources(plan)
    bindings[str(args.plan)] = sha(args.plan)
    rows = list(source_rows(plan))
    assert len(rows) == plan['resources']['target_contexts']
    counts = Counter(g for g, _, _ in rows)
    assert counts == {r['guide']: r['candidates'] for r in receipt['guides']}
    identities = [(g, r['family'], r['gene_node'], r['gene_a'], r['gene_b']) for g, _, r in rows]
    assert len(identities) == len(set(identities))
    needed = {r[k] for _, _, r in rows for k in ['gene_a', 'gene_b']}
    references = {gene for _, _, row in rows for _, genes in reference_sets(row) for gene in genes}
    needed.update(references)
    ids = ordinals(plan, needed)
    models = {}
    with sqlite3.connect('file:' + str(Path(plan['bridge']).resolve()) + '?mode=ro', uri=True) as conn:
        for gene, mid, version in conn.execute("SELECT taxon_id || '_' || protein_id,model_id,version FROM structures"):
            if gene in references:
                assert gene not in models
                models[gene] = (mid, str(version))
    queries = {pair_key(ids, row[side], gene) for _, _, row in rows
               for _, genes in reference_sets(row) for gene in genes for side in ['gene_a', 'gene_b']}
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    query_path = out / 'queries.hex'
    query_path.write_text(''.join(k + '\n' for k in sorted(queries)))
    binary = out / 'query_membership'
    compile_command = [plan['compiler'], '-O2', '-std=c++17', plan['cpp'], '-o', str(binary)]
    subprocess.run(compile_command, check=True)
    subprocess.run([plan['python'], plan['fixtures'], '--binary', str(binary)], check=True)
    membership, stats = {}, {}
    for stream in plan['streams']:
        guide = stream['guide']
        native = json.loads(Path(stream['receipt']).read_text())
        result = out / (guide + '_membership.tsv')
        command = [str(binary), stream['stream'], str(query_path), str(result)]
        run = subprocess.run(command, capture_output=True, text=True, check=True)
        totals = json.loads(run.stdout)
        assert totals['stream_pairs'] == native['unique_unordered_pairs'] and totals['queries'] == len(queries)
        found = {}
        with result.open() as handle:
            for line in handle:
                key, present = line.rstrip('\n').split('\t')
                assert key not in found and present in ['0', '1']
                found[key] = present
        assert set(found) == queries and sum(map(int, found.values())) == totals['present']
        membership[guide], stats[guide] = found, totals
        print(guide, json.dumps(totals), flush=True)
    fields = ['source_guide', 'source_row_number', 'source_row_sha256', 'family', 'gene_node',
              'taxon_id', 'source_context_status', 'sequence_first_status', 'design',
              'reference_gene', 'reference_model', 'reference_version', 'lexical_choice',
              'parent_context_eligible', 'duplicate_side', 'duplicate_gene', 'native_pair_key',
              'profile_native_ortholog', 'mafft_native_ortholog']
    summary = Counter()
    link_count = tie_count = 0
    with (out / 'reference_membership.tsv').open('x') as table, (out / 'contexts.jsonl').open('x') as contexts:
        writer = csv.DictWriter(table, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for guide, number, row in rows:
            rh = json_hash(row)
            eligible = row['status'] in ['provisional_reference_available', 'no_modeled_nonfocal_sister']
            context = dict(source_guide=guide, source_row_number=number, source_row_sha256=rh,
                           source=row, parent_context_eligible=eligible, designs={})
            assignments = json.loads(row['sequence_first_tied_model_assignments'])
            assert [r['gene'] for r in assignments] == json.loads(row['sequence_first_nearest_genes'])
            for item in assignments:
                assert (item['model_id'], str(item['version'])) == models.get(item['gene'], ('', ''))
            if row['chosen_reference_gene']:
                assert models[row['chosen_reference_gene']] == (row['reference_model'], row['reference_version'])
            for design, genes in reference_sets(row):
                evaluated = []
                for gene in genes:
                    keys = {s: pair_key(ids, row['gene_' + s], gene) for s in ['a', 'b']}
                    mid, version = models.get(gene, ('', ''))
                    native_states = {g: coorthology(membership[g][keys['a']], membership[g][keys['b']])
                                     for g in ['profile', 'mafft']}
                    evaluated.append(dict(reference_gene=gene, reference_model=mid, reference_version=version,
                                          lexical_choice=gene == genes[0], native_pair_keys=keys,
                                          native_coorthology=native_states))
                    tie_count += 1
                    summary[f'{guide}|{design}|eligible={int(eligible)}|model={int(bool(mid))}|'
                            f'profile={native_states["profile"]}|mafft={native_states["mafft"]}'] += 1
                    for side in ['a', 'b']:
                        writer.writerow(dict(source_guide=guide, source_row_number=number, source_row_sha256=rh,
                                             family=row['family'], gene_node=row['gene_node'], taxon_id=row['taxon_id'],
                                             source_context_status=row['status'], sequence_first_status=row['sequence_first_status'],
                                             design=design, reference_gene=gene, reference_model=mid, reference_version=version,
                                             lexical_choice=int(gene == genes[0]), parent_context_eligible=int(eligible),
                                             duplicate_side=side, duplicate_gene=row['gene_' + side], native_pair_key=keys[side],
                                             profile_native_ortholog=membership['profile'][keys[side]],
                                             mafft_native_ortholog=membership['mafft'][keys[side]]))
                        link_count += 1
                context['designs'][design] = evaluated
                if not genes:
                    summary[f'{guide}|{design}|eligible={int(eligible)}|no_reference_gene'] += 1
            contexts.write(json.dumps(context, sort_keys=True, separators=(',', ':')) + '\n')
    verify(bindings)
    result = dict(status='complete_full_reference_native_orthology_pending_readback',
                  plan_sha256=sha(args.plan), source_receipt_sha256=sha(Path(plan['source']) / 'receipt.json'),
                  target_contexts=len(rows), guide_contexts=dict(counts), unique_gene_pair_queries=len(queries),
                  reference_tie_records=tie_count, duplicate_reference_links=link_count,
                  guides=stats, context_reference_native_summary=dict(sorted(summary.items())),
                  artifacts={p.name: sha(p) for p in out.iterdir()}, source_hashes=bindings,
                  compile_command=compile_command, scientific_eligibility=False,
                  scope='All original target contexts and fields, both reference designs, every nearest tie and '
                        'both duplicate sides retained, including unmodeled references and ineligible parents. '
                        'Unique physical gene queries deduplicated without deleting logical links. Both full '
                        'reciprocal native streams scanned. Native membership does not qualify parent context, '
                        'validate biological orthology/rooting, reconstruct ancestors or establish a duplication effect; '
                        'small-family supplements are excluded. Context/tie counts overlap across guides and designs.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
