#!/usr/bin/env python3
"""Independently reconstruct every reference link and binary-search native membership."""
import argparse
import csv
import hashlib
import json
import sqlite3
from collections import Counter
from pathlib import Path

from query_reference_orthology_v2 import load_sources, verify
from readback_background_orthology import membership
from run_ortholog_pair_guide_comparison import sha


def read_rows(root):
    for guide in ['profile', 'mafft']:
        with (root / (guide + '_sequence_first_references.tsv')).open() as handle:
            for index, row in enumerate(csv.DictReader(handle, delimiter='\t'), 1):
                yield guide, index, row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    plan = json.loads(args.plan.read_text())
    source_receipt, bindings = load_sources(plan)
    root, source = Path(plan['output']), Path(plan['source'])
    result = json.loads((root / 'receipt.json').read_text())
    assert result['status'] == 'complete_full_reference_native_orthology_pending_readback'
    assert result['plan_sha256'] == sha(args.plan)
    assert result['source_receipt_sha256'] == sha(source / 'receipt.json')
    assert result['source_hashes'] == {**bindings, str(args.plan): sha(args.plan)}
    bindings.update({str(root / k): v for k, v in result['artifacts'].items()})
    bindings[str(args.plan)] = sha(args.plan)
    bindings[str(root / 'receipt.json')] = sha(root / 'receipt.json')
    verify(bindings)
    genes, reference_genes = set(), set()
    for guide, number, row in read_rows(source):
        genes.update([row['gene_a'], row['gene_b']])
        reference_genes.update(json.loads(row['nearest_reference_genes']))
        reference_genes.update(json.loads(row['sequence_first_nearest_genes']))
    genes.update(reference_genes)
    stream_mapping = plan['streams'][0]
    taxa = {}
    for line in Path(stream_mapping['species_ids']).read_text().splitlines():
        native, label = line.split(': ', 1)
        assert native not in taxa
        taxa[native] = label.rsplit('.', 1)[0]
    assert len(taxa) == len(set(taxa.values()))
    ordinals = {}
    with Path(stream_mapping['sequence_ids']).open() as handle:
        for i, line in enumerate(handle):
            native, protein = line.rstrip('\n').split(': ', 1)
            gene = taxa[native.partition('_')[0]] + '_' + protein
            if gene in genes:
                assert gene not in ordinals and i < 16777216
                ordinals[gene] = i
    assert set(ordinals) == genes
    models = {}
    with sqlite3.connect('file:' + str(Path(plan['bridge']).resolve()) + '?mode=ro', uri=True) as conn:
        for taxon, protein, mid, version in conn.execute('SELECT taxon_id,protein_id,model_id,version FROM structures'):
            gene = taxon + '_' + protein
            if gene in reference_genes:
                assert gene not in models
                models[gene] = (mid, str(version))
    requested = set()
    for _, _, row in read_rows(source):
        for gene in set(json.loads(row['nearest_reference_genes']) + json.loads(row['sequence_first_nearest_genes'])):
            for duplicate in [row['gene_a'], row['gene_b']]:
                a, b = min(ordinals[gene], ordinals[duplicate]), max(ordinals[gene], ordinals[duplicate])
                assert a < b
                requested.add(format(a, '06x') + format(b, '06x'))
    query_keys = sorted(requested)
    assert (root / 'queries.hex').read_text() == ''.join(k + '\n' for k in query_keys)
    flags, native_summaries = {}, {}
    for stream in plan['streams']:
        guide = stream['guide']
        audit = json.loads(Path(stream['receipt']).read_text())
        n = audit['unique_unordered_pairs']
        assert Path(stream['stream']).stat().st_size == n * 28
        found = {}
        with (root / (guide + '_membership.tsv')).open() as exported, Path(stream['stream']).open('rb') as native:
            for key in query_keys:
                record = exported.readline().rstrip('\n').split('\t')
                value = membership(native, n, key.encode())
                assert record == [key, str(value)]
                found[key] = str(value)
            assert exported.readline() == ''
        stats = dict(stream_pairs=n, queries=len(query_keys), present=sum(map(int, found.values())),
                     absent=len(query_keys) - sum(map(int, found.values())))
        assert result['guides'][guide] == stats
        flags[guide], native_summaries[guide] = found, stats
        print(guide, json.dumps(stats), flush=True)
    summary, guide_counts = Counter(), Counter()
    contexts_checked = ties_checked = links_checked = 0
    context_identities = set()
    with (root / 'contexts.jsonl').open() as contexts, (root / 'reference_membership.tsv').open() as links:
        ledger = csv.DictReader(links, delimiter='\t')
        for guide, index, row in read_rows(source):
            source_hash = hashlib.sha256(json.dumps(row, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            observed = json.loads(contexts.readline())
            is_eligible = row['status'] in {'provisional_reference_available', 'no_modeled_nonfocal_sister'}
            expected = dict(source_guide=guide, source_row_number=index, source_row_sha256=source_hash,
                            source=row, parent_context_eligible=is_eligible, designs={})
            identity = (guide, row['family'], row['gene_node'], row['gene_a'], row['gene_b'])
            assert identity not in context_identities
            context_identities.add(identity)
            designs = [('availability', json.loads(row['nearest_reference_genes'])),
                       ('sequence_first', json.loads(row['sequence_first_nearest_genes']))]
            sequence_assignments = json.loads(row['sequence_first_tied_model_assignments'])
            assert [(a['gene'], a['model_id'], str(a['version'])) for a in sequence_assignments] == [
                (gene, *models.get(gene, ('', ''))) for gene in designs[1][1]]
            assert row['sequence_first_chosen_gene'] == (designs[1][1][0] if designs[1][1] else '')
            assert row['chosen_reference_gene'] == (designs[0][1][0] if designs[0][1] else '')
            assert int(row['sequence_first_ties']) == len(designs[1][1])
            if row['chosen_reference_gene']:
                assert (row['reference_model'], row['reference_version']) == models[row['chosen_reference_gene']]
            for design, refs in designs:
                assert refs == sorted(set(refs))
                evaluations = []
                if not refs:
                    summary[f'{guide}|{design}|eligible={int(is_eligible)}|no_reference_gene'] += 1
                for ref in refs:
                    assert ref.partition('_')[0] != row['taxon_id']
                    mid, version = models.get(ref, ('', ''))
                    keys = {}
                    for side in ['a', 'b']:
                        duplicate = row['gene_' + side]
                        low, high = sorted([ordinals[duplicate], ordinals[ref]])
                        assert low < high
                        key = format(low, '06x') + format(high, '06x')
                        keys[side] = key
                        exported = next(ledger)
                        expected_link = dict(source_guide=guide, source_row_number=str(index), source_row_sha256=source_hash,
                                             family=row['family'], gene_node=row['gene_node'], taxon_id=row['taxon_id'],
                                             source_context_status=row['status'], sequence_first_status=row['sequence_first_status'],
                                             design=design, reference_gene=ref, reference_model=mid, reference_version=version,
                                             lexical_choice=str(int(ref == refs[0])), parent_context_eligible=str(int(is_eligible)),
                                             duplicate_side=side, duplicate_gene=duplicate, native_pair_key=key,
                                             profile_native_ortholog=flags['profile'][key], mafft_native_ortholog=flags['mafft'][key])
                        assert exported == expected_link
                        links_checked += 1
                    states = {}
                    for native_guide in ['profile', 'mafft']:
                        a, b = flags[native_guide][keys['a']] == '1', flags[native_guide][keys['b']] == '1'
                        states[native_guide] = ('both' if a and b else 'only_a' if a else 'only_b' if b else 'neither')
                    evaluations.append(dict(reference_gene=ref, reference_model=mid, reference_version=version,
                                            lexical_choice=ref == refs[0], native_pair_keys=keys, native_coorthology=states))
                    summary[f'{guide}|{design}|eligible={int(is_eligible)}|model={int(bool(mid))}|'
                            f'profile={states["profile"]}|mafft={states["mafft"]}'] += 1
                    ties_checked += 1
                expected['designs'][design] = evaluations
            assert observed == expected
            contexts_checked += 1
            guide_counts[guide] += 1
        assert contexts.readline() == '' and next(ledger, None) is None
    assert contexts_checked == result['target_contexts'] == plan['resources']['target_contexts']
    assert guide_counts == result['guide_contexts'] == {r['guide']: r['candidates'] for r in source_receipt['guides']}
    assert ties_checked == result['reference_tie_records'] and links_checked == result['duplicate_reference_links'] == 2 * ties_checked
    assert len(query_keys) == result['unique_gene_pair_queries']
    assert dict(summary) == result['context_reference_native_summary']
    verify(bindings)
    proof = dict(status='passed_full_reference_native_orthology_readback', plan_sha256=sha(args.plan),
                 producer_receipt_sha256=sha(root / 'receipt.json'), target_contexts=contexts_checked,
                 unique_gene_pair_queries=len(query_keys), reference_tie_records=ties_checked,
                 duplicate_reference_links=links_checked, guides=native_summaries,
                 context_reference_native_summary=dict(sorted(summary.items())), source_hashes=bindings,
                 checker_sha256=sha(__file__), scientific_eligibility=False,
                 scope='All native protein ordinals, complete source fields, context identities, reference designs, ties, '
                       'model/version assignments, lexical choices, parent eligibility, duplicate sides and native coorthology '
                       'statuses reconstructed. Every unique queried pair checked by independent binary search in both '
                       'globally audited reciprocal streams, with bracketing for absence. Shares source I/O/proof loader '
                       'and fixed-record reader, not producer ledger/summary/merge logic. Does not validate biology, '
                       'tree rooting, small-family supplements, ancestor states or a structural duplication effect.')
    args.output.write_text(json.dumps(proof, indent=2) + '\n')
    print(json.dumps({k: v for k, v in proof.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
