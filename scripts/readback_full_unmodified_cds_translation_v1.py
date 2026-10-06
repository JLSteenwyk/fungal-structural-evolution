#!/usr/bin/env python3
"""Independently replay all original DNA translations from the pinned codebook."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import time

from independent_codon_translation_v1 import DNA, lookup, translate
from reference_measurement_union_sources import bind, verify


def sequence_digest(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def fasta_records(path):
    """Independent FASTA parser; preserve source sequence case and order."""
    path = Path(path)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
        identifier, chunks = None, []
        for line in handle:
            if line.startswith('>'):
                if identifier is not None:
                    yield identifier, ''.join(chunks)
                identifier, chunks = line[1:].split()[0], []
            elif line.strip():
                assert identifier is not None
                chunks.append(''.join(line.split()))
        if identifier is not None:
            yield identifier, ''.join(chunks)


def records(path):
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            yield json.loads(line)


def check_taxon(entry, root, codebook, report):
    tables = {code: lookup(codebook, code) for code in [1, 3, 4, 5, 6, 12, 16, 26]}
    proteins = {}
    for identifier, sequence in fasta_records(entry['proteome']):
        assert identifier not in proteins
        proteins[identifier] = sequence
    assert len(proteins) == entry['source_products']
    directory = Path(root) / entry['taxon_id']
    output_paths = [directory / name for name in ['target_diagnostics.jsonl.gz',
                    'changed_original_codons.jsonl.gz', 'product_diagnostics.jsonl.gz']]
    assert [str(p) for p in output_paths] == report['artifact_paths']
    expected_linked = defaultdict(list)
    target_counts, inherited_counts, changed_counts = Counter(), Counter(), Counter()
    n_targets = dna_bases = changed_total = 0
    differences = iter(records(output_paths[1]))
    originals = fasta_records(entry['original_cds']) if entry['original_cds'] else iter(())
    for n, triple in enumerate(zip_longest(originals, records(entry['targets']), records(output_paths[0])), 1):
        dna_record, original, observed = triple
        assert all(item is not None for item in triple), (entry['taxon_id'], n, 'target census')
        cid, dna = dna_record
        assert original['ordinal'] == n and original['cds_id'] == cid
        assert len(dna) == original['original_target_length']
        assert sequence_digest(dna) == original['original_target_dna_sha256']
        assert all(symbol.upper() in DNA for symbol in dna)
        pid = original['protein_id']
        assert pid is None or pid in proteins
        protein = proteins[pid] if pid is not None else None
        evidence = original['translation_evidence']
        inherited = None if evidence['translation_table'] in [None, ''] else int(evidence['translation_table'])
        roles = dict(inherited=inherited, snapshot_nuclear=entry['context']['nuclear']['code'] or None,
                     snapshot_mitochondrial=entry['context']['mitochondrial']['code'] or None)
        assert all(c is None or c in tables for c in roles.values())
        values, raw_values = {}, {}
        for code in sorted(set(c for c in roles.values() if c is not None)):
            if len(dna) % 3:
                value = dict(status='not_translated_non_triplet_length', code=code, translated_length=None,
                             translated_sha256=None, terminal_stop_removed=None, internal_stop_count=None)
                raw = None
            else:
                raw = translate(dna, tables[code])
                terminal_stop = bool(raw) and raw[-1] == '*'
                protein_translation = raw[:len(raw) - 1] if terminal_stop else raw
                if protein is None:
                    status = 'translated_no_linked_normalized_protein'
                elif protein_translation == protein:
                    status = 'translated_exact_match_to_normalized_protein'
                else:
                    status = 'translated_mismatch_to_normalized_protein'
                value = dict(status=status, code=code, translated_length=len(protein_translation),
                    translated_sha256=sequence_digest(protein_translation), terminal_stop_removed=terminal_stop,
                    internal_stop_count=protein_translation.count('*'))
            values[str(code)] = value
            raw_values[code] = raw
        site_counts = dict(snapshot_nuclear=0, snapshot_mitochondrial=0)
        if inherited is not None and raw_values[inherited] is not None:
            for role in ['snapshot_nuclear', 'snapshot_mitochondrial']:
                other = roles[role]
                if other is None or other == inherited:
                    continue
                for index in range(len(dna) // 3):
                    old, new = raw_values[inherited][index], raw_values[other][index]
                    if old == new:
                        continue
                    expected_site = dict(taxon_id=entry['taxon_id'], ordinal=n, cds_id=cid, protein_id=pid,
                        original_target_dna_sha256=original['original_target_dna_sha256'], role=role, inherited_code=inherited,
                        alternative_code=other, original_codon_index_1based=index + 1,
                        original_dna_start_1based=index * 3 + 1, codon=dna[index * 3:index * 3 + 3],
                        inherited_amino_acid=old, alternative_amino_acid=new,
                        within_normalized_protein_bounds=protein is not None and index < len(protein),
                        biological_mapping_admitted=False)
                    assert next(differences, None) == expected_site, (entry['taxon_id'], n, role, index)
                    site_counts[role] += 1
                    changed_counts[role] += 1
                    changed_total += 1
        expected = dict(ordinal=n, cds_id=cid, protein_id=pid,
            original_target_dna_sha256=original['original_target_dna_sha256'], original_dna_length=len(dna),
            linked_protein_sha256=sequence_digest(protein) if protein is not None else None,
            linked_protein_length=len(protein) if protein is not None else None,
            original_joined_status=original['joined_status'], inherited_translation_evidence=evidence,
            predeclared_code_roles=roles, translation_diagnostics=values, changed_codon_counts=site_counts,
            dna_modified=False, protein_modified=False, best_code_selected=False,
            scientific_eligibility=False, biological_codon_eligibility=False, original_source_target=original)
        assert observed == expected, (entry['taxon_id'], n, 'target diagnostics')
        summary = {key: expected[key] for key in ['ordinal', 'cds_id', 'protein_id', 'original_target_dna_sha256',
            'original_dna_length', 'linked_protein_sha256', 'linked_protein_length', 'original_joined_status',
            'predeclared_code_roles', 'translation_diagnostics', 'changed_codon_counts']}
        if pid is not None:
            expected_linked[pid].append(summary)
        for role, code in roles.items():
            state = 'unspecified_code_not_tested' if code is None else values[str(code)]['status']
            target_counts[role + ':' + state] += 1
        new_inherited = 'unspecified_code_not_tested' if inherited is None else values[str(inherited)]['status']
        inherited_counts[str(evidence['inherited_status']) + ':' + new_inherited] += 1
        n_targets += 1
        dna_bases += len(dna)
    assert next(differences, None) is None, 'extra changed codon rows'
    assert n_targets == entry['target_records']
    n_products = selected = 0
    classifications, cross = Counter(), Counter()
    for triple in zip_longest(proteins.items(), records(entry['products']), records(output_paths[2])):
        source, original, observed = triple
        assert all(item is not None for item in triple), 'product census'
        pid, protein = source
        assert original['protein_id'] == pid and original['protein_length'] == len(protein)
        assert original['sequence_sha256'] == sequence_digest(protein)
        observations = expected_linked.get(pid, [])
        assert {t['ordinal'] for t in original['original_target_evidence']} == {t['ordinal'] for t in observations}
        expected = dict(original_source_product=original, original_target_diagnostics=observations,
            source_target_replaced=False, derived_target_substituted=False, best_code_selected=False,
            protein_modified=False, scientific_eligibility=False, biological_codon_eligibility=False)
        assert observed == expected, (entry['taxon_id'], pid, 'product diagnostics')
        availability = original['representative_availability']
        category = availability['availability'] if availability is not None else 'alternative_not_assigned'
        for role in ['inherited', 'snapshot_nuclear', 'snapshot_mitochondrial']:
            states = set()
            for observation in observations:
                code = observation['predeclared_code_roles'][role]
                states.add('unspecified_code_not_tested' if code is None else
                           observation['translation_diagnostics'][str(code)]['status'])
            state = '|'.join(sorted(states)) if states else 'no_original_target'
            classifications[role + ':' + state] += 1
            cross[original['joined_status'] + ':' + category + ':' + role + ':' + state] += 1
        n_products += 1
        selected += int(original['selected_representative'])
    expected_report = dict(taxon_id=entry['taxon_id'], source_products=n_products,
        selected_representatives=selected, target_records=n_targets, original_dna_bases=dna_bases,
        changed_codon_rows=changed_total, target_role_status_counts=dict(target_counts),
        inherited_status_cross_counts=dict(inherited_counts), changed_codon_role_counts=dict(changed_counts),
        product_role_status_counts=dict(classifications), coding_model_translation_cross_counts=dict(cross),
        artifact_paths=[str(p) for p in output_paths])
    assert expected_report == report, (entry['taxon_id'], 'taxon aggregate')
    assert n_products == entry['source_products'] and selected == entry['selected_representatives']
    return expected_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['plan', 'producer', 'producer-transport', 'receipt']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert not args.receipt.exists()
    plan, producer, transport = [json.loads(p.read_text()) for p in [args.plan, args.producer, args.producer_transport]]
    assert producer['status'] == 'complete_full_unmodified_cds_fixed_code_translation_pending_independent_readback'
    assert transport['original_tool_terminal_exit_code'] == 0
    assert transport['validation_sha256'] == sequence_digest(args.producer.read_text())
    verify(transport['source_hashes'])
    verify(plan['pins'])
    verify(producer['source_hashes'])
    assert plan['expected_taxa'] == producer['taxa'] == 526
    reports = {r['taxon_id']: r for r in producer['taxa_reports']}
    assert len(reports) == 526 and set(reports) == {e['taxon_id'] for e in plan['entries']}
    codebook = json.loads(Path(plan['context_receipt']).read_text())['codebook']
    root = Path(plan['reader_output'])
    root.mkdir(exist_ok=False)
    start, checked = time.monotonic(), []
    with ProcessPoolExecutor(max_workers=plan['cpu']) as pool:
        futures = [pool.submit(check_taxon, entry, plan['output'], codebook, reports[entry['taxon_id']])
                   for entry in plan['entries']]
        for future in as_completed(futures):
            checked.append(future.result())
            state = dict(checked_taxa=len(checked), expected_taxa=526, elapsed_seconds=time.monotonic() - start)
            temporary = root / 'state.partial'
            temporary.write_text(json.dumps(state, indent=2) + '\n')
            temporary.replace(root / 'state.json')
            print(json.dumps(state), flush=True)
    assert len(checked) == 526
    fields = ['target_role_status_counts', 'inherited_status_cross_counts', 'changed_codon_role_counts',
              'product_role_status_counts', 'coding_model_translation_cross_counts']
    totals = {key: sum(r[key] for r in checked) for key in ['source_products', 'selected_representatives',
              'target_records', 'original_dna_bases', 'changed_codon_rows']}
    assert all(value == producer[key] for key, value in totals.items())
    aggregates = {key: Counter() for key in fields}
    for report in checked:
        for key in fields:
            aggregates[key].update(report[key])
    assert all(dict(value) == producer[key] for key, value in aggregates.items())
    pins = dict(producer['source_hashes'])
    for path in [args.plan, args.producer, args.producer_transport, Path(__file__),
                 Path('scripts/independent_codon_translation_v1.py')]:
        bind(pins, path)
    verify(pins)
    result = dict(status='passed_full_unmodified_cds_fixed_code_translation_independent_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(), taxa=526, **totals,
        **{key: dict(value) for key, value in aggregates.items()}, source_hashes=pins,
        elapsed_seconds=time.monotonic() - start, scientific_eligibility=False,
        biological_codon_eligibility=False, genetic_code_admission=False, best_code_selected=False,
        gpu=False, new_predictions=0, all_eight_aims_incomplete=True,
        scope='Complete original DNA/protein/source/output census and independent FASTA/codebook/IUPAC '
              'retranslation of every target, every changed codon and every original product. No Bio or '
              'producer translation imports. Shared pinned source DNA, protein normalization and NCBI '
              'codebook remain dependencies. No code/compartment, codon alignment, selection or biology admission.')
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
