"""Inventory exact-sequence model availability for every aphelid BUSCO copy."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from ancestral_chain_attempt import write_json
from readback_whole_proteome_catalog import sha


def main():
    root = Path('results/ecology/aphelid-marker-copy-inventory-20260928-v1')
    receipt = json.loads((root / 'receipt.json').read_text())
    copies = root / 'marker_protein_copies.tsv'
    assert sha(copies) == receipt['artifacts'][copies.name]
    with copies.open() as h:
        rows = list(csv.DictReader(h, delimiter='\t'))
    wanted = {r['sequence_sha256'] for r in rows}
    assert len(rows) == 241 and len(wanted) == 235
    sources = {str(root / 'receipt.json'): sha(root / 'receipt.json'), str(copies): sha(copies)}
    models = defaultdict(lambda: defaultdict(list))
    scans = {}
    af = Path('results/structures/whole-proteome-afdb-catalog-20260928-v1')
    audit_path = Path(str(af) + '-readback.json')
    audit = json.loads(audit_path.read_text())
    assert audit['status'] == 'passed_full_proteome_sequence_and_model_selection_readback'
    assert audit['producer_receipt_sha256'] == sha(af / 'receipt.json')
    sources[str(audit_path)] = sha(audit_path)
    esm = Path('results/structures/esmfold-all-completed-20260922-v1')
    for source, base, name in [('afdb', af, 'models.jsonl'), ('esmfold', esm, 'inventory.jsonl')]:
        rec = json.loads((base / 'receipt.json').read_text())
        path = base / name
        digest = sha(path)
        assert digest == rec['artifacts'][name]
        sources[str(path)] = digest
        sources[str(base / 'receipt.json')] = sha(base / 'receipt.json')
        count = 0
        with path.open() as handle:
            for line in handle:
                record = json.loads(line)
                if source == 'esmfold':
                    assert record['status'] == 'verified'
                    batch = record['models']
                else:
                    batch = [record]
                for model in batch:
                    count += 1
                    seq = model['sequence_sha256']
                    if seq in wanted:
                        assert sha(model['path']) == model['sha256']
                        models[seq][source].append(model)
        assert count == (rec['unique_models'] if source == 'afdb' else rec['models'])
        scans[source] = count
    records, counts = [], Counter()
    for row in rows:
        by_source = models[row['sequence_sha256']]
        for batch in by_source.values():
            for model in batch:
                assert int(model['length']) == int(row['length'])
        available = [source for source in ('afdb', 'esmfold') if by_source[source]]
        status = '+'.join(available) if available else 'absent_from_both_frozen_catalogs'
        counts[status] += 1
        records.append(dict(marker=row['marker'], protein_id=row['protein_id'],
                            sequence_sha256=row['sequence_sha256'], length=row['length'],
                            availability=status,
                            afdb_models_json=json.dumps(by_source['afdb'], sort_keys=True),
                            esmfold_models_json=json.dumps(by_source['esmfold'], sort_keys=True)))
    out = Path('results/ecology/aphelid-copy-structure-coverage-20260928-v1')
    out.mkdir(parents=True, exist_ok=False)
    table = out / 'all_copy_model_availability.tsv'
    with table.open('w') as h:
        w = csv.DictWriter(h, fieldnames=list(records[0]), delimiter='\t')
        w.writeheader(); w.writerows(records)
    with table.open() as h:
        assert list(csv.DictReader(h, delimiter='\t')) == records
    result = dict(status='complete_exact_sequence_aphelid_copy_model_availability',
                  proteins=len(rows), distinct_sequences=len(wanted),
                  scanned_models=scans, protein_availability_counts=dict(counts),
                  distinct_sequences_with_any_model=sum(bool(v['afdb'] or v['esmfold']) for v in models.values()),
                  sources=sources, script_sha256=sha(__file__),
                  artifacts={table.name: sha(table)},
                  scope='All marker copy identities retained; exact sequence and length matching, selected coordinate hashes checked. Frozen-catalog absence is not public absence. Availability is not confidence-qualified coverage, and no new inference is launched.')
    write_json(out / 'receipt.json', result)
    write_json(Path('metadata/aphelid_copy_structure_coverage_20260928.json'), result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
