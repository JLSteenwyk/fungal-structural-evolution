#!/usr/bin/env python3
"""Check every projected ASA row against source tables, maps and paired FASTAs."""
import argparse
import csv
import gzip
import json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha


def table(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        yield from csv.DictReader(f, delimiter='\t')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['projection', 'inputs', 'snapshot', 'accessibility', 'audit', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new readback receipt')
    receipts = {k: checked_receipt(getattr(a, k)) for k in ['projection', 'inputs', 'snapshot', 'audit']}
    pr = receipts['projection']; ar = receipts['audit']
    for name in ['inputs', 'snapshot', 'accessibility', 'audit']:
        if pr['source_receipts'][name] != sha(getattr(a, name) / 'receipt.json'):
            raise ValueError('Changed projection source: ' + name)
    if ar['status'] != 'passed_full_accessibility_snapshot' or ar['assessment_receipt_sha256'] != sha(a.accessibility / 'receipt.json') or ar['snapshot_receipt_sha256'] != sha(a.snapshot / 'receipt.json'):
        raise ValueError('Missing matching full accessibility audit')
    asa = json.loads((a.accessibility / 'receipt.json').read_text())
    audited = {r['model_id']: r for r in table(a.audit / 'audited_models.tsv')}
    links = {(r['marker'], r['taxon_id']): r for r in table(a.snapshot / 'marker_structure_links.tsv')}
    mapping = defaultdict(dict)
    for r in table(a.snapshot / 'matrix_to_structure_residues.tsv.gz'):
        key = r['marker'], r['taxon_id']; col = int(r['matrix_column_1based'])
        if col in mapping[key]:
            raise ValueError('Duplicate source map')
        mapping[key][col] = (r['model_id'], r['protein_id'], int(r['protein_residue_1based']), float(r['ca_plddt']))
    alignments = {}; columns = {}; seen = {}; expected_count = 0
    ready = [r for r in table(a.inputs / 'marker_summary.tsv') if r['status'] == 'ready_for_inference']
    for r in ready:
        marker = r['marker']; folder = a.inputs / marker
        aa = {x.id: str(x.seq) for x in SeqIO.parse(folder / 'aa.faa', 'fasta')}
        di = {x.id: str(x.seq) for x in SeqIO.parse(folder / '3di.faa', 'fasta')}
        if set(aa) != set(di):
            raise ValueError('Paired taxa differ')
        cols = list(table(folder / 'columns.tsv'))
        columns[marker] = [int(x['matrix_column_1based']) for x in cols]
        if [int(x['paired_column_1based']) for x in cols] != list(range(1, len(cols) + 1)):
            raise ValueError('Paired column grid differs')
        for taxon in aa:
            if len(aa[taxon]) != len(di[taxon]) or len(aa[taxon]) != len(cols):
                raise ValueError('Paired dimensions differ')
            if [c == '?' for c in aa[taxon]] != [c == '?' for c in di[taxon]]:
                raise ValueError('Paired masks differ')
            key = marker, taxon
            alignments[key] = aa[taxon], di[taxon]
            seen[key] = bytearray(len(cols))
            expected_count += sum(c != '?' for c in aa[taxon])

    @lru_cache(maxsize=32)
    def residues(sid):
        rp = a.accessibility / (sid + '.receipt.json')
        if sha(rp) != asa['entry_receipts'][sid] or sha(rp) != audited[sid]['entry_receipt_sha256']:
            raise ValueError('Changed model accessibility receipt')
        entry = json.loads(rp.read_text())
        path = a.accessibility / (sid + '.residues.tsv.gz')
        if sha(path) != entry['table_sha256']:
            raise ValueError('Changed raw residue table')
        rows = list(table(path)); result = {int(r['protein_residue_1based']): r for r in rows}
        if len(result) != len(rows):
            raise ValueError('Duplicate raw residue')
        return entry, result

    count = 0; used_models = set(); marker_counts = defaultdict(int)
    for r in table(a.projection / 'paired_site_accessibility.tsv.gz'):
        key = r['marker'], r['taxon_id']; i = int(r['paired_column_1based']) - 1
        aa, di = alignments[key]
        if not 0 <= i < len(aa) or seen[key][i] or aa[i] == '?':
            raise ValueError('Unexpected, duplicate or masked projected cell')
        col = columns[key[0]][i]
        if int(r['matrix_column_1based']) != col or r['amino_acid'] != aa[i] or r['3di_state'] != di[i]:
            raise ValueError('Projection/alignment correspondence differs')
        link = links[key]; mapped = mapping[key][col]
        if (r['model_id'], r['protein_id'], int(r['protein_residue_1based']), float(r['ca_plddt'])) != mapped:
            raise ValueError('Projection/source mapping differs')
        if r['model_id'] != link['model_id'] or r['protein_id'] != link['protein_id']:
            raise ValueError('Projection/link identity differs')
        entry, raw = residues(r['model_id'])
        if entry['model_sha256'] != link['model_sha256'] or entry['sequence_sha256'] != link['sequence_sha256']:
            raise ValueError('Raw accessibility/model provenance differs')
        source = raw[int(r['protein_residue_1based'])]
        for field in ['amino_acid', 'sasa_angstrom_squared', 'ca_plddt', 'heavy_atom_count', 'context']:
            if r[field] != source[field]:
                raise ValueError('Projected value differs from raw residue: ' + field)
        seen[key][i] = 1; count += 1; used_models.add(r['model_id']); marker_counts[key[0]] += 1
    if count != expected_count or count != pr['observed_sites_linked'] or len(used_models) != pr['models_used']:
        raise ValueError('Full observation/model grid differs')
    for key, flags in seen.items():
        if list(flags) != [int(c != '?') for c in alignments[key][0]]:
            raise ValueError('Missing projected observation')
    for r in table(a.projection / 'marker_summary.tsv'):
        if marker_counts[r['marker']] != int(r['observed_sites_linked']):
            raise ValueError('Marker count differs')
    result = {'status': 'passed_full_raw_accessibility_projection_readback', 'rows': count,
              'markers': len(ready), 'models': len(used_models), 'marker_taxon_cells': len(seen),
              'source_receipt_sha256': sha(a.projection / 'receipt.json'),
              'script_sha256': sha(Path(__file__)),
              'scope': 'Every projected cell exactly matches its raw ASA/confidence/atom/context fields, source mapping and paired AA/3Di identity. Full observed grid checked with duplicate/missing/masked-cell rejection. Original full coordinate audit remains required; no independent ASA integration or biological exposure validation.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
