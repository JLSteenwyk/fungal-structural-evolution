#!/usr/bin/env python3
"""Project audited isolated-chain ASA onto observed paired AA/3Di alignment sites."""
import argparse
import csv
import gzip
import json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['inputs', 'snapshot', 'accessibility', 'audit', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    receipts = {k: checked_receipt(getattr(a, k)) for k in ['inputs', 'snapshot', 'audit']}
    audit = receipts['audit']
    asa_receipt_path = a.accessibility / 'receipt.json'
    asa = json.loads(asa_receipt_path.read_text())
    if audit['status'] != 'passed_full_accessibility_snapshot' or audit['assessment_receipt_sha256'] != sha(asa_receipt_path) or audit['snapshot_receipt_sha256'] != sha(a.snapshot / 'receipt.json'):
        raise ValueError('Full matching accessibility audit required')
    if receipts['inputs']['source_receipts']['snapshot']['sha256'] != sha(a.snapshot / 'receipt.json'):
        raise ValueError('Paired input snapshot differs')
    if audit['models_audited'] != asa['models'] or audit['residues_audited'] != asa['residues']:
        raise ValueError('Audit coverage differs')
    audited = {r['model_id']: r for r in read_table(a.audit / 'audited_models.tsv')}
    if set(audited) != set(asa['entry_receipts']):
        raise ValueError('Audited model universe differs')
    links_list = read_table(a.snapshot / 'marker_structure_links.tsv')
    links = {(r['marker'], r['taxon_id']): r for r in links_list}
    if len(links) != len(links_list):
        raise ValueError('Duplicate model links')
    mapped = defaultdict(dict)
    with gzip.open(a.snapshot / 'matrix_to_structure_residues.tsv.gz', 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            key = r['marker'], r['taxon_id']; col = int(r['matrix_column_1based'])
            if col in mapped[key]:
                raise ValueError('Duplicate residue mapping')
            mapped[key][col] = r
    @lru_cache(maxsize=128)
    def load_model(sid):
        rp = a.accessibility / (sid + '.receipt.json')
        if sha(rp) != asa['entry_receipts'][sid] or sha(rp) != audited[sid]['entry_receipt_sha256']:
            raise ValueError('ASA entry changed since audit')
        r = json.loads(rp.read_text()); table = a.accessibility / (sid + '.residues.tsv.gz')
        if sha(table) != r['table_sha256']:
            raise ValueError('ASA table changed')
        with gzip.open(table, 'rt') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        return r, {int(x['protein_residue_1based']): x for x in rows}
    ready = [r for r in read_table(a.inputs / 'marker_summary.tsv') if r['status'] == 'ready_for_inference']
    if len(ready) != receipts['inputs']['ready_markers']:
        raise ValueError('Ready marker count differs')
    a.output.mkdir(parents=True); summaries = []; total = 0; models_used = set()
    fields = ['marker', 'taxon_id', 'protein_id', 'model_id', 'paired_column_1based', 'matrix_column_1based', 'protein_residue_1based', 'amino_acid', '3di_state', 'sasa_angstrom_squared', 'ca_plddt', 'heavy_atom_count', 'context']
    with gzip.open(a.output / 'paired_site_accessibility.tsv.gz', 'wt') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n'); w.writeheader()
        for marker_row in ready:
            marker = marker_row['marker']; folder = a.inputs / marker
            aa = {r.id: str(r.seq) for r in SeqIO.parse(folder / 'aa.faa', 'fasta')}
            di = {r.id: str(r.seq) for r in SeqIO.parse(folder / '3di.faa', 'fasta')}
            cols = read_table(folder / 'columns.tsv')
            if set(aa) != set(di) or len(aa) != int(marker_row['eligible_taxa']) or len(cols) != int(marker_row['retained_columns']):
                raise ValueError('Paired dimensions differ')
            if [int(x['paired_column_1based']) for x in cols] != list(range(1, len(cols) + 1)):
                raise ValueError('Paired column grid differs')
            n = 0; masked = 0
            for taxon in sorted(aa):
                if len(aa[taxon]) != len(cols) or len(di[taxon]) != len(cols):
                    raise ValueError('Sequence dimensions differ')
                link = links[marker, taxon]; sid = link['model_id']; entry, residues = load_model(sid)
                if entry['model_sha256'] != link['model_sha256'] or entry['sequence_sha256'] != link['sequence_sha256']:
                    raise ValueError('Model identity differs')
                models_used.add(sid)
                for i, col in enumerate(cols):
                    if (aa[taxon][i] == '?') != (di[taxon][i] == '?'):
                        raise ValueError('Paired masks differ')
                    if aa[taxon][i] == '?':
                        masked += 1; continue
                    source = mapped[marker, taxon][int(col['matrix_column_1based'])]
                    pos = int(source['protein_residue_1based']); residue = residues[pos]
                    if source['model_id'] != sid or residue['amino_acid'] != aa[taxon][i] or float(source['ca_plddt']) != float(residue['ca_plddt']):
                        raise ValueError('Residue/alignment identity differs')
                    if float(residue['ca_plddt']) < 70:
                        raise ValueError('Observed paired site below focal confidence threshold')
                    row = {k: residue[k] for k in ['amino_acid', 'sasa_angstrom_squared', 'ca_plddt', 'heavy_atom_count', 'context']}
                    row.update(marker=marker, taxon_id=taxon, protein_id=link['protein_id'], model_id=sid, paired_column_1based=i+1, matrix_column_1based=col['matrix_column_1based'], protein_residue_1based=pos, **{'3di_state': di[taxon][i]})
                    w.writerow(row); n += 1
            summaries.append({'marker': marker, 'taxa': len(aa), 'paired_columns': len(cols), 'observed_sites_linked': n, 'masked_cells_omitted': masked})
            total += n
            print(marker, n, flush=True)
    write_table(a.output / 'marker_summary.tsv', summaries)
    result = {'status': 'complete_audited_accessibility_projection', 'markers': len(ready), 'models_used': len(models_used), 'observed_sites_linked': total,
              'source_receipts': {k: sha(getattr(a, k) / 'receipt.json') for k in ['inputs', 'snapshot', 'accessibility', 'audit']},
              'script_sha256': sha(Path(__file__)),
              'interpretation': 'Continuous isolated-chain focal-residue ASA on every observed paired AA/3Di site; masked alignment cells omitted and counted. Native context confidence comes from the paired inputs. ASA is not normalized across amino acids and no categorical core/surface, interface, evolutionary transition or rate is inferred. Predicted domain placement and unobserved partners can affect ASA.',
              'artifacts': {f.name: sha(f) for f in a.output.iterdir()}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
