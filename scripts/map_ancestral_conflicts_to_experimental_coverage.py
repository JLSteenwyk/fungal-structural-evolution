#!/usr/bin/env python3
"""Attach frozen experimental-search and observed-CA coverage to conflict sites."""
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
from Bio import SeqIO


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    pins = {}
    def checked(root, name):
        root = Path(root)
        rp = root / 'receipt.json'
        p = root / name
        assert sha(p) == json.loads(rp.read_text())['artifacts'][name]
        pins.update({str(x): sha(x) for x in [rp, p]})
        return p
    annotations = json.loads(checked('results/ancestral/ancestral-conflict-case-models-20260927-v1',
                                    'extant_model_annotations.json').read_text())
    pairroot = 'results/experimental_structures/whole-domain-case-residue-pairs-20260927-v1'
    dispositions = {r['sequence_id']: r for r in json.loads(checked(pairroot, 'sequence_dispositions.json').read_text())}
    contexts = collections.defaultdict(list)
    with gzip.open(checked(pairroot, 'residue_correspondences.jsonl.gz'), 'rt') as stream:
        for row in map(json.loads, stream):
            contexts[row['sequence_id']].append(row)
    fasta = Path('results/structural_comparisons/case-independent-control-inputs-20260927-v1/missing_independent_sequences.faa')
    source = json.loads((Path(pairroot) / 'receipt.json').read_text())
    assert sha(fasta) == source['source_hashes'][str(fasta)]
    pins[str(fasta)] = sha(fasta)
    queries = {r.id: str(r.seq) for r in SeqIO.parse(fasta, 'fasta')}
    ca_root = Path('results/experimental_structures/whole-domain-case-ca-mapping-20260927-v1')
    ca_receipt = json.loads((ca_root / 'receipt.json').read_text())
    pins[str(ca_root / 'receipt.json')] = sha(ca_root / 'receipt.json')
    cache = {}
    rows = []
    for annotation in annotations:
        row = dict(annotation)
        sid = 'S' + row['sequence_sha256']
        assert hashlib.sha256(queries[sid].encode()).hexdigest() == row['sequence_sha256']
        position = row['protein_position']
        assert queries[sid][position - 1] == row['residue']
        assert len(contexts[sid]) == dispositions[sid]['contexts']
        matches = []
        for c in contexts[sid]:
            mapped = [s for q, s in c['query_subject_position_pairs'] if q == position]
            assert len(mapped) <= 1
            match = dict(entity_id=c['entity_id'], context_index=c['context_index'],
                         status='outside_aligned_context', observed_positions=[])
            if position in c['query_positions_aligned_to_gap']:
                match['status'] = 'query_residue_aligned_to_gap'
            if mapped:
                assert match['status'] != 'query_residue_aligned_to_gap'
                entry, entity = c['entity_id'].split('_')
                if entry not in cache:
                    rp = ca_root / (entry + '.receipt.json')
                    assert sha(rp) == ca_receipt['entry_receipt_sha256'][entry]
                    p = ca_root / (entry + '.residues.tsv.gz')
                    assert sha(p) == json.loads(rp.read_text())['table_sha256']
                    pins.update({str(x): sha(x) for x in [rp, p]})
                    with gzip.open(p, 'rt') as stream:
                        cache[entry] = list(csv.DictReader(stream, delimiter='\t'))
                observed = [r for r in cache[entry] if r['entity_id'] == entity and int(r['label_seq_id']) == mapped[0]]
                assert observed
                match.update(status='sequence_correspondence_with_coordinate_dispositions',
                             entity_position=mapped[0], observed_positions=observed)
            matches.append(match)
        row.update(experimental_search_hits=dispositions[sid]['hits'],
                   experimental_contexts=matches,
                   experimental_status='no_hits_in_frozen_search' if not matches else 'contexts_retained')
        rows.append(row)
    out = Path('results/ancestral/ancestral-conflict-experimental-coverage-20260927-v1')
    out.mkdir(exist_ok=False)
    p = out / 'conflict_experimental_coverage.json'
    p.write_text(json.dumps(rows, indent=2) + '\n')
    unique = {(r['sequence_sha256'], r['protein_position'], c['entity_id'], c['entity_position'])
              for r in rows for c in r['experimental_contexts'] if 'entity_position' in c}
    result = dict(status='complete_conflict_experimental_coverage_mapping', annotation_rows=len(rows),
        unique_query_positions=len({(r['sequence_sha256'], r['protein_position']) for r in rows}),
        unique_query_entity_position_correspondences=len(unique),
        no_hit_annotation_rows=sum(r['experimental_status']=='no_hits_in_frozen_search' for r in rows),
        pins=pins, script_sha256=sha(__file__), artifacts={p.name:sha(p)},
        scope='All extant conflict annotations retained. Search-context correspondence and observed coordinate dispositions only; absent search hits do not prove database absence. Repeated boundary/alignment contexts and experimental chains are dependent. No ancestral state, site function, homology beyond the recorded alignment, or mechanism is validated.')
    rp = out / 'receipt.json'
    rp.write_text(json.dumps(result, indent=2) + '\n')
    result.update(completed_receipt_path=str(rp), completed_receipt_sha256=sha(rp))
    Path('metadata/ancestral_conflict_experimental_coverage_completed_20260927.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['pins','artifacts','scope']}))


if __name__ == '__main__':
    main()
