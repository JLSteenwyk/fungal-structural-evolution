#!/usr/bin/env python3
"""Measure complete expanded partitions without materializing candidate pairs."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import time
from orthofinder.tools import mcl


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--merged', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    start = time.time()
    receipt_path = a.merged / 'receipt.json'
    receipt_hash = sha(receipt_path)
    receipt = json.loads(receipt_path.read_text())
    a.output.mkdir(parents=True, exist_ok=False)
    summaries = []
    for guide in receipt['guides']:
        for name, digest in guide['artifacts'].items():
            if sha(a.merged / name) != digest:
                raise ValueError('Changed merged partition artifact: ' + name)
        name = guide['guide']
        groups = mcl.GetPredictedOGs(str(a.merged / name / 'clusters_id_pairs.txt'))
        with (a.merged / name / 'family_sources.tsv').open() as handle:
            crosswalk = list(csv.DictReader(handle, delimiter='\t'))
        if len(groups) != len(crosswalk):
            raise ValueError('Family universe mismatch')
        totals = Counter()
        taxa_seen = set()
        table = a.output / (name + '_family_workload.tsv')
        fields = ['family', 'source_type', 'genes', 'taxa', 'maximum_copies_per_taxon',
                  'unordered_cross_species_candidate_pairs', 'occupied_unordered_taxon_pairs']
        top = []
        with table.open('w') as handle:
            writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for i, genes in enumerate(groups):
                source = crosswalk[i]
                family = f'OG{i:07d}'
                digest = hashlib.sha256(('\n'.join(sorted(genes)) + '\n').encode()).hexdigest()
                if (source['new_family'] != family or source['membership_sha256'] != digest
                        or int(source['proteins']) != len(genes) or not genes):
                    raise ValueError('Native membership differs from crosswalk: ' + family)
                copies = Counter(g.split('_', 1)[0] for g in genes)
                taxa_seen.update(copies)
                n, k = len(genes), len(copies)
                pairs = (n*n - sum(v*v for v in copies.values())) // 2
                seen = independent = 0
                for value in copies.values():
                    independent += seen * value
                    seen += value
                if pairs != independent:
                    raise ValueError('Candidate-pair arithmetic mismatch')
                row = dict(family=family, source_type=source['source_type'], genes=n,
                           taxa=k, maximum_copies_per_taxon=max(copies.values()),
                           unordered_cross_species_candidate_pairs=pairs,
                           occupied_unordered_taxon_pairs=k*(k-1)//2)
                writer.writerow(row)
                totals['genes'] += n
                totals['families'] += 1
                totals['tree_eligible_families'] += n >= 3
                totals['single_taxon_families'] += k == 1
                totals['unordered_cross_species_candidate_pairs'] += pairs
                totals['occupied_unordered_family_taxon_pairs'] += k*(k-1)//2
                totals[source['source_type'] + '_candidate_pairs'] += pairs
                top.append(row)
                if len(top) >= 100:
                    top = sorted(top, key=lambda r: (-r['unordered_cross_species_candidate_pairs'], r['family']))[:25]
        top = sorted(top, key=lambda r: (-r['unordered_cross_species_candidate_pairs'], r['family']))[:25]
        if totals['genes'] != 5815847 or len(taxa_seen) != 526:
            raise ValueError('Expanded full-protein scope mismatch')
        directed = 2 * totals['unordered_cross_species_candidate_pairs']
        summary = dict(guide=name, **totals, taxa=len(taxa_seen),
                       directed_cross_species_candidate_pairs=directed,
                       flat_directed_pair_storage_scenarios=[dict(bytes_per_pair=b, gib=directed*b/2**30) for b in [64,128,256]],
                       top25_families=top, artifacts={table.name:sha(table)})
        summaries.append(summary)
        print(json.dumps({k:v for k,v in summary.items() if k != 'top25_families'}), flush=True)
        del groups, crosswalk
    if sha(receipt_path) != receipt_hash:
        raise ValueError('Merge receipt changed during census')
    result = dict(status='completed_expanded_reconciliation_workload', guides=summaries,
                  merge_receipt_sha256=receipt_hash, script_sha256=sha(Path(__file__)),
                  elapsed_seconds=time.time()-start,
                  interpretation='Exact full-partition candidate dimensions, not inferred orthology. Candidate pairs are not materialized. Storage scenarios describe hypothetical flat directed pair tables, not native grouped output predictions. A final runtime/memory execution plan remains required.')
    (a.output/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
