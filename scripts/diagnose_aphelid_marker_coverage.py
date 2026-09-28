"""Trace aphelid coverage from BUSCO copy status to qualified alignment masks."""
from collections import Counter,defaultdict
import csv
import json
import math
from pathlib import Path
from ancestral_chain_attempt import sha,write_json


def rows(path):
    with open(path) as h:return list(csv.DictReader(h,delimiter='\t'))


def main():
    taxon='F1243177';base=Path('results/qc/busco-gene-copies-v1')
    receipt=json.loads((base/'receipt.json').read_text())
    for name,h in receipt['artifacts'].items():assert sha(base/name)==h
    raw=Path('results/busco')/taxon/'run_eukaryota_odb12.2/full_table.tsv'
    binding=next(r for r in receipt['sources'] if r['taxon_id']==taxon)
    assert sha(raw)==binding['busco_table_sha256']
    by_marker=defaultdict(list)
    for line in raw.read_text().splitlines():
        if line and not line.startswith('#'):
            cells=line.split('\t');by_marker[cells[0]].append(cells)
    assert len(by_marker)==125
    statuses={m:{r[1] for r in hits} for m,hits in by_marker.items()}
    assert all(len(s)==1 for s in statuses.values())
    counts=Counter(next(iter(s)) for s in statuses.values())
    assert counts==dict(Complete=1,Duplicated=119,Missing=5)
    copies=[r for r in rows(base/'marker_gene_copies.tsv') if r['taxon_id']==taxon]
    assert len(copies)==125
    for row in copies:
        assert row['busco_status']==next(iter(statuses[row['marker']]))
        if row['busco_status']=='Duplicated':assert int(row['distinct_resolved_gene_ids'])>1
    mapping=Path('results/phylogeny/markers-full-v1/protein_mapping.tsv')
    links=[r for r in rows(mapping) if r['taxon_id']==taxon]
    assert len(links)==1 and links[0]['marker']=='5001734at2759'
    paired=Path('results/phylogeny/paired-inputs-esmfold-all-completed-20260922-v1')
    pr=json.loads((paired/'receipt.json').read_text())
    assert sha(paired/'taxon_coverage.tsv')==pr['artifacts']['taxon_coverage.tsv']
    coverage=[r for r in rows(paired/'taxon_coverage.tsv') if r['taxon_id']==taxon]
    assert len(coverage)==125
    focal=next(r for r in coverage if r['marker']==links[0]['marker'])
    assert focal['model_name'].startswith('ESM-')
    assert int(focal['required_observed_columns'])==max(50,math.ceil(.3*int(focal['original_marker_columns'])))
    assert int(focal['observed'])==53 and int(focal['required_observed_columns'])==54 and focal['taxon_eligible']=='False'
    categories=['observed','noncanonical_or_missing_sequence','no_structural_mapping','invalid_native_feature','low_feature_plddt','high_feature_pae']
    assert sum(int(focal[k]) for k in categories)==int(focal['original_marker_columns'])
    assert all(r['taxon_eligible']=='False' for r in coverage)
    result=dict(status='complete_aphelid_marker_coverage_cause_trace',taxon_id=taxon,
        raw_busco_marker_counts=dict(counts),duplicated_markers_with_multiple_resolved_gene_ids=119,
        single_copy_mapping=links[0],sole_selected_marker_qualification=focal,
        sources={str(p):sha(p) for p in [raw,base/'receipt.json',base/'marker_gene_copies.tsv',mapping,paired/'receipt.json',paired/'taxon_coverage.tsv']},
        script_sha256=sha(__file__),next_required_work='Assess all copies using existing gene identities and copy-aware phylogenetic/structural analyses; do not select an arbitrary paralog or lower confidence eligibility to include this taxon.',
        scope='Trace through frozen single-copy marker and qualification pipelines. Duplicated hits are not proof of biological gene duplication, hybrid origin or assembly artifacts. Missing selected markers are not missing whole-proteome sequences. No thresholds or production selections changed.')
    write_json(Path('metadata/aphelid_marker_coverage_cause_20260928.json'),result)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
