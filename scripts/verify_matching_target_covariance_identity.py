#!/usr/bin/env python3
"""Verify target-node/residual redundancy from every original selection link."""
import argparse
from collections import defaultdict
import csv
from datetime import datetime,timezone
import gzip
import json
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-completion',type=Path,default=Path('metadata/full_matching_case_index_completed_20261002.json'))
    p.add_argument('--counts',type=Path,default=Path('results/phylogeny/full-expanded-model-inputs-20261002-v1/setting_counts.tsv'))
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
    c=json.loads(a.case_completion.read_text())
    assert c['status']=='complete_verified_full_matching_logical_case_index' and c['exact_process_journals_checked']==2
    archive=Path(c['full_hash_archive']);assert sha(archive)==c['full_hash_archive_sha256']
    proof=json.loads(archive.read_text())
    assert proof['status']=='complete_verified_full_matching_logical_case_index_archive' and len(proof['services'])==2
    assert proof['summary']['selected_records']==c['selected_records']
    links=Path(c['producer_receipt']).parent/'selection_case_links.tsv.gz';expected=proof['source_hashes'][str(links)]
    assert sha(links)==expected;del proof
    ids={};groups=defaultdict(set);last=0
    with gzip.open(links,'rt') as f:
        for number,r in enumerate(csv.DictReader(f,delimiter='\t'),1):
            assert int(r['source_row_ordinal'])==number
            token=ids.setdefault(r['target_id'],len(ids));group=r['guide'],r['policy'],r['scenario_id']
            assert token not in groups[group],(group,r['target_id'],'duplicate selected target')
            groups[group].add(token);last=number
    assert last==c['selected_records']
    rows=0
    with a.counts.open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            assert int(r['unique_targets'])==int(r['quality_retained_records']);rows+=1
    result=dict(status='verified_target_node_identity_kernel_in_every_original_matching_stratum',checked_utc=datetime.now(timezone.utc).isoformat(),
        closed_case_completion=str(a.case_completion),closed_case_completion_sha256=sha(a.case_completion),closed_case_archive=str(archive),closed_case_archive_sha256=c['full_hash_archive_sha256'],
        closed_selection_links=str(links),closed_selection_links_sha256=expected,selected_records_checked=last,original_strata=len(groups),distinct_original_target_ids=len(ids),duplicate_targets_within_stratum=0,
        supplemental_input_producer_count_rows=rows,supplemental_count_table_sha256=sha(a.counts),supplemental_counts_independent_acceptance_pending=True,
        mathematical_result='Each original stratum has one selected control per target, hence unit target-node incidence in distinct columns gives Z_target Z_target.T=I. Every mask/quality subset retains this identity. Under a uniform residual diagonal the target-node variance is confounded with residual scale and is not separately identified.',
        source_hashes={str(Path(__file__)):sha(__file__)},
        scope='Complete direct scan of closed original selection links, verified archive status/two-journal service count/summary and archive/link checksums. Supplemental count table agreement is not full independent input acceptance. This does not imply independent physical predictions/genes/backgrounds or identity of the global pooled-case kernel. Nonuniform diagonals and other weighting/loadings require separate qualification. No variances/effects estimated.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','selected_records_checked','original_strata','distinct_original_target_ids','duplicate_targets_within_stratum','supplemental_input_producer_count_rows']}))


if __name__=='__main__':main()
