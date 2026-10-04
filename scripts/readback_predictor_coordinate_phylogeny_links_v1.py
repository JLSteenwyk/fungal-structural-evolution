#!/usr/bin/env python3
"""Rebuild all tree-path projections as taxon sets and verify complete control joins."""
import argparse
from collections import Counter,defaultdict
import csv
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
from pathlib import Path

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify

METRICS=['ca_superposition_rmsd_angstrom','all_distance_rms_change_angstrom',
         'local_distance_rms_change_angstrom','sequence_local_distance_rms_change_angstrom']
STATUSES=['insufficient_matched_taxa','no_observations_on_one_side','terminal_projection',
          'unique_internal_projection','shared_internal_projection']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['producer','transport','receipt']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists()
    prod=json.loads(a.producer.read_text());transport=json.loads(a.transport.read_text())
    assert prod['status']=='complete_full_predictor_coordinate_phylogeny_links_pending_independent_readback'
    assert transport['validation_sha256']==sha(a.producer) and transport['original_tool_terminal_exit_code']==0
    pins=dict(prod['source_hashes'])
    for q in [a.producer,a.transport,Path(__file__)]:bind(pins,q)
    verify(pins)
    root=Path('results/phylogeny/matched-predictor-branch-inputs-20261003-v1')
    axes=json.loads((root/'input_axes.json').read_text());positions=axes['positions']
    with gzip.open(root/'comparison_cases.jsonl.gz','rt') as f:cases=[json.loads(line) for line in f]
    expected_views=[]
    for view in axes['views']:
        original=[]
        for branch in view['branches']:
            mask=int(branch['split_mask_hex'],16)
            original.append({taxon for taxon,position in positions.items() if mask&(1<<position)})
        expected_views.append(original)
    contexts={(r['marker'],r['taxon']):r for r in csv.DictReader(Path('results/phylogeny/paired-source-model-context-20260926-v1/source_model_context.tsv').open(),delimiter='\t')}
    primary={r['model_pair_id']:r for r in csv.DictReader(Path('results/phylogeny/overlap-coordinate-benchmark-20261004-v1/all_model_pair_coordinates.tsv').open(),delimiter='\t') if r['plddt_cutoff']=='70' and r['pae_cutoff']=='10'}
    distributions={}
    for r in csv.DictReader(Path('results/phylogeny/matched-predictor-distribution-summary-20261004-v1/paired_branch_distributions.tsv').open(),delimiter='\t'):
        distributions[r['input_id'],r['split_mask_hex'],r['model'],r['mode']]=r
    aa={}
    for cid in {c['input_id'] for c in cases if c['input_id']}:
        role=json.loads((Path('results/phylogeny/matched-predictor-branch-fits-20261003-v1/roles')/(cid+'-aa.json')).read_text())
        aa[cid]={b['split_mask_hex']:b for b in role['point_estimate']['branches']}
    out=Path('results/phylogeny/predictor-coordinate-phylogeny-links-20261004-v1')
    case_counts=Counter();branch_counts=Counter();path_count=0;slots=0;links_count=0;shared=0;single=0;seen=set()
    with gzip.open(out/'all_case_coordinate_links.jsonl.gz','rt') as cf,gzip.open(out/'all_internal_path_coordinate_links.jsonl.gz','rt') as pf:
        for source in cases:
            saved=json.loads(next(cf));view=axes['views'][source['view_index']];taxa=set(source['taxa'])
            assert taxa<=set(view['taxa'])
            tag=(source['marker'],source['view_index']);assert tag not in seen;seen.add(tag)
            for k in ['marker','view_index','cohort','view','input_id','taxa','original_internal_branches']:
                assert saved[k]==source[k]
            assert saved['case_status']==source['status'] and saved['scientific_eligibility'] is False
            assert saved['coordinate_mask']==dict(plddt=70,pae=10,scope='Whole-protein confidence mask; not marker-fit columns')
            counter={k:0 for k in STATUSES};groups=defaultdict(list)
            for i,original in enumerate(expected_views[source['view_index']]):
                side=original&taxa;other=taxa-side
                if len(taxa)<4:counter['insufficient_matched_taxa']+=1
                elif not side or not other:counter['no_observations_on_one_side']+=1
                elif min(len(side),len(other))==1:counter['terminal_projection']+=1
                else:
                    canon=min([side,other],key=lambda s:(len(s),sum(1<<positions[t] for t in s)))
                    groups[sum(1<<positions[t] for t in canon)].append(i)
            for indices in groups.values():
                counter['unique_internal_projection' if len(indices)==1 else 'shared_internal_projection']+=len(indices)
            assert counter==source['branch_status_counts']==saved['branch_status_counts']
            assert sum(counter.values())==source['original_internal_branches']==len(view['branches'])
            assert saved['mapped_internal_paths']==len(groups)
            slots+=sum(counter.values());branch_counts.update(counter);case_counts[source['status']]+=1
            links=saved['matched_taxon_coordinate_links'];assert [link['taxon'] for link in links]==source['taxa']
            assert len(links)==len(taxa)
            for link in links:
                ctx=contexts[source['marker'],link['taxon']]
                key=tuple(ctx[k] for k in ['reference_model_id','reference_model_version','local_model_id','local_model_version'])
                pid=hashlib.sha256(json.dumps(key,separators=(',',':')).encode()).hexdigest();coord=primary[pid]
                assert link['model_pair_id']==pid and link['coordinate_status']==coord['status']
                assert link['retained_residues']==int(coord['retained_residues'])
                for metric in METRICS:assert link[metric]==(float(coord[metric]) if coord[metric] else None)
                links_count+=1
            for mask,indices in sorted(groups.items()):
                path=json.loads(next(pf));split=hex(mask)
                for k in ['marker','view_index','cohort','view','input_id']:assert path[k]==source[k]
                assert path['observed_split_mask_hex']==split and path['original_branch_indices']==indices
                assert path['original_branch_count']==len(indices)
                if len(indices)==1:single+=1;projection='unique_internal_projection'
                else:shared+=1;projection='shared_internal_projection'
                assert path['projection']==projection
                expected_mapping=next(m for m in source['internal_branch_mapping'] if m['split_mask_hex']==split)
                assert expected_mapping['original_branch_indices']==indices
                units={view['branches'][i]['branch_length_unit'] for i in indices};assert len(units)==1
                assert path['original_path_length_unit']==next(iter(units))
                assert path['original_path_length_sum']==math.fsum(view['branches'][i]['branch_length'] for i in indices)
                assert path['aa_point_branch_length']==aa[source['input_id']][split]['length']
                assert path['aa_point_unit']=='LG+F+G4 expected amino-acid substitutions/site'
                side=[link for link in links if mask&(1<<positions[link['taxon']])]
                other=[link for link in links if not mask&(1<<positions[link['taxon']])]
                assert min(len(side),len(other))>=2
                for label,entries in [('side',side),('complement',other)]:
                    assert path[label+'_taxa']==[link['taxon'] for link in entries]
                    assert path[label+'_model_pair_ids']==[link['model_pair_id'] for link in entries]
                    assert path[label+'_coordinate_compared']==sum(link['coordinate_status']=='coordinates_compared' for link in entries)
                records=path['paired_predictor_distributions'];assert len(records)==6
                expected={(m,mode) for m in ['af','af_empirical','llm'] for mode in ['iid_sites','circular_blocks10']}
                assert {(r['model'],r['mode']) for r in records}==expected
                for r in records:
                    original=distributions[source['input_id'],split,r['model'],r['mode']]
                    assert original['internal']=='True' and all(r[k]==original[k] for k in r)
                    assert r['calibrated_effect_interval']==r['scientific_eligibility']=='False'
                assert path['branch_localization_accepted'] is False and path['scientific_eligibility'] is False
                path_count+=1
        assert next(cf,None) is None and next(pf,None) is None
    assert seen=={(m,i) for m in axes['markers'] for i in range(70)}
    assert len(seen)==8750 and slots==4523750 and path_count==31290
    assert dict(case_counts)==prod['case_status_counts'] and dict(branch_counts)==prod['branch_status_counts']
    assert links_count==prod['matched_taxon_coordinate_link_occurrences'];verify(pins)
    result=dict(status='passed_full_independent_predictor_coordinate_phylogeny_link_readback',checked_utc=datetime.now(timezone.utc).isoformat(),
        producer_receipt=str(a.producer),producer_receipt_sha256=sha(a.producer),comparison_cases=8750,
        original_internal_branch_slots=4523750,mapped_internal_paths=31290,single_original_branch_paths=single,
        shared_original_branch_paths=shared,matched_taxon_coordinate_link_occurrences=links_count,
        paired_distribution_occurrences=187740,source_hashes=pins,scientific_eligibility=False,
        all_eight_aims_incomplete=True,
        scope='Every original branch reconstructed independently with taxon-set intersections/complements, canonical '
              'observed bipartitions and original-index groups; all lost/terminal/insufficient/unique/shared dispositions '
              'and all coordinate/AA/paired-resampling joins checked. No producer branch-mapping helper used. Both-side '
              'model identities, coverage failures, original branch-length units and merged paths remain explicit. '
              'No geometry-based case selection, single-branch attribution, physical-rate conversion or calibrated inference.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
