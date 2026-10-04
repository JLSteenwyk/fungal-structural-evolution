#!/usr/bin/env python3
"""Link full predictor coordinate controls to every original marker/tree-view disposition."""
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
from pathlib import Path

from ancestral_chain_attempt import sha
from matched_predictor_branch_inputs import branch_mapping
from reference_measurement_union_sources import bind,verify

METRICS=['ca_superposition_rmsd_angstrom','all_distance_rms_change_angstrom',
         'local_distance_rms_change_angstrom','sequence_local_distance_rms_change_angstrom']


def identifier(row):
    key=tuple(row[k] for k in ['reference_model_id','reference_model_version','local_model_id','local_model_version'])
    return hashlib.sha256(json.dumps(key,separators=(',',':')).encode()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists() and not a.receipt.exists();pins={}
    def read(q):bind(pins,q);return json.loads(Path(q).read_text())
    complete=read('metadata/matched_predictor_branch_inputs_completed_20261003_v1.json')
    assert complete['status']=='complete_verified_full_matched_predictor_branch_inputs'
    archive=read(complete['full_hash_archive']);assert sha(complete['full_hash_archive'])==complete['full_hash_archive_sha256']
    for q,h in archive['source_hashes'].items():bind(pins,q,h)
    coords=read('metadata/overlap_coordinate_benchmark_20261004_v1.json')
    coordinate_reader=read('metadata/overlap_coordinate_readback_20261004_v1.json')
    coordinate_transport=read('metadata/overlap_coordinate_readback_transport_20261004_v1.json')
    assert coordinate_reader['status']=='passed_full_independent_matched_predictor_coordinate_readback'
    assert coordinate_reader['producer_receipt_sha256']==sha('metadata/overlap_coordinate_benchmark_20261004_v1.json')
    assert coordinate_transport['validation_sha256']==sha('metadata/overlap_coordinate_readback_20261004_v1.json')
    assert coordinate_transport['original_tool_terminal_exit_code']==0
    for q,h in coordinate_reader['source_hashes'].items():bind(pins,q,h)
    distributions=read('metadata/matched_predictor_distribution_summary_20261004_v1.json')
    distribution_reader=read('metadata/matched_predictor_distribution_readback_20261004_v1.json')
    assert distribution_reader['status']=='passed_full_independent_paired_predictor_distribution_summary_readback'
    assert distribution_reader['source_receipt_sha256']==sha('metadata/matched_predictor_distribution_summary_20261004_v1.json')
    for q,h in distribution_reader['source_hashes'].items():bind(pins,q,h)
    for q in [Path(__file__),Path('scripts/matched_predictor_branch_inputs.py')]:bind(pins,q)
    root=Path('results/phylogeny/matched-predictor-branch-inputs-20261003-v1')
    axes=read(root/'input_axes.json')
    cases_path=root/'comparison_cases.jsonl.gz';bind(pins,cases_path)
    with gzip.open(cases_path,'rt') as f:cases=[json.loads(line) for line in f]
    assert len(cases)==8750 and len(axes['views'])==70 and len(axes['markers'])==125 and len(axes['positions'])==526
    context_path=Path('results/phylogeny/paired-source-model-context-20260926-v1/source_model_context.tsv');bind(pins,context_path)
    contexts={(r['marker'],r['taxon']):r for r in csv.DictReader(context_path.open(),delimiter='\t')};assert len(contexts)==673
    primary={r['model_pair_id']:r for r in csv.DictReader(Path(coords['table']).open(),delimiter='\t')
             if r['plddt_cutoff']=='70' and r['pae_cutoff']=='10'};assert len(primary)==643
    distribution_path=Path('results/phylogeny/matched-predictor-distribution-summary-20261004-v1/paired_branch_distributions.tsv');bind(pins,distribution_path)
    paired={}
    for r in csv.DictReader(distribution_path.open(),delimiter='\t'):
        key=(r['input_id'],r['split_mask_hex'],r['model'],r['mode']);assert key not in paired;paired[key]=r
    assert len(paired)==13170
    aa={}
    for input_id in sorted({c['input_id'] for c in cases if c['input_id']}):
        q=Path('results/phylogeny/matched-predictor-branch-fits-20261003-v1/roles')/(input_id+'-aa.json')
        role=read(q);assert role['input_id']==input_id and role['role']=='aa' and role['exit_code']==0
        aa[input_id]={b['split_mask_hex']:b for b in role['point_estimate']['branches']}
    assert len(aa)==133;verify(pins);a.output.mkdir()
    case_path=a.output/'all_case_coordinate_links.jsonl.gz';path_path=a.output/'all_internal_path_coordinate_links.jsonl.gz'
    counts=Counter();branch_counts=Counter();path_count=0;branch_slots=0;matched_links=0;seen=set()
    with gzip.open(case_path,'xt',compresslevel=6) as cf,gzip.open(path_path,'xt',compresslevel=6) as pf:
        for case in cases:
            key=(case['marker'],case['view_index']);assert key not in seen;seen.add(key)
            view=axes['views'][case['view_index']];assert (view['cohort'],view['view'])==(case['cohort'],case['view'])
            mapping,expected_counts=branch_mapping(view,case['taxa'],axes['positions'])
            assert mapping==case['internal_branch_mapping'] and expected_counts==case['branch_status_counts']
            branch_slots+=case['original_internal_branches'];branch_counts.update(expected_counts);counts[case['status']]+=1
            links=[]
            for taxon in case['taxa']:
                ctx=contexts[case['marker'],taxon];pid=identifier(ctx);coord=primary[pid]
                links.append(dict(taxon=taxon,model_pair_id=pid,coordinate_status=coord['status'],
                    retained_residues=int(coord['retained_residues']),**{k:float(coord[k]) if coord[k] else None for k in METRICS}))
                matched_links+=1
            record=dict(marker=case['marker'],view_index=case['view_index'],cohort=case['cohort'],view=case['view'],
                case_status=case['status'],input_id=case['input_id'],taxa=case['taxa'],
                original_internal_branches=case['original_internal_branches'],branch_status_counts=expected_counts,
                coordinate_mask=dict(plddt=70,pae=10,scope='Whole-protein confidence mask; not marker-fit columns'),
                matched_taxon_coordinate_links=links,mapped_internal_paths=len(mapping),scientific_eligibility=False)
            cf.write(json.dumps(record,separators=(',',':'),allow_nan=False)+'\n')
            for group in mapping:
                split=group['split_mask_hex'];mask=int(split,16);original=group['original_branch_indices']
                assert original and aa[case['input_id']][split]['internal']
                sides=[[link for link in links if bool(mask&(1<<axes['positions'][link['taxon']]))==side] for side in [True,False]]
                assert min(map(len,sides))>=2
                native=[]
                for model in ['af','af_empirical','llm']:
                    for mode in ['iid_sites','circular_blocks10']:
                        source=paired[case['input_id'],split,model,mode];assert source['internal']=='True'
                        assert source['scientific_eligibility']=='False' and source['calibrated_effect_interval']=='False'
                        native.append({name:source[name] for name in ['model','mode','resampling_group','unit',
                            'point_AlphaFold','point_ESMFold','point_ESMFold_minus_AlphaFold','delta_q025','delta_median','delta_q975',
                            'positive_draws','negative_draws','zero_draws','calibrated_effect_interval','scientific_eligibility']})
                units={view['branches'][i]['branch_length_unit'] for i in original};assert len(units)==1
                out=dict(marker=case['marker'],view_index=case['view_index'],cohort=case['cohort'],view=case['view'],
                    input_id=case['input_id'],observed_split_mask_hex=split,original_branch_indices=original,
                    original_branch_count=len(original),projection='unique_internal_projection' if len(original)==1 else 'shared_internal_projection',
                    original_path_length_sum=math.fsum(view['branches'][i]['branch_length'] for i in original),
                    original_path_length_unit=next(iter(units)),aa_point_branch_length=aa[case['input_id']][split]['length'],
                    aa_point_unit='LG+F+G4 expected amino-acid substitutions/site',
                    side_taxa=[link['taxon'] for link in sides[0]],complement_taxa=[link['taxon'] for link in sides[1]],
                    side_coordinate_compared=sum(link['coordinate_status']=='coordinates_compared' for link in sides[0]),
                    complement_coordinate_compared=sum(link['coordinate_status']=='coordinates_compared' for link in sides[1]),
                    side_model_pair_ids=[link['model_pair_id'] for link in sides[0]],
                    complement_model_pair_ids=[link['model_pair_id'] for link in sides[1]],
                    paired_predictor_distributions=native,branch_localization_accepted=False,scientific_eligibility=False)
                pf.write(json.dumps(out,separators=(',',':'),allow_nan=False)+'\n');path_count+=1
    assert len(seen)==8750 and path_count==31290 and branch_slots==4523750
    assert dict(counts)==complete['case_status_counts'] and dict(branch_counts)==complete['branch_status_counts']
    for q in [case_path,path_path]:bind(pins,q)
    verify(pins)
    result=dict(status='complete_full_predictor_coordinate_phylogeny_links_pending_independent_readback',checked_utc=datetime.now(timezone.utc).isoformat(),
        marker_slots=125,tree_views=70,project_taxon_entries=526,comparison_cases=8750,mapped_internal_paths=31290,
        original_internal_branch_slots=4523750,case_status_counts=dict(counts),branch_status_counts=dict(branch_counts),
        matched_taxon_coordinate_link_occurrences=matched_links,paired_distribution_occurrences=path_count*6,
        source_hashes=pins,artifacts={q.name:sha(q) for q in [case_path,path_path]},
        scientific_eligibility=False,new_native_fits=0,gpu=False,all_eight_aims_incomplete=True,
        scope='Full125marker/70candidate-tree-view/526-taxon source universe retained. Every ready/nonready disposition and '
              'all4,523,750original branch slots audited;31,290observed internal paths link six conditional predictor '
              'resampling summaries, fixed-topology AA point estimates and both-side whole-protein coordinate controls. '
              'Many-to-one original paths remain explicit; no single-branch localization, physical rate, branch-time conversion, '
              'independence, causality, calibrated discovery or posterior adequacy accepted. Only21selected fungal taxa '
              'have predictor overlap; outgroups/fullfungal coverage remain upstream requirements. No cases filtered by geometry.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes']},indent=2))


if __name__=='__main__':main()
