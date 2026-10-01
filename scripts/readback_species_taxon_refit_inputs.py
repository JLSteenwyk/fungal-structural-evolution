#!/usr/bin/env python3
"""Independent complete native-refit matrix/membership/site-map reconstruction."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def fasta(path):
    records={};label=None;parts=[]
    with Path(path).open() as handle:
        for line in handle:
            text=line.strip()
            if not text:continue
            if text.startswith('>'):
                if label is not None:assert label not in records;records[label]=''.join(parts)
                label=text[1:].split()[0];parts=[]
            else:assert label is not None;parts.append(text)
    if label is not None:assert label not in records;records[label]=''.join(parts)
    return records


def rows(path):
    with Path(path).open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def run(plan_path,output):
    plan_path=Path(plan_path);plan=json.loads(plan_path.read_text());root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text());bindings=dict(plan['pins']);bind(bindings,plan_path);bind(bindings,rp)
    assert r['status']=='complete_species_taxon_refit_inputs_pending_independent_readback' and r['plan_sha256']==sha(plan_path) and r['scientific_eligibility'] is False
    for path,digest in r['source_hashes'].items():bind(bindings,path,digest)
    for name,digest in r['artifacts'].items():bind(bindings,root/name,digest)
    manifest=rows(plan['manifest']);taxa={v['taxon_id']:v for v in manifest};assert len(taxa)==len(manifest)==526
    base={label:fasta(Path(folder)/'matrix.faa') for label,folder in plan['matrices'].items()};assert all(set(v)==set(taxa) for v in base.values())
    columns={label:len(next(iter(v.values()))) for label,v in base.items()};assert all(len(seq)==columns[label] for label,index in base.items() for seq in index.values())
    coverage=rows(plan['coverage']);qc={v['taxon_id']:v for v in coverage};assert len(qc)==len(coverage)==526 and set(qc)==set(taxa)
    for label,index in base.items():
        for tid,sequence in index.items():
            observed=sum(sequence.count(a) for a in 'ACDEFGHIKLMNPQRSTVWY');assert observed==int(qc[tid][label+'_unambiguous_residues']) and abs(observed/columns[label]-float(qc[tid][label+'_fraction']))<1e-14
    hybrids={v['taxon_id'] for v in json.loads(Path(plan['hybrids']).read_text())['records']};label_rows=rows(plan['labels']);uncertain={v['taxon_id'] for v in label_rows if 'incompletely_identified_species_label' in v['flags'].split(';')}
    policies=['exclude_sparse_boundary_taxon','exclude_below10_both_alignments','exclude_curated_hybrids','exclude_hybrids_and_uncertain_labels'];assert r['policies']==policies and len(hybrids)==2 and len(uncertain)==21
    # Reconstruct the10% rule from full raw character strings; avoid rounded fraction thresholds.
    eligible={tid for tid in taxa if all(10*sum(base[label][tid].count(a) for a in 'ACDEFGHIKLMNPQRSTVWY')>=columns[label] for label in base)}
    kept={policies[0]:set(taxa)-{plan['sparse_boundary_taxon']},policies[1]:eligible,policies[2]:set(taxa)-hybrids,policies[3]:set(taxa)-hybrids-uncertain}
    membership={}
    for row in rows(root/'taxon_membership.tsv'):
        key=row['policy'],row['taxon_id'];assert key not in membership;membership[key]=row
    assert set(membership)=={(p,tid) for p in policies for tid in taxa}
    for (p,tid),row in membership.items():
        t=taxa[tid];wanted=dict(policy=p,taxon_id=tid,species_name=t['species_name'],study_role=t['study_role'],lineage=t['lineage'],retained=int(tid in kept[p]),curated_hybrid=int(tid in hybrids),uncertain_label=int(tid in uncertain),profile_canonical_fraction=qc[tid]['profile_fraction'],mafft_canonical_fraction=qc[tid]['mafft_fraction']);assert row=={k:str(v) for k,v in wanted.items()}
    retention={}
    for row in rows(root/'lineage_retention.tsv'):
        key=row['policy'],row['study_role'],row['lineage'];assert key not in retention;retention[key]=row
    groups={(v['study_role'],v['lineage']) for v in manifest};assert set(retention)=={(p,role,lineage) for p in policies for role,lineage in groups}
    for (p,role,lineage),row in retention.items():
        ids={tid for tid,t in taxa.items() if (t['study_role'],t['lineage'])==(role,lineage)};remaining=len(ids&kept[p]);wanted=dict(policy=p,study_role=role,lineage=lineage,original_taxa=len(ids),retained_taxa=remaining,excluded_taxa=len(ids)-remaining,entire_group_lost=int(remaining==0));assert row=={k:str(v) for k,v in wanted.items()}
    seen=set();retained_cells=0;checked_sequence_cells=0;roles={}
    for item in r['matrices']:
        label,policy=item['alignment'],item['policy'];assert (label,policy) not in seen;seen.add((label,policy));folder=Path(item['path']);assert folder==root/(label+'-'+policy)
        mrp=folder/'receipt.json';bind(bindings,mrp,item['matrix_receipt_sha256']);mr=json.loads(mrp.read_text());assert mr['status']=='prepared_native_species_taxon_sensitivity_matrix' and mr['scientific_eligibility'] is False
        observed=fasta(folder/'matrix.faa');expected={tid:base[label][tid] for tid in kept[policy]};assert observed==expected
        mrows=rows(folder/'taxa.tsv');assert len(mrows)==len({v['taxon_id'] for v in mrows})==len(expected) and {v['taxon_id']:v for v in mrows}=={tid:taxa[tid] for tid in expected}
        counts=dict(Counter(taxa[tid]['study_role'] for tid in expected));roles[policy]=counts
        assert len(expected)==item['taxa']==mr['taxa']==plan['expected_retention'][policy]['taxa'] and item['columns']==mr['columns']==columns[label]
        assert counts==item['roles']==mr['roles']==plan['expected_retention'][policy]['roles'] and mr['markers']==125
        source=Path(plan['matrices'][label]);assert mr['baseline_matrix_receipt_sha256']==sha(source/'receipt.json') and mr['baseline_matrix_sha256']==sha(source/'matrix.faa') and mr['manifest_sha256']==sha(plan['manifest'])
        for name,digest in mr['artifacts'].items():bind(bindings,folder/name,digest)
        assert item['matrix_sha256']==sha(folder/'matrix.faa')
        for name in ['partitions.nex','site_mapping.tsv']:assert (folder/name).read_bytes()==(source/name).read_bytes()
        retained_cells+=len(expected);checked_sequence_cells+=len(expected)*columns[label]
    assert seen=={(label,p) for label in base for p in policies} and len(seen)==8 and retained_cells==r['retained_taxon_matrix_cells']
    assert len(membership)==r['membership_rows']==2104 and len(retention)==r['lineage_retention_rows'];verify(bindings)
    result=dict(status='passed_complete_species_taxon_refit_input_readback',plan_sha256=sha(plan_path),producer_receipt_sha256=sha(rp),source_taxa=526,policies=policies,matrices_checked=len(seen),membership_rows_checked=len(membership),lineage_retention_rows_checked=len(retention),retained_taxon_matrix_cells=retained_cells,retained_sequence_character_cells_checked=checked_sequence_cells,policy_roles=roles,entire_lineage_groups_lost=sum(int(row['entire_group_lost']) for row in retention.values()),source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2));return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();run(args.plan,args.output)
