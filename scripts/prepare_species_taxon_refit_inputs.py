#!/usr/bin/env python3
"""Prepare complete native taxon sensitivities without replacing the 526-taxon baseline."""
import argparse
import csv
import json
import shutil
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

POLICIES=['exclude_sparse_boundary_taxon','exclude_below10_both_alignments','exclude_curated_hybrids','exclude_hybrids_and_uncertain_labels']


def table(path,rows):
    with Path(path).open('w') as handle:
        writer=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def run(plan_path):
    plan_path=Path(plan_path);plan=json.loads(plan_path.read_text());bindings=dict(plan['pins']);bind(bindings,plan_path);verify(bindings)
    for path,status,journals in [(plan['baseline_completion'],'complete_verified_full_four_run_pmsf_ML_consensus_sensitivity',2),(plan['boundary_completion'],'complete_verified_full_pmsf_role_boundary_and_character_diagnostics',1)]:
        proof=json.loads(Path(path).read_text());assert proof['status']==status and len(proof['services'])==journals and proof['scientific_eligibility'] is False;bind(bindings,path)
        for p,digest in proof['source_hashes'].items():bind(bindings,p,digest)
    manifest=list(csv.DictReader(Path(plan['manifest']).open(),delimiter='\t'));taxa={r['taxon_id']:r for r in manifest};assert len(taxa)==len(manifest)==526
    assert Counter(r['study_role'] for r in manifest)=={'ingroup':501,'outgroup':25}
    coverage=list(csv.DictReader(Path(plan['coverage']).open(),delimiter='\t'));qc={r['taxon_id']:r for r in coverage};assert len(qc)==len(coverage)==526 and set(qc)==set(taxa)
    curation=json.loads(Path(plan['hybrids']).read_text());hybrids={r['taxon_id'] for r in curation['records']};assert len(hybrids)==2
    labels=list(csv.DictReader(Path(plan['labels']).open(),delimiter='\t'));uncertain={r['taxon_id'] for r in labels if 'incompletely_identified_species_label' in r['flags'].split(';')};assert len(labels)==22 and len(uncertain)==21
    excluded={'exclude_sparse_boundary_taxon':{plan['sparse_boundary_taxon']},'exclude_below10_both_alignments':{tid for tid,r in qc.items() if float(r['profile_fraction'])<.1 or float(r['mafft_fraction'])<.1},'exclude_curated_hybrids':hybrids,'exclude_hybrids_and_uncertain_labels':hybrids|uncertain}
    # The established policy requires at least10% canonical coverage in each alignment.
    assert list(excluded)==POLICIES and all(v<=set(taxa) for v in excluded.values())
    expected=plan['expected_retention'];assert all(len(taxa)-len(excluded[p])==expected[p]['taxa'] for p in POLICIES)
    groups=sorted({(r['study_role'],r['lineage']) for r in manifest});membership=[];retention=[]
    for policy in POLICIES:
        keep=set(taxa)-excluded[policy];roles=Counter(taxa[tid]['study_role'] for tid in keep);assert dict(roles)==expected[policy]['roles']
        for r in manifest:membership.append(dict(policy=policy,taxon_id=r['taxon_id'],species_name=r['species_name'],study_role=r['study_role'],lineage=r['lineage'],retained=int(r['taxon_id'] in keep),curated_hybrid=int(r['taxon_id'] in hybrids),uncertain_label=int(r['taxon_id'] in uncertain),profile_canonical_fraction=qc[r['taxon_id']]['profile_fraction'],mafft_canonical_fraction=qc[r['taxon_id']]['mafft_fraction']))
        for role,lineage in groups:
            ids={r['taxon_id'] for r in manifest if (r['study_role'],r['lineage'])==(role,lineage)};retention.append(dict(policy=policy,study_role=role,lineage=lineage,original_taxa=len(ids),retained_taxa=len(ids&keep),excluded_taxa=len(ids-keep),entire_group_lost=int(not ids&keep)))
    out=Path(plan['output']);assert shutil.disk_usage(out.parent).free>=plan['resources']['minimum_free_disk_gib']*2**30;out.mkdir(exist_ok=False)
    table(out/'taxon_membership.tsv',membership);table(out/'lineage_retention.tsv',retention);matrices=[]
    for label,folder in plan['matrices'].items():
        source=Path(folder);rp=source/'receipt.json';r=json.loads(rp.read_text());assert r['manifest_sha256']==sha(plan['manifest']) and r['taxa']==526 and r['markers']==125;bind(bindings,rp)
        for name,digest in r['artifacts'].items():bind(bindings,source/name,digest)
        records=list(SeqIO.parse(source/'matrix.faa','fasta'));assert len(records)==526 and {v.id for v in records}==set(taxa) and all(len(v.seq)==r['columns'] for v in records)
        for v in records:
            counts=Counter(str(v.seq));n=sum(counts[a] for a in 'ACDEFGHIKLMNPQRSTVWY');assert n==int(qc[v.id][label+'_unambiguous_residues']) and abs(n/r['columns']-float(qc[v.id][label+'_fraction']))<1e-14
        for policy in POLICIES:
            folder=out/(label+'-'+policy);folder.mkdir();kept=[v for v in records if v.id not in excluded[policy]];SeqIO.write(kept,folder/'matrix.faa','fasta');table(folder/'taxa.tsv',[taxa[v.id] for v in kept])
            for name in ['partitions.nex','site_mapping.tsv']:shutil.copyfile(source/name,folder/name)
            assert {v.id:str(v.seq) for v in SeqIO.parse(folder/'matrix.faa','fasta')}=={v.id:str(v.seq) for v in kept}
            mr=dict(status='prepared_native_species_taxon_sensitivity_matrix',alignment=label,policy=policy,taxa=len(kept),columns=r['columns'],markers=r['markers'],roles=expected[policy]['roles'],baseline_matrix_receipt_sha256=sha(rp),baseline_matrix_sha256=sha(source/'matrix.faa'),manifest_sha256=sha(plan['manifest']),artifacts={p.name:sha(p) for p in folder.iterdir()},scientific_eligibility=False)
            with (folder/'receipt.json').open('x') as handle:handle.write(json.dumps(mr,indent=2)+'\n')
            matrices.append(dict(alignment=label,policy=policy,path=str(folder),matrix_receipt_sha256=sha(folder/'receipt.json'),matrix_sha256=sha(folder/'matrix.faa'),taxa=len(kept),columns=r['columns'],roles=expected[policy]['roles']))
    verify(bindings)
    result=dict(status='complete_species_taxon_refit_inputs_pending_independent_readback',plan_sha256=sha(plan_path),source_taxa=526,policies=POLICIES,matrices=matrices,membership_rows=len(membership),lineage_retention_rows=len(retention),retained_taxon_matrix_cells=sum(r['taxa'] for r in matrices),source_hashes=bindings,artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},scientific_eligibility=False,scope=plan['scope'])
    with (out/'receipt.json').open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2));return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);run(parser.parse_args().plan)
