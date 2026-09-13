#!/usr/bin/env python3
"""Classify deposited entity sequences against exact accession-linked model sequences."""
import argparse,json,hashlib
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha
from prepare_paired_phylogenetic_inputs import write_table
AA=set('ACDEFGHIKLMNPQRSTVWY')


def classify(entity,target):
    if not entity or not target:return 'missing_sequence',[]
    if not set(entity)<=AA or not set(target)<=AA:return 'noncanonical_sequence_requires_review',[]
    if entity==target:return 'exact_full_sequence',[1]
    if len(entity)<len(target):
        starts=[i+1 for i in range(len(target)-len(entity)+1) if target.startswith(entity,i)]
        if starts:return ('exact_fragment_unique' if len(starts)==1 else 'exact_fragment_ambiguous'),starts
    if target in entity:return 'target_contained_in_longer_construct',[]
    return 'requires_alignment_or_variant_review',[]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--metadata',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--allow-partial',action='store_true');a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use an immutable screen output')
    config=json.loads((a.metadata/'config.json').read_text());complete=a.metadata/'receipt.json'
    if not a.allow_partial and not complete.exists():raise ValueError('Full metadata retrieval required unless partial explicitly requested')
    mapping=ROOT/'results/structural_markers/gdm-expanded-v1';checked_receipt(mapping)
    inputs=ROOT/'data/domains/marker-inputs-v1';checked_receipt(inputs)
    models=json.loads((mapping/'model_provenance.json').read_text());by_accession=defaultdict(list)
    for m in models:by_accession[m['uniprot_accession']].append(m)
    sequences={hashlib.sha256(str(r.seq).encode()).hexdigest():str(r.seq) for r in SeqIO.parse(inputs/'sequences.faa','fasta')}
    for m in models:
        if m['sequence_sha256'] not in sequences or len(sequences[m['sequence_sha256']])!=m['length']:raise ValueError('Model sequence identity differs')
    inventory=ROOT/'data/experimental_structures/accession-inventory-v1';checked_receipt(inventory)
    if config['inventory_receipt_sha256']!=sha(inventory/'receipt.json'):raise ValueError('Inventory lineage differs')
    expected=(inventory/'polymer_entity_ids.txt').read_text().splitlines();rows=[];unlinked=[];pins={};done=0
    for identifier in expected:
        rp=a.metadata/'polymer_entity'/(identifier+'.receipt.json');ep=a.metadata/'polymer_entity'/(identifier+'.json')
        if not rp.exists():
            if a.allow_partial:continue
            raise ValueError('Missing entity response')
        r=json.loads(rp.read_text())
        if r['config_sha256']!=sha(a.metadata/'config.json') or r['response_sha256']!=sha(ep):raise ValueError('Changed entity response')
        entity=json.loads(ep.read_text())
        if entity['rcsb_id']!=identifier:raise ValueError('Entity ID differs')
        pins[identifier]=sha(rp);done+=1
        poly=entity.get('entity_poly',{});seq=''.join(poly.get('pdbx_seq_one_letter_code_can','').split())
        refs=entity.get('rcsb_polymer_entity_container_identifiers',{}).get('reference_sequence_identifiers',[])
        accessions={r['database_accession'] for r in refs if r.get('database_name')=='UniProt'}&set(by_accession)
        if not accessions:unlinked.append({'entity_id':identifier,'reason':'No currently matching UniProt reference in retrieved entity metadata'});continue
        for accession in sorted(accessions):
            for m in by_accession[accession]:
                target=sequences[m['sequence_sha256']];status,starts=classify(seq,target)
                rows.append({'entity_id':identifier,'uniprot_accession':accession,'model_id':m['model_id'],'model_sequence_sha256':m['sequence_sha256'],'entity_sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'entity_length':len(seq),'target_length':len(target),'sequence_class':status,'target_fragment_starts_1based':';'.join(map(str,starts)),'polymer_type':poly.get('type','unknown'),'reported_mutation_count':poly.get('rcsb_mutation_count','unknown'),'reported_nonstandard_monomers':poly.get('rcsb_non_std_monomer_count','unknown'),'reported_artifact_monomers':poly.get('rcsb_artifact_monomer_count','unknown'),'benchmark_eligibility':'pending_coordinates_quality_context_and_training_overlap'})
    if not done:raise ValueError('No verified entity metadata available')
    a.output.mkdir(parents=True)
    if rows:write_table(a.output/'sequence_correspondence.tsv',rows)
    if unlinked:write_table(a.output/'unlinked_entities.tsv',unlinked)
    result={'status':'complete_frozen_partial_entity_sequence_screen' if done<len(expected) else 'complete_entity_sequence_screen','entities_screened':done,'entities_expected':len(expected),'entities_pending':len(expected)-done,'entity_model_rows':len(rows),'unlinked_entities':len(unlinked),'sequence_class_counts':dict(Counter(r['sequence_class'] for r in rows)),'mapping_receipt_sha256':sha(mapping/'receipt.json'),'sequence_input_receipt_sha256':sha(inputs/'receipt.json'),'metadata_config_sha256':sha(a.metadata/'config.json'),'entity_receipt_sha256':pins,'script_sha256':sha(Path(__file__)),'interpretation':'Accession-linked deposited sequence comparisons only. Unique fragments retain exact offsets; ambiguous fragments, noncanonical sequences and longer constructs remain explicit. Exact canonical sequence does not establish chemical identity, complete observed residues, experimental quality, matching biological state or independence from prediction training. Partial screens follow retrieval order and are not representative samples.','artifacts':{p.name:sha(p) for p in a.output.iterdir() if p.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='entity_receipt_sha256'},indent=2))


if __name__=='__main__':main()
