#!/usr/bin/env python3
"""Map observed conflict-column residues to sequence-identical extant models."""
import csv,hashlib,json
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils import seq1


def main():
    pins={}
    def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    def checked(root,name):
        root=Path(root);rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name
        assert sha(p)==r['artifacts'][name];pins[str(rp)]=sha(rp);pins[str(p)]=sha(p);return p
    base=Path('results/structural_comparisons')
    contexts=json.loads(checked('results/ancestral/ancestral-conflict-observations-20260927-v1','observed_contexts.json').read_text())
    families={r['family'] for r in contexts}
    with checked(base/'whole-domain-case-dossiers-20260927-v1','whole_reference_links.tsv').open() as h:
        links={r['family']:r for r in csv.DictReader(h,delimiter='\t') if r['guide']=='profile' and r['family'] in families}
    with checked(base/'whole-protein-common-residues-20260927-v1','model_triads.tsv').open() as h:
        triads={r['triad_id']:r for r in csv.DictReader(h,delimiter='\t') if r['triad_id'] in {x['triad_id'] for x in links.values()}}
    roles={};wanted={}
    for family,link in links.items():
        triad=triads[link['triad_id']]
        for role,gene_field in [('a','gene_a'),('b','gene_b'),('reference','reference_gene')]:
            key=(triad[role+'_model'],int(triad[role+'_version']))
            roles[family,role]=(link[gene_field],key);wanted[key]=triad[role+'_sequence_sha256']
    found={}
    for folder in ['duplication-alignment-inputs-20260926-v1','duplication-reference-alignment-inputs-20260926-v1']:
        path=checked(base/folder,'inputs.jsonl')
        with path.open() as h:
            for line in h:
                r=json.loads(line);key=(r['model_id'],r['version'])
                if key not in wanted or r['mask']!='full':continue
                assert r['status']=='ready' and hashlib.sha256(r['sequence'].encode()).hexdigest()==wanted[key]
                assert r['original_positions']==list(range(1,r['original_length']+1))
                assert sha(r['path'])==r['sha256'];pins[r['path']]=r['sha256'];found[key]=r
    assert set(found)==set(wanted) and len(found)==6
    ca={}
    for key,r in found.items():
        atoms={}
        for line in Path(r['path']).read_text().splitlines():
            if line.startswith('ATOM') and line[12:16].strip()=='CA':
                pos=int(line[22:26]);assert pos not in atoms
                atoms[pos]=(seq1(line[17:20]),float(line[60:66]),[float(line[i:i+8]) for i in [30,38,46]])
        assert ''.join(atoms[i][0] for i in range(1,len(atoms)+1))==r['sequence']
        ca[key]=atoms
    rows=[];seen=set();sequences={}
    for c in contexts:
        identity=(c['family'],c['coordinate_signature_sha256'],c['boundary'],c['whole_method'],c['whole_column'])
        if identity in seen:continue
        seen.add(identity)
        family=c['family'];observations={r['gene']:r for r in c['observations']}
        if family not in sequences:
            p=Path('results/ancestral/case-alignments-20260927-v1')/(family+'-mafft')/'alignment.faa'
            sequences[family]={r.id:str(r.seq).replace('-','') for r in SeqIO.parse(p,'fasta')};pins[str(p)]=sha(p)
        for role in ['a','b','reference']:
            gene,key=roles[family,role];model=found[key];assert sequences[family][gene]==model['sequence']
            obs=observations[gene];position=obs['protein_position']
            row=dict(family=family,coordinate_signature_sha256=c['coordinate_signature_sha256'],boundary=c['boundary'],whole_method=c['whole_method'],whole_column=c['whole_column'],role=role,gene=gene,model_id=key[0],version=key[1],sequence_sha256=wanted[key],protein_position=position,residue=obs['residue'],domain_coordinate_category=obs['category'],status='alignment_gap_no_residue' if position is None else 'exact_sequence_CA_mapped',plddt=None,ca_xyz=None)
            if position is not None:
                aa,confidence,xyz=ca[key][position];assert aa==obs['residue'];row.update(plddt=confidence,ca_xyz=xyz)
            rows.append(row)
    output=Path('results/ancestral/ancestral-conflict-case-models-20260927-v1');output.mkdir(exist_ok=False)
    p=output/'extant_model_annotations.json';p.write_text(json.dumps(rows,indent=2)+'\n')
    result=dict(status='complete_extant_case_model_conflict_annotations',models=len(found),annotation_rows=len(rows),mapped_residue_rows=sum(r['status']=='exact_sequence_CA_mapped' for r in rows),retained_domain_rows=sum(r['domain_coordinate_category']=='retained_domain_coordinate' for r in rows),pins=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p)},scope='Exact full-sequence equality and CA identity checked for six existing extant AlphaFold models. pLDDT is local prediction confidence, not validation of homology, ancestral state or function. Outside-domain residues remain flagged and must not be represented as matching focal-domain positions. CA-only coordinates do not establish surface, pocket or catalytic roles. No ancestral structures predicted.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');result.update(completed_receipt_path=str(output/'receipt.json'),completed_receipt_sha256=sha(output/'receipt.json'))
    Path('metadata/ancestral_conflict_case_models_completed_20260927.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['pins','artifacts']}))


if __name__=='__main__':main()
