#!/usr/bin/env python3
"""Check every refreshed AlphaFold residue link against alignments and mmCIF."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import time
from Bio import AlignIO, SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import numpy as np
from assess_small_family_output_exposure import sha


def projected_positions(alignment, retained_columns, site_pairs):
    occupied=np.fromiter((c not in '.-' for c in alignment),dtype=bool)
    ordinal=np.cumsum(occupied)
    result=[]
    for profile_column,matrix_column in site_pairs:
        if not 1<=profile_column<=len(retained_columns):raise ValueError('Invalid profile column')
        position=retained_columns[profile_column-1]-1
        if not 0<=position<len(alignment):raise ValueError('Invalid Stockholm column')
        if occupied[position]:result.append((matrix_column,int(ordinal[position])))
    return result


def check_row(actual, identity, matrix_column, residue, confidence):
    if actual is None:raise ValueError('Missing residue row')
    expected=dict(identity,matrix_column_1based=str(matrix_column),protein_residue_1based=str(residue))
    if any(actual.get(k)!=v for k,v in expected.items()) or float(actual['ca_plddt'])!=confidence:
        raise ValueError('Residue identity, position or confidence differs')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',required=True,type=Path)
    a=p.parse_args();plan=json.loads(a.plan.read_text())
    source=json.loads(Path(plan['preparation_plan']).read_text())
    pins={str(a.plan):sha(a.plan),**plan['pins'],**source['pins']}
    def verify():
        for path,digest in pins.items():
            if sha(path)!=digest:raise ValueError('Changed pinned input: '+path)
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);start=time.time()
    def state(status,**details):
        temp=out/'state.tmp';temp.write_text(json.dumps(dict(status=status,elapsed_seconds=time.time()-start,**details))+'\n');temp.replace(out/'state.json')
    try:
        state('waiting_for_exact_mapping_controller')
        proc=Path('/proc')/str(plan['predecessor_pid'])/'stat'
        while proc.exists():
            try:fields=proc.read_text().rsplit(')',1)[1].split()
            except FileNotFoundError:break
            if fields[19]!=str(plan['predecessor_start_ticks']) or fields[0]=='Z':break
            time.sleep(20)
        controller=Path(source['controller_output'])/'mapping/receipt.json'
        completed=json.loads(controller.read_text())
        if completed['status']!='complete_refreshed_mapping_and_catalog_agreement_pending_residue_readback':
            raise ValueError('Mapping did not complete')
        mapping=Path(source['mapping_output']);mr=json.loads((mapping/'receipt.json').read_text())
        if completed['result_receipt_sha256']!=sha(mapping/'receipt.json') or completed['plan_sha256']!=sha(plan['preparation_plan']):
            raise ValueError('Mapping controller binding differs')
        pins[str(controller)]=sha(controller);pins[str(mapping/'receipt.json')]=sha(mapping/'receipt.json')
        pins.update({str(mapping/name):digest for name,digest in mr['artifacts'].items()})
        verify()
        if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
        records=list(csv.DictReader(Path('results/phylogeny/markers-full-v1/protein_mapping.tsv').open(),delimiter='\t'))
        by_marker={}
        for r in records:
            if r['taxon_id'] in by_marker.setdefault(r['marker'],{}):raise ValueError('Duplicate source marker/taxon')
            by_marker[r['marker']][r['taxon_id']]=r
        catalog=list(csv.DictReader((Path(source['catalog'])/'marker_model_links.tsv').open(),delimiter='\t'))
        selected={(r['marker'],r['taxon_id']):r for r in catalog}
        if len(selected)!=len(catalog):raise ValueError('Duplicate catalog link')
        links=list(csv.DictReader((mapping/'marker_structure_links.tsv').open(),delimiter='\t'))
        mapped={(r['marker'],r['taxon_id']):r for r in links}
        if len(mapped)!=len(links) or mapped.keys()!=selected.keys():raise ValueError('Mapping/catalog link universe differs')
        model_rows=json.loads((mapping/'model_provenance.json').read_text())
        models={r['path']:r for r in model_rows}
        if len(models)!=len(model_rows):raise ValueError('Repeated model path')
        site_map={}
        with Path('results/phylogeny/profile-matrix-50-v1/site_mapping.tsv').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                site_map.setdefault(row['marker'],[]).append((int(row['alignment_column_1based']),int(row['matrix_column_1based'])))
        cache={};n=0;checked_links=0;marker_counts=[]
        with gzip.open(mapping/'matrix_to_structure_residues.tsv.gz','rt') as handle:
            output=iter(csv.DictReader(handle,delimiter='\t'))
            for marker,taxa in sorted(by_marker.items()):
                raw=Path('results/phylogeny/markers-full-v1/unaligned')/(marker+'.faa')
                seqs={r.id:str(r.seq) for r in SeqIO.parse(raw,'fasta')}
                profile=Path('results/phylogeny/profile-alignments-full-v1')/marker
                pr=json.loads(profile.with_suffix('.receipt.json').read_text())
                aln={r.id:str(r.seq) for r in AlignIO.read(profile.with_suffix('.sto'),'stockholm')}
                marker_n=0
                for taxon,row in sorted(taxa.items()):
                    if (marker,taxon) not in selected:continue
                    choice=selected[marker,taxon];link=mapped[marker,taxon];model=models[choice['model_path']]
                    sequence=seqs[taxon]
                    if (hashlib.sha256(sequence.encode()).hexdigest()!=row['sequence_sha256']
                            or aln[taxon].upper().replace('-','').replace('.','')!=sequence):
                        raise ValueError('Source sequence/alignment mismatch')
                    if model['path'] not in cache:
                        if sha(model['path'])!=model['sha256']:raise ValueError('Changed coordinates')
                        cif=MMCIF2Dict(model['path'])
                        polymer=[''.join(x.split()) for x in cif['_entity_poly.pdbx_seq_one_letter_code_can']]
                        ca=[(int(pos),float(value)) for atom,pos,value in zip(cif['_atom_site.label_atom_id'],cif['_atom_site.label_seq_id'],cif['_atom_site.B_iso_or_equiv']) if atom=='CA']
                        positions=[x[0] for x in ca]
                        if polymer!=[sequence] or len(positions)!=len(sequence) or set(positions)!=set(range(1,len(sequence)+1)):
                            raise ValueError('Nonunique/incomplete coordinates or wrong polymer')
                        values=dict(ca)
                        confidence=np.array([values[i] for i in range(1,len(sequence)+1)])
                        if not np.isfinite(confidence).all() or (confidence<0).any() or (confidence>100).any():raise ValueError('Invalid pLDDT')
                        cache[model['path']]=(sequence,confidence)
                    saved,confidence=cache[model['path']]
                    if saved!=sequence:raise ValueError('Shared coordinate assigned to different sequence')
                    projected=projected_positions(aln[taxon],pr['retained_stockholm_columns_1based'],sorted(site_map[marker]))
                    identity=dict(marker=marker,taxon_id=taxon,protein_id=row['protein_id'],model_id=model['model_id'],model_version=str(model['version']))
                    observed=[]
                    for column,residue in projected:
                        value=float(confidence[residue-1]);check_row(next(output,None),identity,column,residue,value)
                        observed.append(value);n+=1;marker_n+=1
                    if len(projected)!=int(link['retained_marker_residues']):raise ValueError('Link residue count differs')
                    if observed:
                        if (not math.isclose(sum(observed)/len(observed),float(link['mean_retained_ca_plddt']),abs_tol=1e-10)
                                or not math.isclose(sum(v>=70 for v in observed)/len(observed),float(link['fraction_retained_ca_plddt_ge70']),abs_tol=1e-12)):
                            raise ValueError('Link confidence summary differs')
                    elif link['mean_retained_ca_plddt'] or link['fraction_retained_ca_plddt_ge70']:
                        raise ValueError('Unexpected empty-link confidence summary')
                    checked_links+=1
                marker_counts.append(dict(marker=marker,residue_links=marker_n))
                state('checking_full_residue_mapping',markers=len(marker_counts),residue_links=n,models=len(cache))
            if next(output,None) is not None:raise ValueError('Extra residue rows')
        if n!=mr['matrix_residue_links'] or checked_links!=mr['marker_proteins_linked'] or set(cache)!=set(models):
            raise ValueError('Incomplete mapped universe')
        verify()
        receipt=dict(status='passed_full_refreshed_afdb_residue_mapping_readback',models=len(cache),marker_links=checked_links,
            markers=len(marker_counts),matrix_residue_links=n,mapping_receipt_sha256=sha(mapping/'receipt.json'),
            plan_sha256=sha(a.plan),script_sha256=sha(__file__),elapsed_seconds=time.time()-start,
            scope='Every exported matrix/sequence/coordinate position and pLDDT checked, with complete model/link coverage and per-link confidence summaries. Cumulative non-gap projection is separate from producer code; BioPython Stockholm/mmCIF parsers and upstream alignments are shared. No PAE, coordinate accuracy or biological inference validation.')
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');state(receipt['status'])
    except Exception as exc:
        state('failed_requires_review',error=repr(exc));raise


if __name__=='__main__':main()
