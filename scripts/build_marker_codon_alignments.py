#!/usr/bin/env python3
"""Project verified marker codons onto unchanged MAFFT protein columns."""
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from audit_busco_gene_copies import ROOT,sha,read_table


def project_codons(dna,aligned,code,columns):
    if len(dna)%3:raise ValueError('non_triplet_cds')
    translated=str(Seq(dna).translate(table=code));stop=translated.endswith('*')
    if stop:translated=translated[:-1];dna=dna[:-3]
    if '*' in translated:raise ValueError('internal_stop')
    if translated!=aligned.replace('-',''):raise ValueError('full_protein_translation_mismatch')
    if columns!=sorted(set(columns)) or any(i<1 or i>len(aligned) for i in columns):raise ValueError('invalid_column_mask')
    chosen=set(columns);pieces=[];position=0
    for col,aa in enumerate(aligned,1):
        codon='---' if aa=='-' else dna[3*position:3*position+3]
        position+=aa!='-'
        if col in chosen:pieces.append(codon)
    return ''.join(pieces),stop


def checked(folder,status):
    r=json.loads((folder/'receipt.json').read_text())
    if r['status']!=status:raise ValueError('Incomplete prerequisite '+str(folder))
    for name,h in r.get('artifacts',{}).items():
        if sha(folder/name)!=h:raise ValueError('Changed prerequisite artifact')
    return r


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable codon alignment')
    index=ROOT/'results/cds/marker-source-index-v1';ix=checked(index,'complete_full_marker_cds_identity_index')
    bounds=ROOT/'results/cds/ncbi-marker-boundaries-v1';br=checked(bounds,'complete_full_ncbi_marker_boundary_code_audit')
    if br['source_index_receipt_sha256']!=sha(index/'receipt.json'):raise ValueError('Boundary/index mismatch')
    projection=ROOT/'results/cds/published-outgroup-codon-projection-v2';pr=checked(projection,'complete_published_outgroup_genome_codon_projection_audit')
    if ix['projection_receipt_sha256']!=sha(projection/'receipt.json'):raise ValueError('External projection/index mismatch')
    audit=ROOT/'results/phylogeny/mafft-audit-full-v1';ar=json.loads((audit/'receipt.json').read_text())
    matrix=ROOT/'results/phylogeny/mafft-matrix-50-v1';mr=json.loads((matrix/'receipt.json').read_text())
    if ar['status']!='complete_audit' or mr['audit_receipt_sha256']!=sha(audit/'receipt.json') or mr['mask_sha256']!=sha(audit/'column_masks.json'):raise ValueError('Changed protein alignment masks')
    masks=json.loads((audit/'column_masks.json').read_text());alignments={r['marker']:r['alignment_sha256'] for r in ar['alignment_receipts']}
    links=read_table(index/'marker_cds_source_index.tsv');by_marker=defaultdict(list)
    for r in links:by_marker[r['marker']].append(r)
    if set(by_marker)!=set(alignments) or len(by_marker)!=125:raise ValueError('Full marker universe differs')
    boundary={(r['marker'],r['taxon_id']):r for r in read_table(bounds/'marker_annotation_translation.tsv')}
    expected_ncbi={(r['marker'],r['taxon_id']) for r in links if r['cds_source_kind']=='ncbi_unmodified_cds'}
    if set(boundary)!=expected_ncbi:raise ValueError('NCBI boundary coverage differs')
    wanted={(r['taxon_id'],r['protein_id']) for r in links};external={}
    with (projection/'codon_projection_audit.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if (r['taxon_id'],r['protein_id']) in wanted:external[r['taxon_id'],r['protein_id']]=r
    dna={r.id:str(r.seq) for r in SeqIO.parse(index/'marker_source_cds.fna','fasta')}
    a.output.mkdir(parents=True);summaries=[];rows=[];source_columns={};artifacts={}
    for marker,local in by_marker.items():
        path=ROOT/'results/phylogeny/alignments-full-v1'/(marker+'.faa')
        if sha(path)!=alignments[marker]:raise ValueError('Changed MAFFT alignment')
        records=list(SeqIO.parse(path,'fasta'));aas={r.id:str(r.seq).upper() for r in records}
        if len(aas)!=len(records) or set(aas)!={r['taxon_id'] for r in local}:raise ValueError('Protein alignment identity mismatch')
        columns=masks[marker]['0.5'];source_columns[marker]=columns;accepted={};localrows=[]
        for link in local:
            t=link['taxon_id'];pid=link['protein_id'];aligned=aas[t]
            if hashlib.sha256(aligned.replace('-','').encode()).hexdigest()!=link['sequence_sha256']:raise ValueError('Marker protein hash mismatch')
            row={'marker':marker,'taxon_id':t,'protein_id':pid,'status':'','translation_table':'','code_provenance':'','source_kind':link['cds_source_kind'],'annotation_flags':'','initial_phase_or_omitted_bases':'','terminal_partial_bases':'','terminal_stop_removed':'','noncanonical_codon_columns':'','representative_decision':link['representative_decision'],'gene_mapping_status':link['gene_mapping_status']}
            reasons=[];code=1
            if link['cds_source_kind']=='ncbi_unmodified_cds':
                b=boundary[marker,t]
                if b['protein_id']!=pid:raise ValueError('Boundary protein identity mismatch')
                row.update(code_provenance=b.get('code_provenance',''),annotation_flags=b.get('annotation_flags',''),initial_phase_or_omitted_bases=b.get('initial_phase',''))
                if b['translation_status']!='exact_translation':reasons.append(b['translation_status'])
                if b.get('annotation_issues'):reasons.append(b['annotation_issues'])
                if b['annotation_status']!='consistent_ordered_cds_features':reasons.append(b['annotation_status'])
                if b.get('initial_phase')!='0':reasons.append('initial_phase_requires_coordinate_review')
                if set(b.get('annotation_flags','').split(';'))&{'exception','transl_except','pseudo','pseudogene'}:reasons.append('annotation_exception_requires_review')
                if b.get('translation_table'):code=int(b['translation_table'])
                else:reasons.append('unresolved_translation_code')
            elif link['cds_source_kind']=='published_genome_projected_codon_span':
                b=external[t,pid]
                if b['projection_status']!='exact_genome_linked_codon_translation':reasons.append('unverified_external_projection')
                row.update(code_provenance='external_table_1_genome_projection',initial_phase_or_omitted_bases=b['initial_partial_bases'],terminal_partial_bases=b['terminal_partial_bases'],annotation_flags='partial_boundaries' if int(b['initial_partial_bases']) or int(b['terminal_partial_bases']) else '')
            else:row.update(code_provenance='external_table_1_verified_extraction',annotation_flags='external_extraction_policy_requires_downstream_review')
            row['translation_table']=code
            if link['source_status']!='unique_source_cds':reasons.append(link['source_status'])
            if reasons:row['status']='excluded:'+ '|'.join(sorted(set(reasons)))
            else:
                try:
                    aligned_dna,stop=project_codons(dna[marker+'|'+t],aligned,code,columns)
                    accepted[t]=aligned_dna;row.update(status='translation_verified_codon_alignment',terminal_stop_removed=stop,noncanonical_codon_columns=sum(c!='---' and not set(c)<=set('ACGT') for c in [aligned_dna[i:i+3] for i in range(0,len(aligned_dna),3)]))
                except ValueError as e:row['status']='excluded:'+str(e)
            localrows.append(row);rows.append(row)
        target=a.output/(marker+'.fna')
        with target.open('w') as f:
            for t,seq in sorted(accepted.items()):f.write('>'+t+'\n'+seq+'\n')
        artifacts[target.name]=sha(target)
        summaries.append({'marker':marker,'input_taxa':len(local),'accepted_taxa':len(accepted),'codon_columns':len(columns),'nucleotide_columns':len(columns)*3,'translation_codes':dict(Counter(str(r['translation_table']) for r in localrows if r['status']=='translation_verified_codon_alignment')),'status_counts':dict(Counter(r['status'] for r in localrows))})
        print(marker,len(summaries),len(accepted),flush=True)
    for name,data in [('sequence_audit.tsv',rows)]:
        target=a.output/name
        with target.open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
        artifacts[name]=sha(target)
    target=a.output/'protein_source_columns.json';target.write_text(json.dumps(source_columns)+'\n');artifacts[target.name]=sha(target)
    result={'status':'complete_full_marker_codon_alignment_projection','markers':summaries,'marker_links':len(rows),'status_counts':dict(Counter(r['status'] for r in rows)),'source_index_receipt_sha256':sha(index/'receipt.json'),'boundary_receipt_sha256':sha(bounds/'receipt.json'),'external_projection_receipt_sha256':sha(projection/'receipt.json'),'protein_alignment_audit_receipt_sha256':sha(audit/'receipt.json'),'protein_mask_sha256':sha(audit/'column_masks.json'),'script_sha256':sha(Path(__file__)),'artifacts':artifacts,'interpretation':'Diagnostic codon alignments under exact full protein correspondence, using unchanged MAFFT 50-percent protein masks. Terminal stop codon removed only after verified translation. Partial CDSs and unresolved/alternative gene products remain flagged. Codes may differ within a marker: do not run a single-code selection model across mixed-code taxa. Alignment robustness, clade-specific divergence/saturation, copy/gene policies and selection model eligibility remain pending.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
