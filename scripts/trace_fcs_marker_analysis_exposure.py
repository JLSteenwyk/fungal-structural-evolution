#!/usr/bin/env python3
"""Trace exact flagged marker proteins into frozen structural and codon inputs."""
import argparse,csv,json,hashlib
from collections import Counter
from pathlib import Path
from Bio import SeqIO


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['mapping','audit','fits','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable output')
    receipt=json.loads((a.mapping/'receipt.json').read_text());audit=json.loads((a.audit/'receipt.json').read_text())
    if audit['mapping_receipt_sha256']!=sha(a.mapping/'receipt.json') or audit['status']!='passed_full_region_cds_intersection_and_protein_marker_join_audit':raise ValueError('Missing full overlap audit')
    path=a.mapping/'marker_overlap_review.tsv'
    if sha(path)!=receipt['artifacts'][path.name]:raise ValueError('Changed overlap table')
    flags={(r['marker'],r['taxon_id']):r for r in table(path) if r['fcs_overlap_status']=='cds_overlaps_fcs_region'}
    source_specs=[('esmfold_earlier','results/structural_markers/esmfold-partial-v1','results/phylogeny/paired-inputs-esmfold-partial-v1'),('esmfold_combined','results/structural_markers/esmfold-combined-frozen-v1','results/phylogeny/paired-inputs-esmfold-combined-v1'),('alphafold_expanded','results/structural_markers/gdm-expanded-v1','results/phylogeny/paired-inputs-gdm-expanded-v1')]
    rows=[];sources={}
    for label,snapshot,inputs in source_specs:
        snapshot,inputs=Path(snapshot),Path(inputs);sr=json.loads((snapshot/'receipt.json').read_text());ir=json.loads((inputs/'receipt.json').read_text())
        lp=snapshot/'marker_structure_links.tsv'
        if sha(lp)!=sr['artifacts'][lp.name]:raise ValueError('Changed structural links')
        links={(r['marker'],r['taxon_id']):r for r in table(lp)}
        sources[label]={'snapshot_receipt_sha256':sha(snapshot/'receipt.json'),'input_receipt_sha256':sha(inputs/'receipt.json')}
        taxa_cache={}
        for key,flag in sorted(flags.items()):
            marker,taxon=key;link=links.get(key);inpaired=False
            if link and (link['protein_id']!=flag['protein_id'] or link['sequence_sha256']!=flag['sequence_sha256']):raise ValueError('Structural protein identity differs')
            if marker not in taxa_cache:
                name=marker+'/aa.faa';f=inputs/name
                if name in ir['artifacts']:
                    if sha(f)!=ir['artifacts'][name]:raise ValueError('Changed paired input')
                    taxa_cache[marker]={r.id for r in SeqIO.parse(f,'fasta')}
                else:taxa_cache[marker]=set()
            inpaired=taxon in taxa_cache[marker]
            if inpaired and not link:raise ValueError('Paired protein lacks structural mapping')
            rows.append({'dataset':label,'marker':marker,'taxon_id':taxon,'protein_id':flag['protein_id'],'fcs_actions':flag['fcs_actions'],'exact_structure_link_present':bool(link),'present_in_ready_paired_alignment':inpaired,'model_id':link['model_id'] if link else ''})
    fits=json.loads((a.fits/'receipt.json').read_text());codon=[]
    for cr in fits['case_receipts']:
        case=cr['case_id'];marker=case.split('__')[1];wanted={t:r for (m,t),r in flags.items() if m==marker}
        if not wanted:continue
        root=a.fits/case;rp=root/'receipt.json';r=json.loads(rp.read_text());cp=root/'config.json';c=json.loads(cp.read_text())
        if sha(rp)!=cr['receipt_sha256'] or sha(cp)!=r['config_sha256']:raise ValueError('Changed codon provenance')
        cmd=c['command'];alignment=Path(cmd[cmd.index('--alignment')+1])
        if sha(alignment)!=c['alignment_sha256']:raise ValueError('Changed codon input')
        taxa={r.id for r in SeqIO.parse(alignment,'fasta')}
        for taxon in sorted(taxa & set(wanted)):
            flag=wanted[taxon];codon.append({'case_id':case,'marker':marker,'taxon_id':taxon,'protein_id':flag['protein_id'],'fcs_actions':flag['fcs_actions'],'alignment_sha256':sha(alignment),'fit_receipt_sha256':sha(rp)})
    a.output.mkdir(parents=True)
    for name,data in [('structural_dataset_exposure.tsv',rows),('codon_case_exposure.tsv',codon)]:
        with (a.output/name).open('w',newline='') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    r={'status':'complete_exact_flagged_marker_input_exposure_trace','mapping_receipt_sha256':sha(a.mapping/'receipt.json'),'audit_receipt_sha256':sha(a.audit/'receipt.json'),'source_fit_receipt_sha256':sha(a.fits/'receipt.json'),'structural_sources':sources,'script_sha256':sha(Path(__file__)),'flagged_marker_taxon_links':len(flags),'structural_ready_alignment_exposure_counts':{label:dict(Counter(x['taxon_id'] for x in rows if x['dataset']==label and x['present_in_ready_paired_alignment'])) for label,_,_ in source_specs},'codon_case_exposures':len(codon),'codon_exposure_by_taxon':dict(Counter(x['taxon_id'] for x in codon)),'artifacts':{p.name:sha(p) for p in a.output.iterdir()},'interpretation':'Protein-ID/sequence-hash joins verify flagged marker identities in structural snapshots; exact taxon presence checked in ready structural alignments and completed codon fits. Region overlap is a review signal, not confirmed foreign origin. These inputs need source review and explicit exclusion sensitivities before biological claims; existing results are preserved and no running inputs changed.'}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ['artifacts','structural_sources']},indent=2))


if __name__=='__main__':main()
