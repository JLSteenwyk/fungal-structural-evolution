#!/usr/bin/env python3
"""Inventory every existing BUSCO/HMMER alternative for FCS-overlapping markers."""
import argparse,csv,json,hashlib,math
from collections import defaultdict,Counter
from pathlib import Path
from Bio import SeqIO


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def merged_length(intervals):
    total=0;right=0
    for lo,hi in sorted(intervals):
        total+=max(0,hi-max(right,lo-1));right=max(right,hi)
    return total


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable output')
    paths={'mapping':Path('results/qc/fcs-cds-overlap-v1/receipt.json'),'audit':Path('results/qc/fcs-cds-overlap-audit-v1/receipt.json'),'markers':Path('results/phylogeny/markers-full-v1/receipt.json'),'representatives':Path('metadata/gene_representatives_receipt.json'),'proteins':Path('metadata/qc_input_receipts.json'),'cutoffs':Path('data/busco_downloads/lineages/eukaryota_odb12.2/scores_cutoff')}
    mapping=json.loads(paths['mapping'].read_text());audit=json.loads(paths['audit'].read_text())
    if audit['mapping_receipt_sha256']!=sha(paths['mapping']):raise ValueError('Mismatched source audit')
    root=paths['mapping'].parent
    for name,digest in mapping['artifacts'].items():
        if sha(root/name)!=digest:raise ValueError('Changed overlap artifact')
    flags=[r for r in table(root/'marker_overlap_review.tsv') if r['fcs_overlap_status']=='cds_overlaps_fcs_region'];overlaps={(r['taxon_id'],r['protein_id']):r for r in table(root/'protein_overlap_review.tsv')}
    marker_sources={r['taxon_id']:r for r in json.loads(paths['markers'].read_text())['source_tables']};reps={r['taxon_id']:r for r in json.loads(paths['representatives'].read_text())['taxa']};proteins={r['taxon_id']:r for r in json.loads(paths['proteins'].read_text())};cutoffs={s.split()[0]:float(s.split()[1]) for s in paths['cutoffs'].read_text().splitlines() if s.strip()}
    cache={};rows=[];summaries=[];proofs=[]
    for flag in flags:
        taxon,marker=flag['taxon_id'],flag['marker'];source=marker_sources[taxon]
        if taxon not in cache:
            bp=Path(source['path']);dp=Path(reps[taxon]['path']).with_suffix('.decisions.tsv');pp=Path(proteins[taxon]['input_path'])
            if sha(bp)!=source['sha256'] or sha(dp)!=reps[taxon]['decisions_sha256'] or sha(pp)!=proteins[taxon]['sha256']:raise ValueError('Changed BUSCO/protein/decision source')
            final=defaultdict(dict)
            for line in bp.read_text().splitlines():
                if line and not line.startswith('#'):
                    x=line.split('\t')
                    if x[1]!='Missing':final[x[0]][x[2]]=x[1]
            cache[taxon]=(final,{x['protein_id']:x for x in table(dp)},{x.id:str(x.seq) for x in SeqIO.parse(pp,'fasta')})
        final,decisions,seqs=cache[taxon]
        if flag['protein_id'] not in final[marker]:raise ValueError('Selected protein absent from original BUSCO call')
        raw=Path(source['path']).parent/'hmmer_output/initial_run_results'/(marker+'.out');text=raw.read_text()
        if '# Program:         hmmsearch' not in text or not text.rstrip().endswith('# [ok]'):raise ValueError('Incomplete HMMER output')
        model=Path('data/busco_downloads/lineages/eukaryota_odb12.2/hmms')/(marker+'.hmm')
        if str(model.resolve()) not in text or str(Path(proteins[taxon]['input_path']).resolve()) not in text:raise ValueError('HMMER header source differs')
        hits=defaultdict(list)
        for line in text.splitlines():
            if not line.strip() or line.startswith('#'):continue
            x=line.split(maxsplit=22)
            if len(x)!=23 or x[3]!=marker.split('at')[0]:raise ValueError('Wrong HMMER query or column count')
            hits[x[0]].append(x)
        if not set(final[marker])<=set(hits):raise ValueError('Final BUSCO hit absent from raw HMMER report')
        local=[]
        for pid,domains in sorted(hits.items()):
            if pid not in seqs or pid not in decisions:raise ValueError('Unknown raw hit protein')
            fixed={(x[2],x[5],x[6],x[7]) for x in domains}
            if len(fixed)!=1:raise ValueError('Conflicting sequence-level HMMER fields')
            tlen,qlen,evalue,score=next(iter(fixed));tlen,qlen=int(tlen),int(qlen);score=float(score);evalue=float(evalue)
            if tlen!=len(seqs[pid]) or not math.isfinite(score) or not math.isfinite(evalue):raise ValueError('Invalid HMMER length/score')
            hmm=[];aligned=[]
            for x in domains:
                hs,he,ts,te=map(int,x[15:19])
                if not 1<=hs<=he<=qlen or not 1<=ts<=te<=tlen:raise ValueError('Invalid domain coordinates')
                hmm.append((hs,he));aligned.append((ts,te))
            d=decisions[pid];ov=overlaps.get((taxon,pid));isalt=pid!=flag['protein_id'];unflagged=ov is None
            row={'taxon_id':taxon,'marker':marker,'protein_id':pid,'selected_original_marker':not isalt,'in_final_busco_table':pid in final[marker],'final_busco_status':final[marker].get(pid,''),'recorded_fcs_cds_overlap':not unflagged,'fcs_actions':ov['fcs_actions'] if ov else '', 'full_sequence_bit_score':score,'dataset_score_cutoff':cutoffs[marker],'score_minus_dataset_cutoff':score-cutoffs[marker],'full_sequence_evalue':evalue,'profile_length':qlen,'protein_length':tlen,'reported_domain_rows':len(domains),'profile_union_coverage_fraction':merged_length(hmm)/qlen,'protein_union_coverage_fraction':merged_length(aligned)/tlen,'gene_ids_json':d['gene_ids_json'],'gene_mapping_status':d['status'],'representative_decision':d['decision'],'sequence_sha256':hashlib.sha256(seqs[pid].encode()).hexdigest(),'alternative_no_recorded_overlap_at_or_above_score_cutoff':isalt and unflagged and score>=cutoffs[marker]}
            rows.append(row);local.append(row)
        summaries.append({'taxon_id':taxon,'marker':marker,'selected_protein_id':flag['protein_id'],'raw_hit_proteins':len(local),'alternative_raw_hits':sum(not r['selected_original_marker'] for r in local),'alternative_final_calls_without_overlap':sum(not r['selected_original_marker'] and r['in_final_busco_table'] and not r['recorded_fcs_cds_overlap'] for r in local),'alternative_raw_hits_without_overlap':sum(not r['selected_original_marker'] and not r['recorded_fcs_cds_overlap'] for r in local),'alternative_raw_hits_without_overlap_at_or_above_cutoff':sum(r['alternative_no_recorded_overlap_at_or_above_score_cutoff'] for r in local)})
        proofs.append({'taxon_id':taxon,'marker':marker,'raw_hmmer_path':str(raw),'raw_hmmer_sha256':sha(raw),'model_sha256':sha(model),'busco_table_sha256':source['sha256']})
    a.output.mkdir(parents=True)
    for name,data in [('all_raw_hit_review.tsv',rows),('marker_alternative_summary.tsv',summaries),('score_passing_unflagged_alternatives.tsv',[r for r in rows if r['alternative_no_recorded_overlap_at_or_above_score_cutoff']])]:
        with (a.output/name).open('w',newline='') as f:
            w=csv.DictWriter(f,list(rows[0]) if 'summary' not in name else list(summaries[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    selected=[r for r in rows if r['alternative_no_recorded_overlap_at_or_above_score_cutoff']]
    result={'status':'complete_full_flagged_marker_existing_hit_review','affected_marker_taxon_observations':len(flags),'raw_marker_protein_hits':len(rows),'alternative_final_calls_without_overlap':sum(r['alternative_final_calls_without_overlap'] for r in summaries),'alternative_raw_hits_without_overlap':sum(r['alternative_raw_hits_without_overlap'] for r in summaries),'alternative_raw_hits_without_overlap_at_or_above_dataset_score_cutoff':len(selected),'affected_observations_with_score_passing_unflagged_alternative':len({(r['taxon_id'],r['marker']) for r in selected}),'candidate_counts_by_taxon':dict(Counter(r['taxon_id'] for r in selected)),'source_hashes':{k:sha(v) for k,v in paths.items()},'script_sha256':sha(Path(__file__)),'raw_source_checks':proofs,'artifacts':{p.name:sha(p) for p in a.output.iterdir()},'interpretation':'All existing raw HMMER hits and final BUSCO calls reviewed for every FCS-overlapping selected marker. Score cutoff comparison alone does not reproduce BUSCO acceptance, establish coverage/orthology, prove a hit free of contamination, or authorize replacement. Candidates require full gene-copy/phylogenetic, CDS and alignment review; no marker sequence was changed or inferred.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['raw_source_checks','artifacts','source_hashes']},indent=2))


if __name__=='__main__':main()
