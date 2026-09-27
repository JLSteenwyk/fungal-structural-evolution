#!/usr/bin/env python3
"""Extract both Pfam boundary definitions for every candidate focal-domain hit."""
import csv,hashlib,json
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from Bio import SeqIO
from prepare_case_ancestral_neighborhoods import sha,read


def main():
    seqroot=Path('results/ancestral/case-sequences-20260927-v1')
    domains=Path('results/ancestral/case-domain-context-20260927-v1')
    pins={}
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name
        assert sha(p)==r['artifacts'][name];pins[str(rp)]=sha(rp);pins[str(p)]=sha(p);return p
    proof=Path('metadata/ancestral_case_domain_context_completed_20260927.json')
    pr=json.loads(proof.read_text());assert pr['source_receipt_sha256']==sha(domains/'receipt.json')
    assert pr['all_four_policies_agree_for_every_selected_protein'];pins[str(proof)]=sha(proof)
    annotations=json.loads(checked(domains,'candidate_annotations.json').read_text())
    casesroot=Path('results/structural_comparisons/whole-domain-case-dossiers-20260927-v1')
    cases={r['family']:r for r in read(checked(casesroot,'case_dossiers.tsv'))}
    sequences={}
    for family in cases:
        path=checked(seqroot,family+'.faa')
        for record in SeqIO.parse(path,'fasta'):
            assert record.id not in sequences;sequences[record.id]=str(record.seq)
    assert len(sequences)==len(annotations)==1025
    output=Path('results/ancestral/case-domain-sequences-20260927-v1');output.mkdir(exist_ok=False)
    extracted=defaultdict(dict);records=[];coordinate_rows=[]
    for item in annotations:
        family,gene=item['family'],item['gene'];seq=sequences[gene]
        assert item['sequence_id']=='S'+hashlib.sha256(seq.encode()).hexdigest()
        payload=item['candidate_architectures'];target=cases[family]['pfam_accession']
        sets=[]
        for state in payload['policies'].values():
            hits=payload['alternatives'][state['alternative_index']]['annotations']
            sets.append([h for h in hits if h['pfam_accession']==target])
        assert len(sets)==4 and all(h==sets[0] for h in sets)
        hits=sets[0]
        # Preserve a disposition for every protein; never select one repeated
        # domain copy implicitly if future inputs include multiple hits.
        for boundary in ['alignment','envelope']:
            status='no_retained_target_annotation' if not hits else 'multiple_target_hits_unresolved' if len(hits)>1 else 'extracted_single_target_hit'
            row=dict(family=family,gene=gene,boundary=boundary,target_pfam=target,target_hits=len(hits),status=status,source_length=len(seq),source_sha256=hashlib.sha256(seq.encode()).hexdigest(),hit_id='',start='',end='',domain_length='',domain_sha256='',hmm_coverage='',partial_hmm_below_070='')
            if len(hits)==1:
                hit=hits[0];start,end=(int(hit[boundary+'_'+k]) for k in ['start','end'])
                assert 1<=start<=end<=len(seq)
                domain=seq[start-1:end];extracted[family,boundary][gene]=domain
                row.update(hit_id=hit['hit_id'],start=start,end=end,domain_length=len(domain),domain_sha256=hashlib.sha256(domain.encode()).hexdigest(),hmm_coverage=hit['hmm_coverage'],partial_hmm_below_070=int(Decimal(hit['hmm_coverage'])<Decimal('.70')))
                coordinate_rows.extend(dict(family=family,boundary=boundary,gene=gene,domain_position=i+1,protein_position=start+i,residue=aa) for i,aa in enumerate(domain))
            records.append(row)
    assert len(records)==2050
    summaries=[]
    for (family,boundary),data in sorted(extracted.items()):
        assert cases[family]['gene_a'] in data and cases[family]['gene_b'] in data
        name=family+'-'+boundary+'.faa';path=output/name
        path.write_text(''.join('>'+g+'\n'+s+'\n' for g,s in sorted(data.items())))
        assert {r.id:str(r.seq) for r in SeqIO.parse(path,'fasta')}==data
        lengths=list(map(len,data.values()))
        summaries.append(dict(family=family,boundary=boundary,proteins=len(data),taxa=len({g.split('_',1)[0] for g in data}),residues=sum(lengths),minimum_length=min(lengths),maximum_length=max(lengths),fasta=name,sha256=sha(path)))
    assert len(summaries)==26
    for name,data in [('protein_dispositions.tsv',records),('residue_coordinates.tsv',coordinate_rows),('input_summary.tsv',summaries)]:
        with (output/name).open('w') as f:w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
        assert read(output/name)==[{k:str(v) for k,v in row.items()} for row in data]
    # Check all exported residues by direct full-protein indexing, and recover
    # every extracted string from the serialized coordinate table.
    rebuilt=defaultdict(list)
    for row in read(output/'residue_coordinates.tsv'):
        assert row['residue']==sequences[row['gene']][int(row['protein_position'])-1]
        key=(row['family'],row['boundary'],row['gene'])
        assert int(row['domain_position'])==len(rebuilt[key])+1;rebuilt[key].append(row['residue'])
    expected={(f,b,g):s for (f,b),data in extracted.items() for g,s in data.items()}
    assert {k:''.join(v) for k,v in rebuilt.items()}==expected
    for p,h in pins.items():assert sha(p)==h
    result=dict(status='complete_case_domain_sequence_inputs_with_residue_readback',families=13,input_fastas=26,proteins_with_single_target_hit=len({r['gene'] for r in records if r['status']=='extracted_single_target_hit'}),protein_boundary_dispositions=len(records),domain_residue_records=len(coordinate_rows),status_counts={s:sum(r['status']==s for r in records) for s in sorted({r['status'] for r in records})},source_hashes=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in output.iterdir()},scope='All 13 cases and both Pfam alignment/envelope boundaries. All single hits including partial HMM matches retained, no deduplication or residue replacement. Every excluded/no-hit protein explicit. Coordinates checked against original proteins. Inputs for domain sensitivity only, not domain gain/loss, homology certainty or ancestral inference.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']},indent=2))


if __name__=='__main__':main()
