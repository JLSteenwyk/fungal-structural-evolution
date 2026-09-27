#!/usr/bin/env python3
"""Verify sequence-context distinctions against all selected full-mask PDBs."""
import argparse,csv,json,hashlib
from collections import Counter
from pathlib import Path
from Bio.Data.IUPACData import protein_letters_3to1


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not a.output.exists();p=json.loads(a.plan.read_text());root=Path(p['output']);r=json.loads((root/'receipt.json').read_text())
    assert r['plan_sha256']==sha(a.plan)
    for f,h in p['pins'].items():assert sha(f)==h
    for f,h in r['artifacts'].items():assert sha(root/f)==h
    source=read(p['candidates']);output=read(root/'candidate_sequence_context.tsv');keys=['family','gene_a','gene_b','pfam_accession']
    source={tuple(x[k] for k in keys):x for x in source};indexed={tuple(x[k] for k in keys):x for x in output};assert len(indexed)==len(output)==len(source)==r['candidates'] and set(indexed)==set(source)
    details=read(root/'interval_pair_context.tsv');d={x['triad_key']:x for x in details};assert len(d)==len(details)==r['triads']
    expected={t for x in source.values() for t in json.loads(x['triad_keys_json'])};assert set(d)==expected
    triads={x['triad_key']:x for x in map(json.loads,Path(p['triads']).open()) if x['triad_key'] in expected}
    ids={x[f'interval_{role}'] for x in d.values() for role in ['a','b']};assert len(ids)==r['intervals']
    intervals={x['interval_id']:x for x in map(json.loads,Path(p['intervals']).open()) if x['interval_id'] in ids}
    inp={x['interval_id']:x for x in map(json.loads,Path(p['inputs']).open()) if x['interval_id'] in ids and x['mask']=='full'}
    sequences={}
    for iid,x in inp.items():
        assert sha(x['path'])==x['sha256'];letters=[];positions=[]
        for line in Path(x['path']).read_text().splitlines():
            if line.startswith('ATOM  '):
                assert line[12:16].strip()=='CA';letters.append(protein_letters_3to1[line[17:20].strip().title()]);positions.append(int(line[22:26]))
        sequence=''.join(letters);assert sequence==x['sequence'] and positions==x['original_positions'];assert len(sequence)==intervals[iid]['length'];sequences[iid]=sequence
    for tk,x in d.items():
        for role in ['a','b']:
            iid=triads[tk]['interval_'+role];assert x['interval_'+role]==iid;i=intervals[iid]
            for k in ['model_id','version','start','end','length','original_length','sequence_sha256','source_sha256']:assert x[role+'_'+k]==str(i[k])
            assert x[role+'_interval_sequence_sha256']==hashlib.sha256(sequences[iid].encode()).hexdigest()
        assert int(x['complete_intervals_identical'])==int(sequences[x['interval_a']]==sequences[x['interval_b']])
        assert int(x['full_chain_sequence_hashes_identical'])==int(intervals[x['interval_a']]['sequence_sha256']==intervals[x['interval_b']]['sequence_sha256'])
        assert int(x['same_model_version'])==int((x['a_model_id'],x['a_version'])==(x['b_model_id'],x['b_version']))
    counts=Counter()
    for key,x in indexed.items():
        original=source[key]
        for f,v in original.items():assert x[f]==v
        pool=[d[t] for t in json.loads(original['triad_keys_json'])]
        core=int(float(original['all_alternatives_min_sequence_identity_ab'])==1);n=len(pool);equal=sum(int(t['complete_intervals_identical']) for t in pool)
        full={int(t['full_chain_sequence_hashes_identical']) for t in pool};same={int(t['same_model_version']) for t in pool};assert len(full)==len(same)==1
        expected_values=dict(all_aligned_cores_identical=core,complete_interval_pairs=n,identical_complete_interval_pairs=equal,all_complete_intervals_identical=int(equal==n),full_chain_sequence_hashes_identical=next(iter(full)),same_model_version=next(iter(same)))
        for f,v in expected_values.items():assert int(x[f])==v
        for role in ['a','b']:
            for f in ['model_id','version','original_length','sequence_sha256']:assert {t[role+'_'+f] for t in pool}=={x[role+'_'+f]}
        counts[(x['candidate_class'],x['study_role'],core,int(equal==n),int(x['full_chain_sequence_hashes_identical']))]+=1
    expected_counts={(x['candidate_class'],x['study_role'],x['all_aligned_cores_identical'],x['all_complete_intervals_identical'],x['full_chain_sequence_hashes_identical']):x['count'] for x in r['strata']};assert dict(counts)==expected_counts
    result=dict(status='passed_full_candidate_sequence_context_readback',candidates=len(output),triads=len(d),full_mask_pdbs_checked=len(sequences),all_aligned_cores_identical=sum(int(x['all_aligned_cores_identical']) for x in output),all_complete_intervals_identical=sum(int(x['all_complete_intervals_identical']) for x in output),full_chain_hashes_identical=sum(int(x['full_chain_sequence_hashes_identical']) for x in output),identical_intervals_but_not_all_aligned_cores_identical=sum(x['all_complete_intervals_identical']=='1' and x['all_aligned_cores_identical']=='0' for x in output),plan_sha256=sha(a.plan),producer_receipt_sha256=sha(root/'receipt.json'),checker_sha256=sha(__file__),scope='Every candidate field and triad membership checked. All selected full-mask PDB sequence letters, residue positions and hashes read independently; interval equality and classification strata reconstructed. Full-chain equality uses previously validated source sequence hashes, not new full-chain coordinate parsing.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
