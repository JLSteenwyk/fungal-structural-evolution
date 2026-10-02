#!/usr/bin/env python3
"""Check every raw native MSA using separate FASTA parsing and residue counters."""
import argparse
import base64
from collections import Counter
import hashlib
import io
import itertools
import json
from pathlib import Path
import sqlite3
import zlib

from Bio import SeqIO
from full_triad_sequence_sources import load_catalog
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def reconstruct(raw,sequences):
    records=list(SeqIO.parse(io.StringIO(raw.decode('ascii')),'fasta'))
    assert len(records)==3 and {r.id for r in records}=={'m0','m1','m2'}
    strings={r.id:str(r.seq).upper() for r in records}
    assert len({len(s) for s in strings.values()})==1 and len(strings['m0'])>0
    offsets=[]
    for i,original in enumerate(sequences):
        text=strings['m'+str(i)];assert text.replace('-','')==original
        n=0;positions=[]
        for character in text:
            assert character=='-' or character in original
            if character=='-':positions.append(None)
            else:n+=1;positions.append(n)
        assert n==len(original);offsets.append(positions)
    columns=list(zip(*offsets));assert all(any(p is not None for p in t) for t in columns)
    triples=[list(t) for t in columns if all(p is not None for p in t)]
    return dict(alignment_columns=len(columns),common_residues=len(triples),common_full_triples=triples)


def inspect(plan_path,output):
    plan=json.loads(Path(plan_path).read_text());sets,source,bindings=load_catalog(plan,plan_path)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_full_triad_sequence_alignment_dispositions_pending_independent_readback'
    assert receipt['plan_sha256']==sha(plan_path)
    for name,digest in receipt['artifacts'].items():bind(bindings,root/name,digest)
    byid={r['sequence_set_id']:r for r in sets}
    db=sqlite3.connect('file:'+str(root/'native_alignments.sqlite')+'?mode=ro',uri=True)
    assert db.execute('PRAGMA integrity_check').fetchone()==('ok',)
    seen=set();counts=Counter();rows=0;common=0
    for sid,method,permutation,payload,h in db.execute('SELECT * FROM alignments'):
        raw=zlib.decompress(payload);assert hashlib.sha256(raw).hexdigest()==h;r=json.loads(raw)
        key=(sid,method,permutation);assert key not in seen;seen.add(key)
        assert sid in byid and method in ['mafft_auto','famsa_default']
        order=tuple(map(int,permutation));assert order in list(itertools.permutations(range(3)))
        assert tuple(r[k] for k in ['sequence_set_id','method','permutation'])==key
        record=byid[sid];assert r['models']==record['models'] and r['original_lengths']==record['original_lengths']
        fasta=''.join('>m%d\n%s\n'%(i,record['sequences'][i]) for i in order)
        assert r['input_fasta']==fasta and r['input_sha256']==hashlib.sha256(fasta.encode()).hexdigest()
        path=root/'native_work'/sid/(method+'-'+permutation)/'input.faa'
        cmd=([plan['tools']['mafft'],'--amino','--anysymbol','--thread','1','--auto',str(path)] if method=='mafft_auto'
             else [plan['tools']['famsa'],'-t','1',str(path),'STDOUT'])
        assert r['native_command']==cmd
        stdout=base64.b64decode(r['stdout_base64'],validate=True);stderr=base64.b64decode(r['stderr_base64'],validate=True)
        assert hashlib.sha256(stdout).hexdigest()==r['stdout_sha256'] and hashlib.sha256(stderr).hexdigest()==r['stderr_sha256']
        metrics=None
        if r['timed_out']:status='native_timeout'
        elif r['returncode']:status='native_nonzero_exit'
        else:
            try:metrics=reconstruct(stdout,record['sequences']);status='valid_sequence_alignment'
            except (AssertionError,UnicodeError,ValueError,IndexError):status='invalid_native_alignment'
        assert r['status']==status
        for field in ['alignment_columns','common_residues','common_full_triples']:
            assert r[field]==(metrics[field] if metrics is not None else None)
        assert isinstance(r['native_pid'],int) and r['native_pid']>0 and r['native_created']>0 and r['elapsed_seconds']>=0
        if r['timed_out']:assert r['returncode']!=0
        counts[method+':'+status]+=1;common+=r['common_residues'] or 0;rows+=1
    db.close()
    expected={(sid,method,''.join(map(str,order))) for sid in byid for method in ['mafft_auto','famsa_default'] for order in itertools.permutations(range(3))}
    assert seen==expected and rows==source['native_alignment_states']
    summary=dict(ordered_model_triads=source['ordered_model_triads'],source_ready_triads=source['source_ready_triads'],
        unique_sequence_model_sets=len(sets),native_alignment_states=rows,counts=dict(counts),common_residue_occurrences=common)
    assert all(receipt[k]==v for k,v in summary.items())
    bind(bindings,rp);verify(bindings)
    result=dict(status='passed_full_triad_sequence_alignment_raw_msa_readback',**summary,
        plan_sha256=sha(plan_path),producer_receipt_sha256=sha(rp),source_hashes=bindings,
        scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();inspect(a.plan,a.output)
