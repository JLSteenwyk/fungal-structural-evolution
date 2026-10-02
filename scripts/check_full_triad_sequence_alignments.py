#!/usr/bin/env python3
"""Run actual pinned aligners over indel/repeat/identical/unknown-letter fixtures."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import sqlite3
import zlib

import run_full_triad_sequence_alignments as producer
import readback_full_triad_sequence_alignments as reader
from triad_sequence_native import execute,parse_alignment
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--tools',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=False)
    sequences=[['ACDEFGHIK','ACEFGHIK','AACDEFGHIK'],['ACACACDE','ACACDE','ACACACACDE'],
               ['ACDEFG','ACDEFG','ACDEFG'],['ACDUOX','ACDOX','AACDUOX']]
    sets=[dict(sequence_set_id='fixture'+str(i),models=[['M'+str(i)+str(j),1] for j in range(3)],
        sequences=seq,original_lengths=list(map(len,seq)),sequence_sha256=[hashlib.sha256(s.encode()).hexdigest() for s in seq]) for i,seq in enumerate(sequences)]
    source=dict(ordered_model_triads=5,source_ready_triads=4,native_alignment_states=48)
    def load(*args):return sets,source,{}
    producer.load_catalog=load;reader.load_catalog=load
    plan=dict(tools=json.loads(a.tools.read_text()),output=str(out/'native'),workers=2,native_timeout_seconds=30,
        resources=dict(minimum_free_disk_gib=0,output_allowance_gib=2),scope='Actual native software fixtures only')
    pp=out/'plan.json';pp.write_text(json.dumps(plan,indent=2)+'\n');producer.run(pp);reader.inspect(pp,out/'readback.json')
    root=Path(plan['output']);dbpath=root/'native_alignments.sqlite';rp=root/'receipt.json'
    original_db=dbpath.read_bytes();original_receipt=rp.read_bytes();rejected=0
    fields=['input_fasta','input_sha256','models','original_lengths','native_command','method','permutation',
            'returncode','timed_out','status','alignment_columns','common_residues','common_full_triples','stdout_sha256']
    for field in fields:
        db=sqlite3.connect(dbpath);key,payload=db.execute('SELECT rowid,payload FROM alignments WHERE sequence_set_id="fixture0" LIMIT 1').fetchone()
        r=json.loads(zlib.decompress(payload));value=r[field]
        r[field]=['changed'] if isinstance(value,list) else (not value if isinstance(value,bool) else value+1 if isinstance(value,(int,float)) else 'changed')
        raw=json.dumps(r,separators=(',',':')).encode();db.execute('UPDATE alignments SET payload=?,payload_sha256=? WHERE rowid=?',(zlib.compress(raw),hashlib.sha256(raw).hexdigest(),key));db.commit();db.close()
        receipt=json.loads(original_receipt);receipt['artifacts'][dbpath.name]=sha(dbpath);rp.write_text(json.dumps(receipt))
        try:reader.inspect(pp,out/('changed-'+str(rejected)+'.json'))
        except (AssertionError,KeyError,TypeError,ValueError):rejected+=1
        else:raise AssertionError('Altered native checkpoint accepted '+field)
        finally:dbpath.write_bytes(original_db);rp.write_bytes(original_receipt)
    db=sqlite3.connect(dbpath);db.execute('DELETE FROM alignments WHERE rowid=(SELECT rowid FROM alignments LIMIT 1)');db.commit();db.close()
    receipt=json.loads(original_receipt);receipt['artifacts'][dbpath.name]=sha(dbpath);rp.write_text(json.dumps(receipt))
    try:reader.inspect(pp,out/'missing-row.json')
    except AssertionError:rejected+=1
    else:raise AssertionError('Missing native order accepted')
    finally:dbpath.write_bytes(original_db);rp.write_bytes(original_receipt)
    assert rejected==15
    # An interrupted uncommitted task may leave its exact input in the owned
    # working directory; replay that task without overwriting foreign content.
    record=sets[0];restart=out/'interrupted-task';folder=restart/'native_work'/record['sequence_set_id']/'famsa_default-012'
    folder.mkdir(parents=True);folder.joinpath('input.faa').write_text(''.join('>m%d\n%s\n'%(i,s) for i,s in enumerate(record['sequences'])))
    replay=execute(plan,record,'famsa_default',(0,1,2),restart);assert replay['status']=='valid_sequence_alignment'
    # Altered or all-gap raw alignments cannot become coordinate correspondences.
    for raw in [b'>m0\nACDEFGHIK\n>m1\nACDEFGHIK\n>m2\nAACDEFGHIK\n',
                b'>m0\n-ACDEFGHIK\n>m1\n-ACDEFGHIK\n>m2\n-AACDEFGHIK\n']:
        for parser in [parse_alignment,reader.reconstruct]:
            try:parser(raw,sets[0]['sequences'])
            except AssertionError:pass
            else:raise AssertionError('Changed raw sequence accepted')
    counts=json.loads((root/'receipt.json').read_text())['counts']
    assert sum(v for k,v in counts.items() if k.endswith(':valid_sequence_alignment'))==48
    files=['scripts/full_triad_sequence_sources.py','scripts/triad_sequence_native.py',
        'scripts/run_full_triad_sequence_alignments.py','scripts/readback_full_triad_sequence_alignments.py',str(Path(__file__)),str(a.tools)]
    result=dict(status='passed_full_triad_sequence_alignment_native_software_contracts',actual_native_states=48,
        methods=2,input_orders=6,fixtures=4,altered_exports_rejected=rejected,
        exact_interrupted_task_replay_passed=True,counts=counts,scientific_eligibility=False,
        source_hashes={s:sha(s) for s in files},scope='Actual installed MAFFT/FAMSA software fixtures, all6orders, indel/repeat/identical/nonstandard-letter sequences. Full raw independent SeqIO/counter readback, changed checkpoints/scope reject, interrupted exact input replay. No biological pilot or production alignment completion.')
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
