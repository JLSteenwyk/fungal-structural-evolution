#!/usr/bin/env python3
"""Run a checksum-bound full overlap comparison of reciprocal ortholog pairs."""
import argparse, hashlib, json, subprocess
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''): h.update(b)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    c=json.loads(a.plan.read_text());out=Path(c['output'])
    if out.exists():raise FileExistsError(out)
    for path,h in c['pins'].items():
        if sha(path)!=h:raise ValueError('Changed pin '+path)
    # Equal files establish equal zero-based protein ordinals and taxon labels.
    for name in ['sequence_ids','species_ids']:
        if sha(c['left'][name])!=sha(c['right'][name]):raise ValueError('Guide ID universes differ')
    receipts=[]
    for label in ['left','right']:
        source=c[label];r=json.loads(Path(source['receipt']).read_text());receipts.append(r)
        audit_plan=json.loads(Path(source['audit_plan']).read_text())
        if sha(source['audit_plan'])!=r['plan_sha256']:raise ValueError('Audit plan mismatch')
        identity_path=audit_plan['identity_plan'];identity=json.loads(Path(identity_path).read_text())
        if sha(identity_path)!=audit_plan['pins'][identity_path]:raise ValueError('Identity plan mismatch')
        if sha(Path(identity['output'])/'receipt.json')!=r['identity_receipt_sha256']:raise ValueError('Identity receipt mismatch')
        snapshot=json.loads(Path(identity['supplement_receipt']).read_text())
        for name in ['sequence_ids','species_ids']:
            if sha(source[name])!=snapshot['input_hashes'][source[name]]:raise ValueError('Source ID mapping mismatch')
        if r['status']!='complete_native_ortholog_pair_multiplicity_audit':raise ValueError('Incomplete audit')
        if any(r[k] for k in ['duplicate_directed_incidences','pairs_missing_reverse','pairs_with_repeated_direction','pairs_with_unequal_multiplicity']):raise ValueError('Invalid pair multiplicities')
        if sha(source['stream'])!=r['sorted_pairs_sha256']:raise ValueError('Changed stream')
    out.mkdir(parents=True)
    binary=out/'compare_pair_streams'
    compile_command=[c['compiler'],'-O3','-std=c++17','-Wall','-Wextra','-pedantic',c['source'],'-o',str(binary)]
    subprocess.run(compile_command,check=True)
    command=[str(binary),c['left']['stream'],c['right']['stream']]
    result=json.loads(subprocess.check_output(command,text=True))
    if result['left_pairs']!=receipts[0]['unique_unordered_pairs'] or result['right_pairs']!=receipts[1]['unique_unordered_pairs']:raise ValueError('Source counts differ')
    assert result['shared_pairs']+result['left_only_pairs']==result['left_pairs']
    assert result['shared_pairs']+result['right_only_pairs']==result['right_pairs']
    for label,r in zip(['left','right'],receipts):
        if sha(c[label]['stream'])!=r['sorted_pairs_sha256']:raise ValueError('Stream changed during comparison')
    for path,h in c['pins'].items():
        if sha(path)!=h:raise ValueError('Pin changed during comparison')
    result.update(status='complete_native_ortholog_pair_guide_overlap',jaccard=result['shared_pairs']/result['union_pairs'] if result['union_pairs'] else None,plan_sha256=sha(a.plan),script_sha256=sha(__file__),compile_command=compile_command,command=command,binary_sha256=sha(binary),source_pins=c['pins'],scope='Exact unordered pair overlap in native tables for identical protein ordinals. Reciprocal multiplicity and sorted uniqueness rechecked during merge. Excludes small-family supplements. Guide-sensitive assignments are not validated errors, duplication events or biological orthology.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
