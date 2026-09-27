#!/usr/bin/env python3
"""Independently check every one/two-residue background_domain RMSD analytically."""
import csv,json,hashlib,math
from pathlib import Path
import numpy as np
from duplication_alignment_numeric_readback import load_pdb


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    dp=Path('metadata/background_domain_rmsd_diagnostic_plan_20260927.json');plan=json.loads(dp.read_text())
    root=Path(plan['output']);r=json.loads((root/'receipt.json').read_text())
    assert r['status']=='complete_background_domain_rmsd_diagnostic_not_scientific_acceptance'
    assert sha(root/'numeric_readback.tsv')==r['artifacts']['numeric_readback.tsv']
    with (root/'numeric_readback.tsv').open() as f:short=[x for x in csv.DictReader(f,delimiter='\t') if int(x['aligned_length'])<3]
    source=json.loads(Path(plan['source_plan']).read_text());native=Path(source['output']);nr=json.loads((native/'receipt.json').read_text())
    assert sha(native/'receipt.json')==r['producer_receipt_sha256']
    assert sha(native/'checkpoint_manifest.tsv')==nr['artifacts']['checkpoint_manifest.tsv']
    with (native/'checkpoint_manifest.tsv').open() as f:hashes={x['path']:x['sha256'] for x in csv.DictReader(f,delimiter='\t')}
    records=[];needed=set()
    for row in short:
        path=f"pairs/{row['pair_key'][:2]}/{row['pair_key']}-{row['mask']}-{row['order']}.json"
        assert sha(native/path)==hashes[path]
        record=json.loads((native/path).read_text());records.append(record)
        needed.update((x['interval_id'],row['mask']) for x in record['inputs'])
    inputs={}
    for spec in [{'inputs':source['inputs']}]:
        folder=Path(spec['inputs']);ir=json.loads((folder/'receipt.json').read_text());path=folder/'inputs.jsonl'
        assert sha(path)==ir['artifacts']['inputs.jsonl']
        with path.open() as f:
            for line in f:
                row=json.loads(line);key=row['interval_id'],row['mask']
                if key in needed:
                    assert key not in inputs;inputs[key]=row
    assert set(inputs)==needed
    results=[]
    for row,record in zip(short,records):
        coords=[load_pdb(inputs[(x['interval_id'],row['mask'])]) for x in record['inputs']]
        strings=[record['metrics']['alignment_left'],record['metrics']['alignment_right']]
        masks=[np.array(list(x))!='-' for x in strings];paired=masks[0]&masks[1]
        indexes=[(np.cumsum(x)-1)[paired] for x in masks]
        xyz=[c[1][i] for c,i in zip(coords,indexes)];n=len(xyz[0])
        assert n==len(xyz[1])==int(row['aligned_length']) and n in [1,2]
        analytic=0. if n==1 else abs(math.dist(*xyz[0])-math.dist(*xyz[1]))/2
        error=abs(analytic-float(row['rmsd_recomputed']));assert error<1e-10
        native_error=abs(analytic-float(row['rmsd_native']))
        assert (native_error>.00501)==(row['rmsd_status']=='outside_printed_rounding')
        results.append(dict(pair_key=row['pair_key'],mask=row['mask'],order=int(row['order']),aligned_length=n,analytic_minimum_rmsd=analytic,diagnostic_agreement_error=error,rmsd_status=row['rmsd_status'],unique_rotation=False))
    assert len(results)==13
    result=dict(status='passed_all_short_background_domain_alignment_analytic_rmsd_checks',alignments=len(results),one_residue=sum(x['aligned_length']==1 for x in results),two_residue=sum(x['aligned_length']==2 for x in results),outside_printed_rounding=sum(x['rmsd_status']=='outside_printed_rounding' for x in results),maximum_agreement_error=max(x['diagnostic_agreement_error'] for x in results),source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),rows=results,scope='All one/two-pair mappings: direct pair-distance formula verifies minimum RMSD without SVD or quaternion fitting. Rotation is nonunique for every such mapping. This does not accept native discrepant values or qualify the remaining cohort scientifically.')
    out=Path('metadata/background_domain_short_geometry_readback_20260927.json')
    with out.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))


if __name__=='__main__':main()
