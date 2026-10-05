#!/usr/bin/env python3
"""Reconstruct all expected proteome/model links from frozen primary inputs."""
import argparse
import csv
import hashlib
import json
import time
import subprocess
from pathlib import Path
import psutil

from run_after_verified_dependencies_v2 import completed


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def fasta(path):
    # Independent parser: do not call the producer's BioPython parser.
    name=None;parts=[]
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:yield name,''.join(parts)
                name=line[1:].split()[0];parts=[]
            else:
                if name is None and line.strip():raise ValueError('Sequence before FASTA header')
                parts.append(''.join(line.split()))
    if name is not None:yield name,''.join(parts)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());started=time.time()
    output=Path(plan['output'])
    if output.exists():raise FileExistsError(output)
    def verify():
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed readback dependency: '+p)
    verify()
    dependency=json.loads(Path(plan['producer_launch']).read_text())
    assert dependency['plan_sha256']==sha(plan['producer_plan'])
    # Outer queued controller has already waited for the original producer.
    # Require its original invocation journal and full native/wrapper terminal
    # payloads; never accept collected systemd success defaults as proof.
    completed(dependency)
    execution_path=Path(plan['producer_execution'])
    adapter_path=Path(plan['producer_adapter_receipt'])
    execution=json.loads(execution_path.read_text());adapter=json.loads(adapter_path.read_text())
    assert execution['status']=='exited_zero_with_receipt' and execution['exit_code']==0
    assert not execution['timed_out'] and execution['receipt_sha256']==sha(adapter_path)
    assert adapter['status']=='completed_full_afdb_catalog_pending_independent_readback'
    assert adapter['scientific_eligibility'] is False
    assert execution['wrapper']==dict(pid=dependency['pid'],created=dependency['created'],cmdline=dependency['cmdline'])
    assert execution['invocation_id']==dependency['invocation_id']
    rows=[json.loads(line) for line in subprocess.check_output(
        ['journalctl','--user','-u',dependency['unit'],'-o','json','--no-pager'],text=True).splitlines()]
    inv=execution['invocation_id']
    original=[r for r in rows if inv in [r.get('_SYSTEMD_INVOCATION_ID'),r.get('USER_INVOCATION_ID')]]
    exact=[r for r in original if r.get('_PID')==str(dependency['pid']) and r.get('_CMDLINE')==' '.join(dependency['cmdline'])]
    assert len(exact)==2
    assert json.loads(exact[0]['MESSAGE'])==dict(original_wrapper=execution['wrapper'],invocation_id=inv)
    assert json.loads(exact[1]['MESSAGE'])=={k:v for k,v in execution.items() if k not in ['source_hashes','artifacts','command','wrapper','child','scope']}
    assert len([r for r in original if r.get('USER_INVOCATION_ID')==inv and 'Started ' in r.get('MESSAGE','')])==1
    assert len([r for r in original if r.get('USER_INVOCATION_ID')==inv and r.get('CPU_USAGE_NSEC')])==1
    pins={}
    for mapping in [plan['pins'],execution['source_hashes'],execution['artifacts'],adapter['source_hashes']]:
        for p,h in mapping.items():
            if p in pins:assert pins[p]==h
            assert sha(p)==h,p
            pins[p]=h
    for p in [args.plan,execution_path,adapter_path,Path(__file__)]:pins[str(p)]=sha(p)
    verify()
    producer=json.loads(Path(plan['producer_plan']).read_text());catalog=Path(producer['output'])
    receipt=json.loads((catalog/'receipt.json').read_text())
    if receipt['status']!='complete_whole_representative_proteome_exact_sequence_catalog' or receipt['plan_sha256']!=sha(plan['producer_plan']):
        raise ValueError('Producer not successfully complete or changed')
    for p,h in producer['pins'].items():
        if sha(p)!=h:raise ValueError('Changed primary input')
    for p,h in receipt['artifacts'].items():
        if sha(catalog/p)!=h:raise ValueError('Changed catalog artifact')
    latest={}
    with open(producer['inventory']) as handle:
        for line in handle:
            record=json.loads(line);latest[record['uniprot_accession']]=record
    eligible=[]
    for accession,record in latest.items():
        if record['status']=='verified':
            eligible.extend(dict(m,uniprot_accession=accession) for m in record['models']
                            if (m.get('provider'),m.get('tool'))==(producer['provider'],producer['tool']))
    del latest
    # Stable descending sort independently implements the producer's ranking.
    eligible.sort(key=lambda m:(m['mean_ca_plddt'],m['version'],m['model_id']),reverse=True)
    best={}
    for model in eligible:best.setdefault(model['sequence_sha256'],model)
    del eligible
    representatives=json.loads(Path(producer['representatives']).read_text())
    with open(producer['manifest']) as handle:
        manifest={r['taxon_id']:r for r in csv.DictReader(handle,delimiter='\t')}
    with (catalog/'taxon_coverage.tsv').open() as handle:
        coverage_rows=list(csv.DictReader(handle,delimiter='\t'))
        coverage={r['taxon_id']:r for r in coverage_rows}
    if len(coverage_rows)!=len(coverage):raise ValueError('Repeated taxon coverage row')
    expected_taxa={r['taxon_id'] for r in representatives['taxa']}
    if set(coverage)!=expected_taxa or set(manifest)!=expected_taxa:raise ValueError('Taxon grid differs')
    screened=linked=0;used=set()
    with (catalog/'protein_model_links.tsv').open() as handle:
        emitted=iter(csv.DictReader(handle,delimiter='\t'))
        for entry in representatives['taxa']:
            if sha(entry['path'])!=entry['sha256']:raise ValueError('Changed protein source')
            count=matched=0;seen=set()
            for protein,sequence in fasta(entry['path']):
                if protein in seen:raise ValueError('Repeated protein')
                seen.add(protein);count+=1
                digest=hashlib.sha256(sequence.encode()).hexdigest();model=best.get(digest)
                if model is None:continue
                if len(sequence)!=model['length']:raise ValueError('Length disagreement')
                actual=next(emitted,None)
                expected=dict(taxon_id=entry['taxon_id'],protein_id=protein,sequence_sha256=digest,
                              model_id=model['model_id'],version=str(model['version']),model_path=model['path'])
                if actual!=expected:raise ValueError('Protein/model link differs from primary inputs')
                matched+=1;used.add(digest)
            c=coverage[entry['taxon_id']];m=manifest[entry['taxon_id']]
            if (count!=entry['selected_proteins'] or int(c['representative_proteins'])!=count
                    or int(c['proteins_with_model'])!=matched or int(c['proteins_without_catalog_model'])!=count-matched
                    or c['species_name']!=m['species_name'] or c['study_role']!=m['study_role']):
                raise ValueError('Per-taxon coverage differs')
            screened+=count;linked+=matched
        if next(emitted,None) is not None:raise ValueError('Unexpected extra protein/model link')
    seen=set()
    with (catalog/'models.jsonl').open() as handle:
        for line in handle:
            model=json.loads(line);digest=model['sequence_sha256']
            if digest not in used or digest in seen or model!=best[digest]:raise ValueError('Model selection/provenance differs')
            seen.add(digest)
    if seen!=used or screened!=receipt['proteins_screened'] or linked!=receipt['proteins_linked'] or len(used)!=receipt['unique_models']:
        raise ValueError('Catalog scope differs')
    if screened!=producer['expected_proteins'] or len(expected_taxa)!=producer['expected_taxa']:raise ValueError('Planned scope differs')
    if receipt['taxa']!=len(coverage) or receipt['taxa_with_models']!=sum(int(c['proteins_with_model'])>0 for c in coverage.values()):
        raise ValueError('Receipt taxon totals differ')
    verify()
    result=dict(status='passed_full_proteome_sequence_and_model_selection_readback',
                taxa=len(expected_taxa),proteins_screened=screened,protein_links=linked,models=len(used),
                producer_receipt_sha256=sha(catalog/'receipt.json'),plan_sha256=sha(args.plan),
                elapsed_seconds=time.time()-started,
                scope='Every representative FASTA and protein/model link reconstructed; model selection independently sorted from frozen latest-status records; all coverage rows checked. Coordinate hashes were checked by producer and are not repeated here. No new coordinate-content or confidence validation.')
    for p,h in pins.items():assert sha(p)==h,p
    result.update(source_hashes=pins,scientific_eligibility=False,
                  original_producer_invocation_verified=True,new_downloads=0,new_predictions=0,gpu=False)
    with output.open('x') as handle:handle.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2),flush=True)


if __name__=='__main__':main()
