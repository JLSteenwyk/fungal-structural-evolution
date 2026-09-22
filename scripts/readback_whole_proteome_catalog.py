#!/usr/bin/env python3
"""Reconstruct all expected proteome/model links from frozen primary inputs."""
import argparse
import csv
import hashlib
import json
import time
from pathlib import Path
import psutil


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
    if plan.get('predecessor'):
        dependency=plan['predecessor']
        while True:
            try:
                proc=psutil.Process(dependency['pid'])
                live=proc.create_time()==dependency['create_time'] and proc.status()!=psutil.STATUS_ZOMBIE
            except psutil.NoSuchProcess:live=False
            if not live:break
            time.sleep(20)
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
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
