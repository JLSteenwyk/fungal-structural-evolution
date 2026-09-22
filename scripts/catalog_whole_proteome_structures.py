#!/usr/bin/env python3
"""Catalog exact-sequence structure reuse across every representative proteome."""
import argparse
import csv
import hashlib
import json
import time
from collections import Counter
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(8*1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());started=time.time()
    for path,digest in plan['pins'].items():
        if sha(path)!=digest:raise ValueError('Changed input: '+path)
    output=Path(plan['output'])
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    latest={};records=0
    with open(plan['inventory']) as handle:
        for line in handle:
            row=json.loads(line);latest[row['uniprot_accession']]=row;records+=1
    candidates={}
    for accession,row in latest.items():
        if row['status']!='verified':continue
        for model in row['models']:
            if model.get('provider')!=plan['provider'] or model.get('tool')!=plan['tool']:continue
            key=model['sequence_sha256']
            order=(model['mean_ca_plddt'],model['version'],model['model_id'])
            if key not in candidates or order>candidates[key][0]:
                candidates[key]=(order,dict(model,uniprot_accession=accession))
    statuses=dict(Counter(r['status'] for r in latest.values()));del latest
    print('Frozen log records',records,'eligible sequence candidates',len(candidates),flush=True)
    representatives=json.loads(Path(plan['representatives']).read_text())
    if representatives['status']!='complete' or representatives['processed_taxa']!=plan['expected_taxa']:
        raise ValueError('Incomplete representative-proteome inventory')
    with open(plan['manifest']) as handle:
        taxa={r['taxon_id']:r for r in csv.DictReader(handle,delimiter='\t')}
    if set(taxa)!={r['taxon_id'] for r in representatives['taxa']}:raise ValueError('Taxon universe differs')
    selected={};summary=[];screened=linked=0
    with (output/'protein_model_links.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['taxon_id','protein_id','sequence_sha256','model_id','version','model_path'])
        for entry in representatives['taxa']:
            path=Path(entry['path'])
            if sha(path)!=entry['sha256']:raise ValueError('Changed representative FASTA')
            seen=set();count=matches=0
            with path.open() as fasta:
                for title,sequence in SimpleFastaParser(fasta):
                    protein=title.split()[0]
                    if protein in seen:raise ValueError('Repeated protein identity')
                    seen.add(protein);count+=1
                    digest=hashlib.sha256(sequence.encode()).hexdigest()
                    if digest not in candidates:continue
                    model=candidates[digest][1]
                    if model['length']!=len(sequence):raise ValueError('Sequence length differs')
                    key=model['path']
                    if key in selected and selected[key]!=model:raise ValueError('Conflicting model provenance')
                    selected[key]=model;matches+=1
                    writer.writerow([entry['taxon_id'],protein,digest,model['model_id'],model['version'],key])
            if count!=entry['selected_proteins']:raise ValueError('Representative count differs')
            screened+=count;linked+=matches
            taxon=taxa[entry['taxon_id']]
            summary.append(dict(taxon_id=entry['taxon_id'],species_name=taxon['species_name'],
                                study_role=taxon['study_role'],representative_proteins=count,
                                proteins_with_model=matches,proteins_without_catalog_model=count-matches))
    if screened!=plan['expected_proteins']:raise ValueError('Protein universe differs')
    if len({(m['model_id'],m['version']) for m in selected.values()})!=len(selected):
        raise ValueError('Model identity occurs under multiple paths')
    model_bytes=0
    with (output/'models.jsonl').open('w') as handle:
        for number,(path,model) in enumerate(sorted(selected.items()),1):
            if sha(path)!=model['sha256']:raise ValueError('Changed coordinate file: '+path)
            model_bytes+=Path(path).stat().st_size
            handle.write(json.dumps(model,separators=(',',':'))+'\n')
            if number%10000==0:print('Coordinate hashes verified',number,'of',len(selected),flush=True)
    with (output/'taxon_coverage.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(summary[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(summary)
    for path,digest in plan['pins'].items():
        if sha(path)!=digest:raise ValueError('Input changed during cataloging')
    result=dict(status='complete_whole_representative_proteome_exact_sequence_catalog',
                taxa=len(summary),proteins_screened=screened,proteins_linked=linked,
                unique_models=len(selected),taxa_with_models=sum(r['proteins_with_model']>0 for r in summary),
                coordinate_bytes_verified=model_bytes,inventory_records=records,
                latest_accession_status_counts=statuses,plan_sha256=sha(args.plan),
                elapsed_seconds=time.time()-started,artifacts={p.name:sha(p) for p in output.iterdir()},
                scope='All representative proteins screened by exact sequence hash and length; selected coordinate bytes checked against existing verification records. One source-specific model per sequence. No new CIF-content audit, confidence qualification, PAE assessment, structure clustering, orthology assignment or experimental validation. Missing means absent from this frozen source catalog, not public-database absence.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
