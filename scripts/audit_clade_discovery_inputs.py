#!/usr/bin/env python3
"""Read back every staged clade FASTA and both complete guide partitions."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha(p):
    with p.open('rb') as h:
        return hashlib.file_digest(h,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    plan=json.loads(args.plan.read_text()); root=Path(plan['output'])
    receipt=json.loads((root/'receipt.json').read_text())
    if sha(args.plan)!=receipt['plan_sha256']: raise ValueError('Plan differs')
    for name,h in receipt['artifacts'].items():
        if sha(root/name)!=h: raise ValueError('Artifact differs')
    data=json.loads((root/'clades.json').read_text())
    review=json.loads(Path(plan['review']).read_text())
    if sha(Path(plan['review']))!=plan['pins'][plan['review']]: raise ValueError('Review changed')
    with (root/'files.tsv').open() as h: files=list(csv.DictReader(h,delimiter='\t'))
    actual={str(p.relative_to(root)) for p in (root/'clades').rglob('*') if p.is_file()}
    expected={r['relative_path'] for r in files}
    if actual!=expected or len(expected)!=len(files): raise ValueError('Copied file universe differs')
    ids_by_clade={c:set() for c in data['clades']}; seen_files=set(); sequence_total=0
    for row in files:
        dest=root/row['relative_path']; source=Path(row['path']); sid=int(row['species_id']); clade=row['clade_id']
        if (clade,sid) in seen_files: raise ValueError('Duplicate species input')
        seen_files.add((clade,sid))
        if sha(dest)!=row['sha256'] or sha(source)!=row['sha256']: raise ValueError('Source/copy mismatch')
        if dest.is_symlink() or (dest.stat().st_dev,dest.stat().st_ino)==(source.stat().st_dev,source.stat().st_ino):
            raise ValueError('Source alias')
        ids=set(); residues=0
        with dest.open() as h:
            for line in h:
                if line.startswith('>'):
                    gene=line[1:].strip()
                    if gene in ids or int(gene.split('_')[0])!=sid: raise ValueError('Duplicate/wrong species header')
                    ids.add(gene)
                else: residues+=len(line.strip())
        if len(ids)!=int(row['proteins']) or residues!=int(row['residues']): raise ValueError('Sequence count differs')
        if ids_by_clade[clade].intersection(ids): raise ValueError('Clade contains duplicate proteins')
        ids_by_clade[clade].update(ids); sequence_total+=len(ids)
    for key,c in data['clades'].items():
        if {s for k,s in seen_files if k==key}!=set(c['species']) or len(ids_by_clade[key])!=c['proteins']:
            raise ValueError('Clade membership differs')
    per_guide={}; unions=[]
    for guide in review['guides']:
        name=guide['guide']; mapping=data['guide_clades'][name]; union=set()
        if len(mapping)!=len(guide['chunks']): raise ValueError('Clade count differs')
        for entry,original in zip(mapping,guide['chunks']):
            key=entry['clade_id']
            if entry['native_clade_index']!=original['clade'] or data['clades'][key]['species']!=original['species']:
                raise ValueError('Guide clade membership differs')
            if union.intersection(ids_by_clade[key]): raise ValueError('Within-guide protein duplication')
            union.update(ids_by_clade[key])
        if len(union)!=receipt['distinct_input_proteins']: raise ValueError('Guide misses discovery proteins')
        per_guide[name]=len(union); unions.append(union)
    if any(u!=unions[0] for u in unions[1:]): raise ValueError('Guides have different input protein universes')
    result=dict(status='passed_complete_discovery_input_readback',files=len(files),unique_clades=len(data['clades']),
                protein_records=sequence_total,unique_proteins_per_guide=per_guide,
                receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),
                scope='All file digests, source isolation, sequence identifiers, residues and counts checked; exact guide clade membership and complete identical protein universes verified. No family inference or biological homology validation.')
    args.output.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))


if __name__=='__main__':main()
