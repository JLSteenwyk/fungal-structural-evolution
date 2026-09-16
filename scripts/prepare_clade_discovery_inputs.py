#!/usr/bin/env python3
"""Stage complete discovery inputs, sharing only identical guide-defined clades."""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle,'sha256').hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args(); plan=json.loads(args.plan.read_text())
    for name,h in plan['pins'].items():
        if sha(Path(name))!=h: raise ValueError('Changed source '+name)
    output=Path(plan['output'])
    if output.exists(): raise FileExistsError(output)
    if shutil.disk_usage(output.parent).free<plan['resources']['minimum_free_disk_gib']*2**30:
        raise RuntimeError('Insufficient space')
    review=json.loads(Path(plan['review']).read_text())
    with Path(plan['taxon_coverage']).open() as handle:
        rows={int(r['species_id']):r for r in csv.DictReader(handle,delimiter='\t')}
    hashes=json.loads(Path(plan['intermediate_checksums']).read_text())
    sources={}
    for name,h in hashes.items():
        source=Path(name); sid=int(source.stem.split('Species')[1])
        if sha(source)!=h: raise ValueError('Changed intermediate FASTA')
        headers=0; residues=0
        with source.open() as handle:
            for line in handle:
                if line.startswith('>'): headers+=1
                else: residues+=len(line.strip())
        if headers!=int(rows[sid]['intermediate_unassigned_proteins']):
            raise ValueError('Protein count differs')
        sources[sid]=dict(path=str(source),sha256=h,proteins=headers,residues=residues,bytes=source.stat().st_size)
    if set(sources)!=set(rows): raise ValueError('Incomplete source universe')
    clades={}; guide_map={}
    for guide in review['guides']:
        seen=set(); entries=[]
        for chunk in guide['chunks']:
            species=tuple(sorted(chunk['species']))
            if seen.intersection(species): raise ValueError('Overlapping guide clades')
            seen.update(species)
            key='C'+hashlib.sha256(','.join(map(str,species)).encode()).hexdigest()[:16]
            if key in clades and tuple(clades[key]['species'])!=species: raise ValueError('Clade key collision')
            clades[key]=dict(clade_id=key,species=list(species),taxa=len(species),
                             proteins=sum(sources[s]['proteins'] for s in species),
                             residues=sum(sources[s]['residues'] for s in species),
                             ordered_candidate_pairs_including_self=sum(sources[s]['proteins'] for s in species)**2,
                             single_taxon=len(species)==1)
            entries.append(dict(native_clade_index=chunk['clade'],clade_id=key))
        if seen!=set(rows): raise ValueError('Guide coverage incomplete')
        guide_map[guide['guide']]=entries
    output.mkdir(parents=True)
    files=[]
    for key,clade in sorted(clades.items()):
        folder=output/'clades'/key/'proteomes'; folder.mkdir(parents=True)
        for sid in clade['species']:
            source=Path(sources[sid]['path']); dest=folder/(rows[sid]['taxon_id']+'.faa')
            if dest.exists(): raise ValueError('Duplicate taxon filename')
            shutil.copyfile(source,dest)
            if sha(dest)!=sources[sid]['sha256'] or sha(source)!=sources[sid]['sha256']:
                raise ValueError('Copy checksum differs')
            if dest.is_symlink() or (dest.stat().st_dev,dest.stat().st_ino)==(source.stat().st_dev,source.stat().st_ino):
                raise ValueError('Copy aliases source')
            files.append(dict(clade_id=key,species_id=sid,taxon_id=rows[sid]['taxon_id'],relative_path=str(dest.relative_to(output)),**sources[sid]))
    with (output/'files.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(files[0]),delimiter='\t'); writer.writeheader();writer.writerows(files)
    (output/'clades.json').write_text(json.dumps(dict(clades=clades,guide_clades=guide_map),indent=2)+'\n')
    result=dict(status='prepared_complete_unique_clade_discovery_inputs',guides=len(guide_map),
                clades_per_guide={g:len(v) for g,v in guide_map.items()},unique_clades=len(clades),
                shared_clades=len(set(e['clade_id'] for e in guide_map['profile']) & set(e['clade_id'] for e in guide_map['mafft'])),
                unique_single_taxon_clades=sum(c['single_taxon'] for c in clades.values()),
                species=len(sources),distinct_input_proteins=sum(s['proteins'] for s in sources.values()),
                proteins_in_unique_clade_runs=sum(c['proteins'] for c in clades.values()),
                residues_in_unique_clade_runs=sum(c['residues'] for c in clades.values()),
                copied_files=len(files),copied_bytes=sum(f['bytes'] for f in files),
                plan_sha256=sha(args.plan),script_sha256=sha(Path(__file__)),
                artifacts={name:sha(output/name) for name in ['files.tsv','clades.json']},
                interpretation='Full intermediate discovery input retained for both conditional guides; only identical complete taxon sets share a planned run. Native original protein IDs are preserved as FASTA headers. No sequence filtering, singleton removal, clustering or reconciliation performed. Single-taxon clades remain explicit and require an appropriate inference path before launch.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))


if __name__=='__main__': main()
