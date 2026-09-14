#!/usr/bin/env python3
"""Audit published mean chronograms and exact-name overlap without transferring dates."""
import argparse,csv,hashlib,io,json,math,re,zipfile
from pathlib import Path
from Bio import Phylo


def sha_bytes(b):return hashlib.sha256(b).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    archive=Path('data/raw/timetree_data.zip');blob=archive.read_bytes()
    source=json.loads(Path('metadata/timetree_source.json').read_text())
    files=[x for x in source['files'] if x['name']=='Data.zip']
    if len(files)!=1 or hashlib.md5(blob).hexdigest()!=files[0]['supplied_md5']:raise ValueError('Publisher archive checksum differs')
    z=zipfile.ZipFile(io.BytesIO(blob));members=[]
    for info in z.infolist():
        if not info.is_dir():members.append({'member':info.filename,'bytes':info.file_size,'sha256':sha_bytes(z.read(info.filename))})
    with Path('metadata/analysis_manifest.tsv').open() as f:manifest=list(csv.DictReader(f,delimiter='\t'))
    project={}
    for r in manifest:project.setdefault(r['species_name'].replace(' ','_'),[]).append(r)
    results=[];overlaps=[];member_pins={};checked_edges=0;max_error=0
    for path in sorted(z.namelist()):
        if not path.endswith('.timetree.tree.combined'):continue
        folder=path.rsplit('/',1)[0];scenario=folder.rsplit('/',1)[1]
        name_path=folder+'/short_toFull_speciesNames.tsv'
        names={}
        for line in z.read(name_path).decode().splitlines():
            full,short=line.split('\t')
            if short in names:raise ValueError('Duplicate published tip')
            names[short]=full
        trees=list(Phylo.parse(io.StringIO(z.read(path).decode()),'nexus'))
        if len(trees)!=1:raise ValueError('Expected one annotated summary tree')
        tree=trees[0];tips=[x.name for x in tree.get_terminals()]
        if len(tips)!=len(set(tips)) or not set(tips)<=set(names):raise ValueError('Published name grid differs')
        summary_path=next(n for n in z.namelist() if n.startswith(folder+'/') and n.endswith('.timetree.ages.summary'))
        age_rows=list(csv.DictReader(io.StringIO(z.read(summary_path).decode()),delimiter='\t'));ages={r['Index']:r for r in age_rows}
        if len(ages)!=len(age_rows):raise ValueError('Duplicate summary node index')
        seen=set()
        def index(node):return re.search(r'index=(\d+)',node.comment or '')[1]
        for node in tree.find_clades():
            ident=index(node);seen.add(ident);r=ages[ident]
            vals=[float(r[k]) for k in ['Mean','Variance','Min','Max','95CILower','95CIUpper']]
            if any(not math.isfinite(x) for x in vals) or min(vals)<0:raise ValueError('Invalid age summary')
            if not float(r['Min'])<=float(r['Mean'])<=float(r['Max']):raise ValueError('Invalid age range')
            hpds=re.search(r'age_95%_HPD=\{([^,]+),([^}]+)\}',node.comment)
            if any(abs(float(x)-float(r[k]))>1e-8 for x,k in zip(hpds.groups(),['95CILower','95CIUpper'])):raise ValueError('Tree/summary interval differs')
            for child in node.clades:
                expected=float(r['Mean'])-float(ages[index(child)]['Mean']);err=abs(expected-child.branch_length)
                if expected<0 or err>1e-8:raise ValueError('Mean-age/branch-duration mismatch')
                max_error=max(max_error,err);checked_edges+=1
        if seen!=set(ages):raise ValueError('Node summary grid differs')
        matched=0
        for short,full in sorted(names.items()):
            local=project.get(full,[]);matched+=bool(local) and short in tips
            overlaps.append({'scenario':scenario,'published_tip':short,'published_name':full,'present_in_chronogram':short in tips,'project_taxon_ids':';'.join(r['taxon_id'] for r in local),'project_roles':';'.join(r['study_role'] for r in local),'status':'name_map_entry_not_in_chronogram' if short not in tips else 'exact_name_candidate_requires_identity_review' if len(local)==1 else 'ambiguous_exact_name' if local else 'no_exact_name_match'})
        results.append({'scenario':scenario,'summary_trees':len(trees),'tips':len(tips),'name_map_entries':len(names),'name_map_entries_absent_from_tree':len(set(names)-set(tips)),'nodes':len(ages),'exact_name_matching_published_tips':matched,'posterior_draws_in_this_tree_file':0})
        for n in [path,name_path,summary_path,folder+'/INFO.txt']:member_pins[n]=sha_bytes(z.read(n))
    if len(results)!=4:raise ValueError('Expected four published scenarios')
    a.output.mkdir(parents=True)
    for name,rows in [('archive_members.tsv',members),('scenario_inventory.tsv',results),('taxon_name_overlap.tsv',overlaps)]:
        with (a.output/name).open('w',newline='') as f:
            w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    r={'status':'complete_published_mean_chronogram_inventory','archive_sha256':sha_bytes(blob),'publisher_md5_verified':files[0]['supplied_md5'],'source_doi':source['doi'],'manifest_sha256':sha_bytes(Path('metadata/analysis_manifest.tsv').read_bytes()),'script_sha256':sha_bytes(Path(__file__).read_bytes()),'scenarios':results,'mean_age_edges_checked':checked_edges,'max_mean_age_edge_difference':max_error,'source_members':member_pins,'artifacts':{p.name:sha_bytes(p.read_bytes()) for p in a.output.iterdir()},'interpretation':'Four annotated mean chronograms and marginal age summaries, not joint dated posterior draws. Numeric internal labels are node indices, not branch support. Exact spelling overlap is not accession/strain or clade identity; no dates, priors or durations transferred to project trees. Published pairwise ordering probabilities do not reconstruct full joint posterior branch-duration uncertainty.'}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ['source_members','artifacts']},indent=2))


if __name__=='__main__':main()
