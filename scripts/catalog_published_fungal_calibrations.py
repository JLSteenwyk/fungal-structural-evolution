#!/usr/bin/env python3
"""Catalog published calibration headings as candidates, without approving priors."""
import argparse,csv,hashlib,json,re,subprocess
from pathlib import Path


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    folder=Path('data/literature/fungal-timetree-2025');pdf=folder/'supplementary_sections_1_6.pdf';receipt=json.loads((folder/'supplementary_sections_1_6_receipt.json').read_text())
    if sha(pdf)!=receipt['sha256']:raise ValueError('Source PDF changed')
    process=subprocess.run(['pdftotext','-layout',str(pdf),'-'],capture_output=True,text=True,check=True);text=process.stdout
    with Path('metadata/analysis_manifest.tsv').open() as f:manifest=list(csv.DictReader(f,delimiter='\t'))
    names={}
    for row in manifest:names.setdefault(row['species_name'],[]).append(row['taxon_id'])
    rows=[]
    for m in re.finditer(r'Calibration n[ºo°](\d+[AB]?):',text):
        value=re.match(r'\s*(.*?)\s+(\d+(?:\.\d+)?)(\s*Ma)?[;\s]*(minimum|soft\s+maximum)\s+age',text[m.end():],re.S)
        if value is None or len(value[0])>350:raise ValueError('Unparsed calibration heading '+m[1])
        prefix=' '.join(value[1].split());kind=' '.join(value[4].split())
        stem=re.search(r'(stem|crown)',prefix)
        target=re.split(r'\s*\((?:stem|crown)\)|;\s*(?:stem|crown)',prefix)[0].rstrip(':; ')
        anchors=re.split(r'\s*[–-]\s*',target) if target!='Root' else []
        if target!='Root' and len(anchors)!=2:raise ValueError('Ambiguous published node anchors')
        anchors=[x.replace('_',' ') for x in anchors]
        matches=[names.get(x,[]) for x in anchors]
        rows.append({'calibration_label':m[1],'numeric_label_group':re.match(r'\d+',m[1])[0],
                     'pdf_page_1based':text[:m.start()].count('\f')+1,
                     'anchor_1_published':anchors[0] if anchors else '',
                     'anchor_2_published':anchors[1] if anchors else '',
                     'placement_label':stem[1] if stem else 'root',
                     'bound_type_as_reported':kind.replace(' ','_'),'age_ma':float(value[2]),
                     'ma_explicit_in_heading':bool(value[3]),
                     'anchor_1_exact_project_ids':';'.join(matches[0]) if matches else '',
                     'anchor_2_exact_project_ids':';'.join(matches[1]) if matches else '',
                     'both_anchors_exact_names':len(matches)==2 and all(len(x)==1 for x in matches),
                     'qualification':'candidate_only_specimen_placement_geochronology_prior_and_node_mapping_review_required'})
    if len(rows)!=27 or {r['numeric_label_group'] for r in rows}!={str(i) for i in range(1,25)}:raise ValueError('Expected full 27-heading/24-numbered-group catalogue')
    a.output.mkdir(parents=True)
    path=a.output/'calibration_candidates.tsv'
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    (a.output/'pdf_extraction_stderr.txt').write_text(process.stderr)
    r={'status':'complete_published_calibration_heading_catalogue','source_pdf':receipt,'script_sha256':sha(Path(__file__)),
       'manifest_sha256':sha(Path('metadata/analysis_manifest.tsv')),'calibration_headings':len(rows),'numbered_groups':24,
       'minimum_bounds':sum(r['bound_type_as_reported']=='minimum' for r in rows),'soft_maximum_bounds':sum(r['bound_type_as_reported']=='soft_maximum' for r in rows),
       'both_anchor_names_exactly_matched':sum(r['both_anchors_exact_names'] for r in rows),
       'extraction_warnings_recorded':bool(process.stderr),'artifacts':{p.name:sha(p) for p in a.output.iterdir()},
       'interpretation':'Structured factual transcription of all 27 calibration headings in Supplementary Information 1. A/B suffixes create 27 entries across 24 numbered groups; groups are not assumed to be identical nodes. Ages retain published values, not updated geological bounds. Missing explicit Ma labels use section context and are flagged. No calibration is approved; exact names do not establish accession or stem/crown node equivalence.'}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))


if __name__=='__main__':main()
