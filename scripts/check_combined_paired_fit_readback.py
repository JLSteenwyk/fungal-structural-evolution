#!/usr/bin/env python3
"""Known two-source table merge with altered provenance and duplicate-row rejection."""
import csv,json,tempfile
from pathlib import Path
from readback_combined_paired_fit_sources import verify_collection
from run_ortholog_pair_guide_comparison import sha


def write(p,r):p.write_text(json.dumps(r)+'\n')
def table(p,records):
    with p.open('w') as h:
        w=csv.DictWriter(h,fieldnames=list(records[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(records)

with tempfile.TemporaryDirectory() as td:
    root=Path(td);out=root/'out';out.mkdir();sources=[];pins={};combined={n:[] for n in ['fit_summary.tsv','paired_branches.tsv','warnings.tsv']}
    for marker,label in [('M1','unchanged'),('M2','refitted')]:
        fits=root/('paired-fits-'+marker);audit=root/('paired-audit-'+marker);fits.mkdir();audit.mkdir()
        write(fits/'receipt.json',dict(marker=marker));sources.append(dict(marker=marker,source=label,fit_directory=str(fits/marker)))
        data={'fit_summary.tsv':[dict(marker=marker,fit=f,value='1.234') for f in ['aa','3di_af','3di_af_empirical','3di_llm']],
              'paired_branches.tsv':[dict(marker=marker,branch='x',value='0.0001')],
              'warnings.tsv':[dict(marker=marker,warning='fixture')]}
        for name,records in data.items():
            table(audit/name,records)
            combined[name].extend(dict(r,source_run=label,source_fit_directory=str(fits/marker),source_audit_directory=str(audit)) for r in records)
        write(audit/'receipt.json',dict(status='complete_paired_fit_audit',fit_receipt_sha256=sha(fits/'receipt.json'),artifacts={n:sha(audit/n) for n in data}))
        pins[str(audit/'receipt.json')]=sha(audit/'receipt.json');pins[str(fits/'receipt.json')]=sha(fits/'receipt.json')
    inventory=root/'inventory.json';write(inventory,dict(sources=sources));pins[str(inventory)]=sha(inventory)
    for name,data in combined.items():table(out/name,data)
    receipt=dict(status='complete_source_preserving_paired_fit_table_collection',source_sha256=pins,markers=2,paired_branches=2,warning_rows=2,artifacts={n:sha(out/n) for n in combined});write(out/'receipt.json',receipt)
    result=verify_collection(out,inventory);assert result['tables']['fit_summary.tsv']['rows']==8
    for kind in ['provenance','duplicate']:
        records=[dict(r) for r in combined['paired_branches.tsv']]
        if kind=='provenance':records[0]['source_run']='refitted'
        else:records.append(dict(records[0]))
        table(out/'paired_branches.tsv',records);receipt['artifacts']['paired_branches.tsv']=sha(out/'paired_branches.tsv');write(out/'receipt.json',receipt)
        try:verify_collection(out,inventory)
        except ValueError:pass
        else:raise AssertionError('Accepted changed '+kind)
print('Passed two-source all-field readback; updated-hash wrong provenance and duplicate row rejected.')
