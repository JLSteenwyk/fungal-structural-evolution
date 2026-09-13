#!/usr/bin/env python3
"""Read back all NCBI CDS audit tables and reconcile marker classifications."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
from audit_busco_gene_copies import ROOT,sha,read_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable audit summary')
    receipt_path=a.audit/'receipt.json';r=json.loads(receipt_path.read_text());cp=a.audit/'config.json';config=json.loads(cp.read_text())
    if r['status']!='complete_full_ncbi_strict_cds_translation_audit' or sha(cp)!=r['config_sha256']:raise ValueError('Completed full audit required')
    if sha(ROOT/'scripts/audit_ncbi_cds_translation.py')!=config['script_sha256']:raise ValueError('Changed producer source')
    for key,name in [('cds','cds_download_receipts.json'),('gff','annotation_download_receipts.json'),('protein','qc_input_receipts.json')]:
        if sha(ROOT/'metadata'/name)!=config['input_receipts'][key]:raise ValueError('Changed input inventory')
    acquisition=json.loads((ROOT/'metadata/cds_download_receipts.json').read_text());expected={x['taxon_id']:x for x in acquisition['taxa']}
    if len(r['taxa'])!=519 or {x['taxon_id'] for x in r['taxa']}!=set(expected):raise ValueError('Full taxon coverage differs')
    boundary_path=ROOT/'results/cds/ncbi-marker-boundaries-v1/receipt.json';br=json.loads(boundary_path.read_text())
    for name,h in br['artifacts'].items():
        if sha(boundary_path.parent/name)!=h:raise ValueError('Changed marker boundary audit')
    marker_rows=read_table(boundary_path.parent/'marker_annotation_translation.tsv');markers={(x['taxon_id'],x['protein_id']):x for x in marker_rows}
    # One protein can be nominated for more than one marker; assess distinct protein identities.
    summaries=[];totals=Counter();codes=Counter();marker_comparison=[];missing_total=0
    for taxon in r['taxa']:
        t=taxon['taxon_id'];rp=a.audit/(t+'.receipt.json')
        if json.loads(rp.read_text())!=taxon or taxon['config_sha256']!=r['config_sha256']:raise ValueError('Taxon receipt differs from completed batch')
        for name,h in taxon['artifacts'].items():
            if sha(a.audit/name)!=h:raise ValueError('Changed taxon audit artifact')
        counts=Counter();localcodes=Counter();ids=set();protein_ids=set();without=0
        with (a.audit/(t+'.audit.tsv')).open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                if row['taxon_id']!=t:raise ValueError('Mixed taxon audit')
                status=row['status'];pid=row['protein_id']
                if status=='protein_without_cds_record':without+=1
                else:
                    if not row['cds_id'] or row['cds_id'] in ids:raise ValueError('Duplicate/missing CDS record')
                    ids.add(row['cds_id']);counts[status]+=1
                if pid:protein_ids.add(pid)
                if row['code_source']:localcodes[row['code_source']+':'+str(row['translation_table'] or 'None')]+=1
                if (t,pid) in markers:
                    m=markers[t,pid];numeric_agrees=row['translation_table']==m['translation_table']
                    if not numeric_agrees:raise ValueError('Marker translation code disagreement')
                    disposition='same_strict_translation_status' if status==m['translation_status'] else 'full_audit_annotation_exception_gate' if status=='annotation_exception_requires_review' else 'unexplained_disagreement'
                    if disposition=='unexplained_disagreement':raise ValueError('Marker translation status disagreement')
                    marker_comparison.append({'taxon_id':t,'protein_id':pid,'full_audit_status':status,'marker_translation_status':m['translation_status'],'translation_table':row['translation_table'],'comparison':disposition})
        if len(ids)!=taxon['cds_records'] or len(ids)!=expected[t]['cds_records'] or dict(counts)!=taxon['cds_status_counts'] or without!=taxon['proteins_without_cds'] or dict(localcodes)!=taxon['code_source_counts']:raise ValueError('Recounted taxon totals differ')
        totals.update(counts);codes.update(localcodes);missing_total+=without
        summaries.append({'taxon_id':t,'proteins':taxon['proteins'],'cds_records':len(ids),'proteins_without_cds':without,**dict(counts)})
        print(t,len(summaries),flush=True)
    if {(x['taxon_id'],x['protein_id']) for x in marker_comparison}!=set(markers) or len(marker_comparison)!=len(markers):raise ValueError('Marker reconciliation coverage differs')
    a.output.mkdir(parents=True)
    for name,rows in [('taxon_translation_summary.tsv',summaries),('marker_audit_comparison.tsv',marker_comparison)]:
        fields=list(dict.fromkeys(k for row in rows for k in row))
        with (a.output/name).open('w') as f:
            w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result={'status':'complete_full_ncbi_cds_audit_readback','taxa':len(summaries),'cds_records':sum(totals.values()),'protein_records':sum(x['proteins'] for x in summaries),'proteins_without_cds':missing_total,'cds_status_counts':dict(totals),'code_source_counts':dict(codes),'marker_proteins_reconciled':len(marker_comparison),'marker_comparison_counts':dict(Counter(x['comparison'] for x in marker_comparison)),'source_receipt_sha256':sha(receipt_path),'marker_boundary_receipt_sha256':sha(boundary_path),'script_sha256':sha(Path(__file__)),'interpretation':'All source-inventory/config/producer hashes and output table hashes checked; every CDS identity/status and taxon/code total recounted. All indexed marker protein codes and strict statuses reconciled with the independent boundary audit, retaining its different annotation-exception gate. This does not independently retranslate every nonmarker CDS, reconstruct genomes or establish selection eligibility. The original full audit labels omitted codes as assumptions; marker review separately supports NCBI documented defaults and checks region-code fallback.','artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
