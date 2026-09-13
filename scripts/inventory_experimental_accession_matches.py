#!/usr/bin/env python3
"""Inventory experimental PDB entities linked to all frozen model UniProt accessions."""
import argparse,datetime,fcntl,json,time,urllib.request,urllib.parse
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mapping',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    checked_receipt(a.mapping);models=json.loads((a.mapping/'model_provenance.json').read_text())
    accessions=sorted({m['uniprot_accession'] for m in models if m.get('uniprot_accession')})
    if not accessions:raise ValueError('No reference accessions')
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    config={'mapping_receipt_sha256':sha(a.mapping/'receipt.json'),'script_sha256':sha(Path(__file__)),'endpoint':'https://search.rcsb.org/rcsbsearch/v2/query','batch_size':100,'accessions':accessions,'content_type':'experimental','method':'UniProt-accession candidate nomination; exact sequence/coordinates not verified'}
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Changed inventory configuration')
    cp.write_text(json.dumps(config,indent=2)+'\n');ids=set();batches=[]
    for index,start in enumerate(range(0,len(accessions),100)):
        group=accessions[start:start+100];bp=a.output/f'batch-{index:04d}.json';rp=a.output/f'batch-{index:04d}.receipt.json'
        q={'query':{'type':'terminal','service':'text','parameters':{'attribute':'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession','operator':'in','value':group}},'return_type':'polymer_entity','request_options':{'return_all_hits':True,'results_content_type':['experimental']}}
        if rp.exists():
            r=json.loads(rp.read_text())
            if r['config_sha256']!=sha(cp) or r['query']!=q or r['response_sha256']!=sha(bp):raise ValueError('Changed cached response')
        else:
            url=config['endpoint']+'?'+urllib.parse.urlencode({'json':json.dumps(q)})
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(url,timeout=30) as response:status=response.status;raw=response.read()
                    if status not in (200,204):raise ValueError('Unexpected HTTP status')
                    parsed=json.loads(raw) if status==200 else {'total_count':0,'result_set':[]}
                    hits=parsed.get('result_set',[])
                    if parsed.get('total_count')!=len(hits):raise ValueError('Incomplete result pagination')
                    break
                except Exception:
                    if attempt==3:raise
                    time.sleep(2**attempt)
            bp.write_bytes(raw if status==200 else json.dumps(parsed).encode())
            r={'query':q,'config_sha256':sha(cp),'response_sha256':sha(bp),'http_status':status,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'total_count':len(hits)}
            rp.write_text(json.dumps(r,indent=2)+'\n');time.sleep(.25)
        parsed=json.loads(bp.read_text());hits=parsed.get('result_set',[])
        if len(hits)!=parsed['total_count'] or len({h['identifier'] for h in hits})!=len(hits):raise ValueError('Result count/identity mismatch')
        ids.update(h['identifier'] for h in hits);batches.append({'batch':index,'receipt_sha256':sha(rp),'hits':len(hits)});print(index,len(group),len(hits),'entities',flush=True)
    out=a.output/'polymer_entity_ids.txt';out.write_text('\n'.join(sorted(ids))+'\n')
    result={'status':'complete_experimental_accession_candidate_inventory','config_sha256':sha(cp),'input_model_count':len(models),'unique_accessions':len(accessions),'unique_polymer_entities':len(ids),'batches':batches,'interpretation':'Experimental-filtered RCSB polymer entities nominated by accession links only. Exact sequence, engineered constructs, observed residues, method/quality, oligomeric context and prediction-training overlap remain unverified. Missing accession hits do not establish absence of experimental homologs.','artifacts':{out.name:sha(out)}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='batches'},indent=2))


if __name__=='__main__':main()
