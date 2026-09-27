#!/usr/bin/env python3
"""Bind matched model pairs to prediction provenance and report batch ascertainment."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

p=argparse.ArgumentParser(description=__doc__)
for name in ['context','reference-provenance','local-provenance','output']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();r=json.loads((a.context/'receipt.json').read_text());table=a.context/'source_model_context.tsv'
if r['status']!='complete_paired_source_full_sequence_context' or sha(table)!=r['artifacts'][table.name]:raise ValueError('Changed context')
pins={str(a.context/'receipt.json'):sha(a.context/'receipt.json'),str(table):sha(table)};sources=[]
for path in [a.reference_provenance,a.local_provenance]:
 rp=path.parent/'receipt.json'
 if str(rp) not in r['source_sha256'] or sha(rp)!=r['source_sha256'][str(rp)]:raise ValueError('Wrong source provenance receipt')
 receipt=json.loads(rp.read_text())
 if sha(path)!=receipt['artifacts'][path.name]:raise ValueError('Changed model provenance')
 pins[str(path)]=sha(path);pins[str(rp)]=sha(rp)
 records=json.loads(path.read_text());index={(x['model_id'],str(x['version'])):x for x in records}
 if len(index)!=len(records):raise ValueError('Duplicate provenance identity')
 sources.append(index)
rows=[];seen=set();configs={};counts=Counter();tools=Counter()
for row in csv.DictReader(table.open(),delimiter='\t'):
 key=tuple(row[k] for k in ['reference_model_id','reference_model_version','local_model_id','local_model_version'])
 if key in seen:continue
 seen.add(key);ref,local=sources[0][key[:2]],sources[1][key[2:]]
 for prefix,record in [('reference',ref),('local',local)]:
  if record['sequence_sha256']!=row[prefix+'_sequence_sha256'] or record['sha256']!=row[prefix+'_model_sha256'] or record['path']!=row[prefix+'_model_path']:raise ValueError('Provenance model differs')
 pred=Path(local['prediction_receipt_path']);config=pred.parent/'config.json'
 if sha(pred)!=local['prediction_receipt_sha256'] or sha(config)!=local['prediction_config_sha256']:raise ValueError('Changed prediction/config')
 pins[str(pred)]=sha(pred);pins[str(config)]=sha(config);configdata=json.loads(config.read_text());configs[str(config)]=configdata
 counts[str(config)]+=1;tools[(ref['provider'],ref['tool'],str(ref['version']))]+=1
 rows.append(dict(reference_model_id=key[0],reference_version=key[1],local_model_id=key[2],local_version=key[3],sequence_sha256=ref['sequence_sha256'],reference_provider=ref['provider'],reference_tool=ref['tool'],prediction_receipt=str(pred),prediction_receipt_sha256=sha(pred),local_config=str(config),local_config_sha256=sha(config)))
if len(rows)!=r['distinct_model_pairs']:raise ValueError('Incomplete pair universe')
a.output.mkdir(parents=True,exist_ok=False)
with (a.output/'model_pair_settings.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(a.output/'prediction_configs.json').write_text(json.dumps(configs,indent=2)+'\n')
for path,h in pins.items():
 if sha(path)!=h:raise ValueError('Changed source during inventory')
result=dict(status='complete_overlap_prediction_settings_inventory',model_pairs=len(rows),local_configuration_counts=dict(counts),reference_provider_tool_version_counts=[dict(provider=k[0],tool=k[1],version=k[2],model_pairs=v) for k,v in tools.items()],source_sha256=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in a.output.iterdir()},scope='Model identities and saved local prediction/configuration bytes bound to upstream audited provenance. Batch composition is highly ascertained; no causal settings effect, reproducibility rerun or experimental accuracy claim. AFDB tool/version metadata are not complete native inference configurations.')
(a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['local_configuration_counts'],indent=2));print(result['reference_provider_tool_version_counts'])
