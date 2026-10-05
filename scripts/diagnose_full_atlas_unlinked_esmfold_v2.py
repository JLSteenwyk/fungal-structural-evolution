#!/usr/bin/env python3
"""Trace every source-unlinked ESMFold model to original gene-product policy."""
import argparse
import csv
from datetime import datetime,timezone
import json
import hashlib
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pathlib import Path
import sqlite3

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists()
    plan_path=Path('metadata/full_prediction_atlas_union_plan_20261005_v1.json')
    plan=json.loads(plan_path.read_text());root=Path(plan['output'])
    producer_path=Path('metadata/full_prediction_atlas_union_20261005_v1.json')
    transport_path=Path('metadata/full_prediction_atlas_union_transport_20261005_v1.json')
    producer,transport=[json.loads(q.read_text()) for q in [producer_path,transport_path]]
    assert producer['status']=='completed_full_prediction_atlas_union_pending_independent_readback'
    assert transport['validation_sha256']==sha(producer_path)
    assert transport['original_tool_terminal_exit_code']==0
    pins={};bind(pins,root/'prediction_atlas.sqlite',producer['source_hashes'][str(root/'prediction_atlas.sqlite')])
    db=sqlite3.connect('file:'+str((root/'prediction_atlas.sqlite').resolve())+'?mode=ro',uri=True)
    unused={row[0] for row in db.execute("SELECT sequence_sha256 FROM models WHERE source='ESMFold' EXCEPT SELECT sequence_sha256 FROM proteins WHERE has_esmfold=1")}
    assert len(unused)==producer['esmfold_models_without_representative_links']
    marker_path=Path('results/phylogeny/markers-full-v1/protein_mapping.tsv')
    representatives=json.loads(Path(plan['representatives']).read_text())
    entries={r['taxon_id']:r for r in representatives['taxa']}
    with marker_path.open() as h: origins=[r for r in csv.DictReader(h,delimiter='\t') if r['sequence_sha256'] in unused]
    assert {r['sequence_sha256'] for r in origins}==unused
    diagnostics=[]
    for origin in origins:
        taxon=origin['taxon_id'];entry=entries[taxon]
        qc_path=Path('metadata/qc_input_receipts.json');bind(pins,qc_path)
        qc=next(r for r in json.loads(qc_path.read_text()) if r['taxon_id']==taxon)
        assert qc['sha256']==entry['source_sha256']
        source_path=Path(qc['input_path']);bind(pins,source_path,entry['source_sha256'])
        marker_fasta=Path(f"results/busco/{taxon}/run_eukaryota_odb12.2/busco_sequences/single_copy_busco_sequences/{origin['marker']}.faa")
        bind(pins,marker_fasta,origin['source_fasta_sha256'])
        with marker_fasta.open() as h:marker_records=list(SimpleFastaParser(h))
        assert len(marker_records)==1 and marker_records[0][0].split()[0]==origin['protein_id']
        assert hashlib.sha256(marker_records[0][1].encode()).hexdigest()==origin['sequence_sha256']
        with source_path.open() as h:source_records={title.split()[0]:sequence for title,sequence in SimpleFastaParser(h)}
        assert source_records[origin['protein_id']]==marker_records[0][1]

        decision_path=Path(entry['path']).with_suffix('.decisions.tsv');bind(pins,decision_path,entry['decisions_sha256'])
        bind(pins,entry['path'],entry['sha256'])
        with decision_path.open() as h: decisions={r['protein_id']:r for r in csv.DictReader(h,delimiter='\t')}
        original=decisions[origin['protein_id']]
        genes=set(json.loads(original['gene_ids_json']))
        alternatives=[r for r in decisions.values() if set(json.loads(r['gene_ids_json'])) & genes
                      and r['decision']=='longest_per_gene_lexical_tiebreak']
        assert original['decision']=='alternative_product_retained_in_source'
        assert original['status']=='unique_gene' and len(genes)==1 and len(alternatives)==1
        selected=alternatives[0]
        assert selected['status']=='unique_gene' and set(json.loads(selected['gene_ids_json']))==genes
        baseline=db.execute('SELECT sequence_sha256,length,has_afdb,has_esmfold FROM proteins WHERE taxon_id=? AND protein_id=?',(taxon,selected['protein_id'])).fetchone()
        assert baseline is not None and baseline[1]==int(selected['protein_length'])
        assert hashlib.sha256(source_records[selected['protein_id']].encode()).hexdigest()==baseline[0]
        assert len(source_records[selected['protein_id']])==baseline[1]
        assert int(selected['protein_length'])>int(original['protein_length'])
        assert db.execute('SELECT * FROM proteins WHERE taxon_id=? AND protein_id=?',(taxon,origin['protein_id'])).fetchone() is None
        model=db.execute("SELECT model_id,path,length,prediction_config_sha256 FROM models WHERE source='ESMFold' AND sequence_sha256=?",(origin['sequence_sha256'],)).fetchone()
        assert model is not None and model[2]==int(original['protein_length'])
        diagnostics.append(dict(original_marker_identity=origin,full_source_proteome_sha256=entry['source_sha256'],single_marker_fasta_sha256=origin['source_fasta_sha256'],gene_ids=sorted(genes),
            original_gene_product_decision=original,selected_gene_product_decision=selected,
            baseline_sequence_sha256=baseline[0],baseline_has_afdb=bool(baseline[2]),baseline_has_esmfold=bool(baseline[3]),
            retained_esmfold_model_id=model[0],retained_esmfold_path=model[1],prediction_config_sha256=model[3],
            disposition='annotated_alternative_product_excluded_by_fixed_longest_product_baseline',
            biological_isoform_expression_validated=False,gene_duplication_inferred=False))
    db.close()
    for path in [plan_path,producer_path,transport_path,marker_path,Path(plan['representatives']),Path(__file__)]:bind(pins,path)
    verify(pins)
    result=dict(status='completed_all_source_unlinked_esmfold_gene_product_diagnostic_pending_union_readback',
        checked_utc=datetime.now(timezone.utc).isoformat(),unlinked_model_sequences=len(unused),
        original_marker_links=len(origins),diagnostics=diagnostics,source_hashes=pins,
        scientific_eligibility=False,whole_union_independent_readback_complete=False,
        baseline_policy_changed=False,models_discarded=False,new_predictions=0,gpu=False,
        scope='Every producer-reported source-unlinked ESMFold sequence traced to original marker identity '
              'and frozen gene-product decisions. A longer annotated product explains fixed baseline '
              'exclusion; alternative model/source annotation retained. Neither true isoform expression '
              'nor duplication/annotation accuracy inferred. Full independent union readback remains required.')
    with a.receipt.open('x') as h:json.dump(result,h,indent=2,allow_nan=False);h.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
