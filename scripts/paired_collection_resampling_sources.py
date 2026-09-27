#!/usr/bin/env python3
"""Resolve and verify native fits for resampling a provenance-preserving collection."""
import csv,json,shutil
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha
from combine_audited_paired_fit_sources import checked_receipt_for_label
from readback_combined_paired_fit_sources import verify_collection


def load_sources(inputs, models, collection, inventory_path):
    inputs,models,collection,inventory_path=map(Path,(inputs,models,collection,inventory_path))
    verify_collection(collection,inventory_path)
    inventory=json.loads(inventory_path.read_text())
    if inventory['status']!='passed_paired_fit_source_input_and_settings_compatibility':raise ValueError('Unverified source inventory')
    for p,h in inventory['source_sha256'].items():
        if sha(Path(p))!=h:raise ValueError('Changed inventory source')
    pinned_inputs={Path(p).resolve():h for p,h in inventory['source_sha256'].items()}
    if pinned_inputs.get((inputs/'receipt.json').resolve())!=sha(inputs/'receipt.json'):raise ValueError('Wrong expanded cohort')
    checked_receipt(inputs);checked_receipt(models)
    markers=[r for r in csv.DictReader((inputs/'marker_summary.tsv').open(),delimiter='\t') if r['status']=='ready_for_inference']
    selected={r['marker']:r for r in inventory['sources']}
    if len(selected)!=len(inventory['sources']) or set(selected)!={r['marker'] for r in markers}:raise ValueError('Marker scope differs')
    executable=shutil.which('iqtree3')
    if not executable:raise ValueError('IQ-TREE unavailable')
    configurations={};native={};bindings={}
    for marker,s in selected.items():
        folder=Path(s['fit_directory']).resolve();root=folder.parent
        if root not in configurations:
            config=json.loads((root/'config.json').read_text());receipt=json.loads((root/'receipt.json').read_text())
            if receipt['status']!='complete_matched_topology_point_estimates' or receipt['config_sha256']!=sha(root/'config.json'):raise ValueError('Incomplete native fits')
            if config['model_receipt_sha256']!=sha(models/'receipt.json') or config['executable_sha256']!=sha(Path(executable)):raise ValueError('Different model/executable')
            configurations[root]=config
        for name,digest in s['expanded_input_sha256'].items():
            if sha(inputs/marker/name)!=digest or sha(Path(s['input_directory'])/name)!=digest:raise ValueError('Selected alignment changed')
        config=configurations[root]
        if config['input_receipt_sha256']!=sha(Path(s['input_directory']).parent/'receipt.json'):raise ValueError('Native input receipt differs')
        for label,alignment in [('aa','aa.faa'),('3di_af','3di.faa')]:
            fit=checked_receipt_for_label(folder,label)
            fcpath=folder/(label+'.config.json');fc=json.loads(fcpath.read_text())
            if fit['config_sha256']!=sha(fcpath) or fc['parent_config_sha256']!=sha(root/'config.json') or fc['alignment_sha256']!=s['expanded_input_sha256'][alignment]:raise ValueError('Selected fit provenance differs')
            if label!='aa' and fc['topology_sha256']!=sha(folder/'aa.treefile'):raise ValueError('Selected topology differs')
        native[marker]=folder
        bindings[marker]=dict(source=s['source'],fit_directory=str(folder),native_receipt_sha256=sha(root/'receipt.json'),native_config_sha256=sha(root/'config.json'),aa_receipt_sha256=sha(folder/'aa.receipt.json'),structural_receipt_sha256=sha(folder/'3di_af.receipt.json'),topology_sha256=sha(folder/'aa.treefile'))
    return markers,native,executable,bindings
