#!/usr/bin/env python3
"""Write concrete schemas for every expanded case, directed state and contrast field."""
import argparse
import csv
import gzip
from pathlib import Path
from full_matching_case_sources import CASE_FIELDS
from full_expanded_measurement_catalog_sources_v3 import FIELDS as DIRECTED
from full_expanded_case_measurement_sources import FIELDS as CONTRAST,STATES,STATE_FIELDS,OUTCOMES


def definition(artifact,field):
    ids={'row_identity':'SHA256 of namespace, logical case ID and physical confidence mask; future model/covariance row identity.',
        'case_id':'SHA256 of namespace and ordered target/background node IDs; gene contexts remain distinct.',
        'physical_case_id':'SHA256 of namespace and ordered target/background model-pair IDs; reuse descriptor, not an independent sample.',
        'target_id':'Original fixed-matching target node ID; joins target_nodes.jsonl.',
        'background_id':'Original fixed-matching background node ID; joins background_nodes.jsonl.',
        'pair_key':'SHA256 of the sorted pair of model-ID/version endpoints.',
        'target_pair_key':'Original physical target model-pair ID; role remains target.',
        'background_pair_key':'Original physical background model-pair ID; role remains background.',
        'guide':'Original gene-tree guide: profile or mafft.',
        'gene_node':'Original target gene-tree node; not an accepted reconciled duplication assignment.',
        'focal_taxon':'Original target taxon entry ID.',
        'role':'target or background; shared physical pairs across roles are kept separate.',
        'mask':'Physical coordinate mask: full or plddt70; both is eligibility metadata, not a measurement row.',
        'source_kind':'Frozen native owner label: target_original, background_old, reference_old or background_new_native.',
        'source_checkpoint':'Exact original native checkpoint path, bound in the complete source proof.',
        'source_checkpoint_sha256':'Exact original checkpoint SHA256; verified through its original receipt lineage.',
        'native_status':'Original native disposition, including failed/unavailable states.',
        'numerical_exclusion_reasons':'Ordered semicolon-separated exclusions; original native failure label or RMSD/short/rotation reasons.',
        'rmsd_status':'Original independently checked printed-rounding classification.',
        'geometry_status':'Original independently checked numeric rotation-uniqueness classification.'}
    if field in ids:return 'string','none',ids[field]
    if field in ['policy','scenario_id','family','target_comparison_disposition','background_comparison_disposition']:
        return 'string','none',{'policy':'Original one of four fixed architecture support policies; not selected on structural outcomes.', 'scenario_id':'Original frozen metadata caliper/scoring scenario:S01 throughS54.', 'family':'Original target gene-family ID.', 'target_comparison_disposition':'Original distinct_model_pair or identical_model_no_alignment; not structural conservation.', 'background_comparison_disposition':'Original distinct_model_pair or identical_model_no_alignment; not structural conservation.'}[field]
    if field in ['source_row_ordinal','endpoint_order','eligible_candidates','equal_score_candidates']:
        return 'integer','count or label',{'source_row_ordinal':'One-based ordinal in the entire original selected_pair_coverage stream; exact source association retained.', 'endpoint_order':'Original fixed-matching endpoint map label:0 or1; distinct from native alignment order.', 'eligible_candidates':'Original number of metadata-eligible controls in this scenario.', 'equal_score_candidates':'Original count tied at the selected metadata matching score.'}[field]
    if field in ['score','score_gap_to_second']:return 'float','metadata matching score','Original frozen selection score or gap to runner-up; no structural-outcome rematching. Gap can be blank when no runner-up exists.'
    if field in ['order','source_order']:
        return 'integer','index',('Canonical input direction:0=A-left/B-right,1=B-left/A-right.' if field=='order' else 'Original native checkpoint order label; may differ from canonical order.')
    if field=='numerical_usable':return 'integer','boolean','1 requires aligned native status, printed-RMSD agreement, at least three pairs and unique numeric rotation;0 quarantines raw diagnostics.'
    if field.endswith('_same_model'):return 'integer','boolean','1 means identical model ID/version endpoints; no native measurement or invented zero distance.'
    if field.endswith('_sequence_distance'):return 'float','gene-tree distance','Original fixed-matching sequence divergence; zero remains zero, no log epsilon or time interpretation.'
    if field.endswith('_pass_bits'):
        return 'integer','bit mask','Six original-length screen bits, low to high:n30_c50,n30_c70,n30_c90,n50_c50,n50_c70,n50_c90. Both-mask bits are intersections. Pair-mask bits require both numeric orders.'
    if field.endswith('_family'):return 'string','none','Original gene-family ID in the corresponding graph node; not a structural cluster.'
    if field.startswith('target_gene_') or field.startswith('background_gene_'):return 'string','none','Original annotated gene ID, including logical genes that share a predicted model.'
    if field.startswith('background_taxon_'):return 'string','none','Original background endpoint taxon entry ID.'
    if field in ['selection_records','endpoint_order_bits','policy_bits','scenario_bits']:
        return 'integer','count or bit mask',{'selection_records':'Number of original selections using this logical case; dependent reuse, not sample size.', 'endpoint_order_bits':'Original matching endpoint mapping labels observed:bit0=0,bit1=1; separate from native alignment order.', 'policy_bits':'Bits in plan policy order:alignment_evalue,alignment_bitscore,envelope_evalue,envelope_bitscore.', 'scenario_bits':'Bits in original54-scenario order; bit53 is S54. Full membership table retains exact policy/scenario associations.'}[field]
    if field in ['model_a','model_b']:return 'string','none','Canonical sorted versioned-model endpoint ID; nativeleft/right direction follows order.'
    if field.startswith('version_'):return 'integer','version','Original AlphaFold model version for this endpoint.'
    if field.startswith('original_length_'):return 'integer','residues','Original pre-mask model protein length from the fixed coverage source; distinct from masked alignment-input length.'
    if field=='aligned_length':return 'integer','residue pairs','Original native aligned-pair count; blank for unavailable native states. Raw excluded counts stay quarantined in the directed catalog.'
    if field=='joint_plddt70_pairs':return 'integer','residue pairs','Number of aligned pairs where both original CA confidence values pass70.'
    if field in ['rmsd_recomputed','rmsd_native','rmsd_rounding_error']:return 'float','angstrom',{'rmsd_recomputed':'Independently reconstructed CA RMSD; raw excluded diagnostics require numerical_usable gating.', 'rmsd_native':'Original printed native RMSD; not substituted for independently checked RMSD.', 'rmsd_rounding_error':'Absolute independently recomputed-versus-printed native RMSD difference.'}[field]
    if field in ['sequence_identity_exact','joint_plddt70_fraction']:return 'float','fraction','Original aligned exact sequence identity or joint confidence fraction; not whole-gene identity or prediction accuracy.'
    if field in ['tm_left_native','tm_right_native']:return 'float','native normalized score','Original native TM score normalized by the directed native endpoint length; not independently optimized TM.'
    if field in ['coverage_left','coverage_right']:return 'float','fraction','Aligned count/native alignment-input length; mask-dependent denominator, not original complete length.'
    if field.startswith('original_coverage_'):return 'float','fraction','Raw aligned count/original pre-mask model length; raw excluded values remain quarantined, independent of target QC blanking.'
    if field in ['rank_left','rank_right']:return 'integer','rank','Original numeric rank of the centered native-directed coordinate point set.'
    if field=='determinant_correction':return 'integer','sign','Original proper-rotation determinant correction,+1 or-1.'
    if field.startswith('rms_width'):return 'float','angstrom','Original singular coordinate width for the directed native endpoint.'
    if field.startswith('width2_to_width1') or field in ['relative_rotation_curvature','relative_numeric_tolerance']:return 'float','dimensionless','Original geometry shape/curvature/tolerance diagnostic; not prediction uncertainty.'
    if field.startswith('cross_s') or field=='minimum_rotation_curvature':return 'float','angstrom squared','Original centered cross-product singular/rotation-curvature diagnostic.'
    for role,order in STATES:
        prefix=f'{role}_order{order}_'
        if field.startswith(prefix):
            name=field[len(prefix):]
            if name=='state_key':return 'string','foreign key','role|pair_key|mask|canonical-order in the closed directed catalog; blank for same-model cases.'
            if name in OUTCOMES:return 'float','angstrom' if name=='rmsd' else 'native dimensionless','Usable per-order outcome only; raw excluded values remain in the directed catalog. Native TM dissimilarity=1-(TM_left+TM_right)/2.'
            if name in STATE_FIELDS:return definition('directed_measurements.tsv.gz',name)
    for outcome in OUTCOMES:
        prefix=outcome+'_'
        if field.startswith(prefix):
            name=field[len(prefix):]
            if name=='disposition':return 'string','none','complete_four_order_pairs, incomplete_numerical_orders or identical_model_no_alignment; quality gates remain separate.'
            if name=='usable_order_pair_bits':return 'integer','bit mask','Valid target/control order cells:bit0=00,bit1=01,bit2=10,bit3=11.15 required for complete envelope.'
            if name.startswith('order_pair_'):meaning='Target-minus-background difference for the explicit cross-order pair; null if either numeric state is unavailable/excluded.'
            elif name.startswith('complete_order_envelope_'):meaning='Min/max/span across allfour valid cross-order deltas; null unless allfour exist. Sensitivity envelope, not confidence interval.'
            elif name.endswith('_both_orders_mean'):meaning='Arithmetic mean of both usable orders on the named side; never chooses an order or substitutes a one-order mean.'
            else:meaning='Target both-orders mean minus background both-orders mean; null unless both sides have complete numeric means.'
            return 'float','angstrom' if outcome=='rmsd' else 'native dimensionless',meaning
    raise ValueError((artifact,field))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--membership-source',type=Path,default=Path('results/structural_comparisons/full-matching-logical-cases-20261002-v1/selection_case_links.tsv.gz'));a=p.parse_args()
    with gzip.open(a.membership_source,'rt') as handle:membership=csv.DictReader(handle,delimiter='\t').fieldnames
    with a.output.open('x') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n');writer.writerow(['artifact','field','type','unit','definition','missing_value'])
        for artifact,fields in [('case_index.tsv.gz',CASE_FIELDS),('selection_case_links.tsv.gz',membership),('directed_measurements.tsv.gz',DIRECTED),('case_mask_contrasts.tsv.gz',CONTRAST)]:
            assert len(fields)==len(set(fields))
            for field in fields:
                kind,unit,text=definition(artifact,field);writer.writerow([artifact,field,kind,unit,text,'Blank numeric=unavailable/excluded;zero is a real numeric value. Original status/flags distinguish states.' if kind=='float' else 'See source/status; no silent replacement.'])


if __name__=='__main__':main()
