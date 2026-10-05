"""Unmodified genomic CDS candidates; uncertain ordering is an explicit disposition."""
import hashlib
import re

from Bio.Seq import Seq


def coordinate_status(seqid, start, end, lengths, circular):
    if seqid not in lengths:
        return 'missing_genomic_sequence_id'
    if start < 1 or end < start:
        return 'invalid_coordinate_interval'
    length = lengths[seqid]
    if end <= length:
        return 'within_original_sequence_bounds'
    if seqid not in circular:
        return 'out_of_bounds_without_circular_annotation'
    if end - start + 1 > length or end > 2 * length:
        return 'circular_span_requires_review'
    return 'documented_circular_virtual_coordinates'


def reconstruct(parts, genomes, circular, lengths=None):
    """Concatenate deposited intervals without trimming phases, overlaps or exceptions."""
    assert parts
    if lengths is None:
        lengths = {key: len(value) for key, value in genomes.items()}
    statuses = [coordinate_status(p['seqid'], p['start'], p['end'], lengths, circular) for p in parts]
    report = dict(feature_row_ids=[p['row_id'] for p in parts], coordinate_statuses=statuses,
                  original_phases=[p['phase'] for p in parts], sequence_sha256=None, sequence_length=None,
                  annotation_exceptions=[dict(row_id=p['row_id'], attributes={k:v for k,v in p['attrs'].items()
                    if k in ('exception', 'transl_except', 'partial', 'start_range', 'end_range',
                             'pseudo', 'pseudogene', 'part', 'transl_table')}) for p in parts],
                  order_source=None, ordered_feature_row_ids=[], status=None, flags=[])
    admitted_coordinates = {'within_original_sequence_bounds', 'documented_circular_virtual_coordinates'}
    if any(status not in admitted_coordinates for status in statuses):
        report['status'] = 'coordinate_review_required'
        return report, None
    if any(p['strand'] not in ('+', '-') for p in parts):
        report['status'] = 'unsupported_or_unknown_cds_strand'
        return report, None
    explicit = [p['attrs'].get('part', []) for p in parts]
    if any(explicit):
        positions = []
        for value in explicit:
            match = re.fullmatch(r'([1-9][0-9]*)/([1-9][0-9]*)', value[0]) if len(value) == 1 else None
            if not match:
                report['status'] = 'incomplete_or_invalid_explicit_part_order'
                return report, None
            positions.append(tuple(map(int, match.groups())))
        if {total for _, total in positions} != {len(parts)} or sorted(n for n, _ in positions) != list(range(1, len(parts)+1)):
            report['status'] = 'conflicting_explicit_part_order'
            return report, None
        ordered = [part for _, part in sorted(zip(positions, parts), key=lambda x: x[0])]
        report['order_source'] = 'complete_explicit_part_X_over_Y'
    else:
        if len({(p['seqid'], p['strand']) for p in parts}) != 1:
            report['status'] = 'multiple_sequences_or_strands_without_explicit_order'
            return report, None
        if len(parts) > 1 and parts[0]['seqid'] in circular:
            report['status'] = 'circular_multipart_order_requires_review'
            return report, None
        ordered = sorted(parts, key=lambda p: (p['start'], p['end'], p['row_id']), reverse=parts[0]['strand'] == '-')
        report['order_source'] = 'single_sequence_strand_coordinate_order'
    ordered_positions = sorted(parts, key=lambda p: (p['seqid'], p['start'], p['end']))
    if any(a['seqid'] == b['seqid'] and a['end'] >= b['start'] for a, b in zip(ordered_positions, ordered_positions[1:])):
        report['flags'].append('overlapping_parts_concatenated_without_repair')
    if any(p['phase'] != '0' for p in ordered):
        report['flags'].append('original_nonzero_or_unknown_phase_not_trimmed')
    if any(p['attrs'].get('exception') or p['attrs'].get('transl_except') for p in ordered):
        report['flags'].append('annotation_exception_retained_not_corrected')
    chunks = []
    for part in ordered:
        sequence = genomes[part['seqid']]
        length = len(sequence)
        low, high = part['start'] - 1, part['end']
        if high <= length:
            chunk = sequence[low:high]
        else:
            low %= length
            span = part['end'] - part['start'] + 1
            first = min(span, length-low)
            chunk = sequence[low:low+first] + sequence[:span-first]
        if part['strand'] == '-':
            chunk = str(Seq(chunk).reverse_complement())
        chunks.append(chunk.upper())
    joined = ''.join(chunks)
    report.update(status='unmodified_genomic_cds_candidate', ordered_feature_row_ids=[p['row_id'] for p in ordered],
                  sequence_length=len(joined), sequence_sha256=hashlib.sha256(joined.encode('ascii')).hexdigest())
    return report, joined
