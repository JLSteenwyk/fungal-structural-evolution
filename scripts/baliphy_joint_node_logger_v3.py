"""Versioned same-record ancestral alignment/category logging experiment.

Keep the separate legacy alignment log for comparison. A second alignment
representation is encoded within the category logger action so that its
correspondence to the state records can be checked, rather than assumed.
"""
from baliphy_full_node_logger import BEFORE, AFTER, fresh_seeds, validate_frame as validate_legacy
from independent_native_ancestral_alignment import fasta_records
from independent_short_sampler_outputs import digest

FIELD_BEFORE = '((J.toJSONKey "catStates") .= catStates)'
FIELD_AFTER = '(' + FIELD_BEFORE + ' <> ((J.toJSONKey "alignmentLines") .= (map Text.pack (lines (Text.unpack ancAlignment)))))'
IMPORT_BEFORE = "import qualified Data.Text.IO as T"
IMPORT_AFTER = "import qualified Data.Text as Text\n" + IMPORT_BEFORE
EDITS = [(BEFORE, AFTER), (FIELD_BEFORE, FIELD_AFTER), (IMPORT_BEFORE, IMPORT_AFTER)]


def restore(source):
    for before, after in reversed(EDITS):
        assert source.count(after) == 1
        source = source.replace(after, before)
    return source


def transform(source):
    original = source
    for before, after in EDITS:
        assert source.count(before) == 1 and after not in source
        source = source.replace(before, after)
    assert restore(source) == original
    assert source.count(';let {ancStates  = prop_anc_cat_states properties}') == 1
    assert source.count('observe sequenceData (phyloCTMC tree alignment smodel scale1)') == 1
    return source


def validate_frame(frame, iteration, sequences, tips):
    assert set(frame) == {'iter', 'catStates', 'alignmentLines', 'properties', 'conditions'}
    assert type(frame['alignmentLines']) is list and frame['alignmentLines']
    assert all(type(line) is str and line and not any(ord(c) < 32 or c in '\"\\' for c in line) for line in frame['alignmentLines'])
    encoded = fasta_records(frame['alignmentLines'])
    assert encoded == sequences
    legacy = {k: v for k, v in frame.items() if k != 'alignmentLines'}
    result = validate_legacy(legacy, iteration, sequences, tips)
    result['frame_sha256'] = digest(frame)
    result['same_record_alignment_state_correspondence_checked'] = True
    return result
