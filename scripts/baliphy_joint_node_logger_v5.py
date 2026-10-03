"""Reversible same-record logger with native CJSON numeric property encoding.

Compose the preserved V3 full-node/sequence edits with two pure rendering edits.
The model, distributions, initializer and shared ancestral draw stay byte-identical
after reversal. The retained frame validator's normalization assertion is unchanged.
"""
from baliphy_joint_node_logger_v3 import (
    transform as sequence_transform, restore as sequence_restore,
    validate_frame, fresh_seeds,
)

PROPERTY_BEFORE = '((J.toJSONKey "properties") .= (prop_smodel_properties properties))'
PROPERTY_AFTER = '(E.pair (J.toJSONKey "properties") (J.cjsonToEncoding (J.toCJSON (prop_smodel_properties properties))))'
IMPORT_BEFORE = 'import qualified Data.JSON as J'
IMPORT_AFTER = IMPORT_BEFORE + '\nimport qualified Data.JSON.Encoding as E'
NUMERIC_EDITS = [(IMPORT_BEFORE, IMPORT_AFTER), (PROPERTY_BEFORE, PROPERTY_AFTER)]


def numeric_transform(source):
    original = source
    for before, after in NUMERIC_EDITS:
        assert source.count(before) == 1 and after not in source
        source = source.replace(before, after)
    assert numeric_restore(source) == original
    return source


def numeric_restore(source):
    for before, after in reversed(NUMERIC_EDITS):
        assert source.count(after) == 1
        source = source.replace(after, before)
    return source


def transform(source):
    result = numeric_transform(sequence_transform(source))
    assert restore(result) == source
    return result


def restore(source):
    return sequence_restore(numeric_restore(source))
