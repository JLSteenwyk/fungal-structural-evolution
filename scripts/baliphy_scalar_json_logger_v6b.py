"""Reversible V5-to-V6 scalar logger rendering with explicit nonfinite tags.

The one context action is evaluated once, and tagged before CJSON can replace
nonfinite numbers with null. Parameter values use the native encoder directly.
The TSV, ancestral loggers, model and initializer are byte-preserved by reversal.
"""
import math
from pathlib import Path

HELPER = Path(__file__).resolve().parent.parent / 'config/native-format-probes/project-scalar-json-v6b.hs'
SCHEMA = 'native-cjson-explicit-nonfinite-v6'
QUALITY_KEY = '__project_scalar_v6_quality__'
IMPORT_BEFORE = 'import qualified Data.JSON.Encoding as E'
IMPORT_AFTER = IMPORT_BEFORE + '\nimport qualified Data.JSON.Types.Internal as JSONInternal\nimport MCMC.Types (runContextAction)\nimport Compiler.RealFloat (isNaN, isInfinite, isNegativeZero)'
OPEN_BEFORE = ';logParamsJSON <- if jsonEnabled then jsonLogger '
OPEN_AFTER = ';logParamsJSON <- if jsonEnabled then projectScalarJSONLoggerV6 '


def transform(source):
    helper = HELPER.read_text()
    assert source.count(IMPORT_BEFORE) == source.count(OPEN_BEFORE) == 1
    assert 'projectScalarJSONLoggerV6' not in source and QUALITY_KEY not in source
    assert source.count('\nsample_smodel ') == 1
    result = source.replace(IMPORT_BEFORE, IMPORT_AFTER).replace(OPEN_BEFORE, OPEN_AFTER)
    result = result.replace('\nsample_smodel ', '\n'+helper+'\nsample_smodel ')
    assert restore(result) == source
    return result


def restore(source):
    helper = HELPER.read_text()
    assert source.count(helper) == source.count(IMPORT_AFTER) == source.count(OPEN_AFTER) == 1
    return source.replace('\n'+helper+'\nsample_smodel ', '\nsample_smodel ').replace(
        OPEN_AFTER, OPEN_BEFORE).replace(IMPORT_AFTER, IMPORT_BEFORE)


def _quality(values, quality):
    assert type(values) is dict and type(quality) is dict
    assert set(quality) == {'numericLeafCount', 'nonfinite', 'literalNullPaths'}
    assert type(quality['numericLeafCount']) is int and quality['numericLeafCount'] >= 0
    numeric = 0; nulls = set()
    def walk(value, path):
        nonlocal numeric
        if type(value) is dict:
            for k, v in value.items():
                assert type(k) is str
                walk(v, path+(k,))
        elif type(value) is list:
            for index, v in enumerate(value):walk(v, path+(index,))
        elif value is None:nulls.add(path)
        elif type(value) in [int, float]:
            assert math.isfinite(value);numeric += 1
        else:assert type(value) in [bool, str]
    walk(values, ())
    def path(value):
        assert type(value) is list and value and all(type(x) in [int, str] for x in value)
        assert all(type(x) is not int or x >= 0 for x in value)
        return tuple(value)
    assert type(quality['nonfinite']) is list and type(quality['literalNullPaths']) is list
    tagged = {}
    for entry in quality['nonfinite']:
        assert type(entry) is dict and set(entry) == {'path', 'kind'}
        p = path(entry['path']);assert p not in tagged
        assert entry['kind'] in ['nan', 'positive_infinity', 'negative_infinity']
        tagged[p] = entry['kind']
    literal = [path(p) for p in quality['literalNullPaths']]
    assert len(literal) == len(set(literal)) and set(literal).isdisjoint(tagged)
    assert nulls == set(literal) | set(tagged)
    assert numeric + len(tagged) == quality['numericLeafCount']
    return dict(numeric_leaves=numeric+len(tagged), nonfinite=tagged, literal_null_paths=set(literal))


def validate_record(record):
    assert type(record) is dict and set(record) == {'iter', 'statistics//', 'parameters//', 'numericParameterQuality//'}
    assert type(record['iter']) is int and record['iter'] >= 0
    statistics = dict(record['statistics//']);quality = statistics.pop(QUALITY_KEY)
    context = _quality(statistics, quality)
    parameters = _quality(record['parameters//'], record['numericParameterQuality//'])
    return dict(iteration=record['iter'], context=context, parameters=parameters,
                finite_record=not context['nonfinite'] and not parameters['nonfinite'])
