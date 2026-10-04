"""Independent strict decoder and mapped TSV check for project scalar V6.

This module does not import the generator or its quality validator. Nonfinite
tags retain diagnostic states; they never admit a sample as an adequate posterior.
"""
import csv
import json
import math
from pathlib import Path


def load(text):
    def object_pairs(rows):
        out = {}
        for key, value in rows:
            if key in out:raise ValueError('Duplicate key')
            out[key] = value
        return out
    def invalid(value):raise ValueError('Nonstandard JSON number '+value)
    return json.loads(text, object_pairs_hook=object_pairs, parse_constant=invalid)


def read_record(row):
    if type(row) is not dict or sorted(row) != ['iter', 'numericParameterQuality//', 'parameters//', 'statistics//']:
        raise ValueError('Scalar V6 schema required')
    if type(row['iter']) is not int or row['iter'] < 0:raise ValueError('Invalid iteration')
    contexts = row['statistics//']
    if type(contexts) is not dict or '__project_scalar_v6_quality__' not in contexts:
        raise ValueError('Context quality absent')
    context = {k:v for k,v in contexts.items() if k != '__project_scalar_v6_quality__'}
    sections = [('statistics', context, contexts['__project_scalar_v6_quality__']),
                ('parameters', row['parameters//'], row['numericParameterQuality//'])]
    reviews = [];total = 0
    for section, data, audit in sections:
        if type(data) is not dict or type(audit) is not dict or set(audit) != {'numericLeafCount','nonfinite','literalNullPaths'}:
            raise ValueError('Wrong quality schema')
        if type(audit['numericLeafCount']) is not int or audit['numericLeafCount'] < 0:
            raise ValueError('Wrong count')
        stack = [((), data)];nulls = [];numeric = 0;placeholders = {}
        while stack:
            path, value = stack.pop()
            if isinstance(value, dict):
                for key, child in value.items():
                    if type(key) is not str:raise ValueError('Nonstring key')
                    stack.append((path+(key,), child))
            elif isinstance(value, list):
                stack.extend((path+(i,), child) for i,child in enumerate(value))
            elif value is None:raise ValueError('Unrepresented null')
            elif type(value) is str and value.startswith('__project_scalar_v6__:'):
                nulls.append(path);placeholders[path]=value
            elif type(value) in (int,float):
                if not math.isfinite(value):raise ValueError('Unencoded nonfinite')
                numeric += 1
            elif type(value) not in (str,bool):raise ValueError('Wrong leaf')
        if type(audit['nonfinite']) is not list or type(audit['literalNullPaths']) is not list:
            raise ValueError('Nonlist audit')
        def normalize(p):
            if type(p) is not list or not p or any(type(x) not in (str,int) or (type(x) is int and x<0) for x in p):
                raise ValueError('Bad typed path')
            return tuple(p)
        annotated = []
        for tag in audit['nonfinite']:
            if type(tag) is not dict or set(tag) != {'kind','path'} or tag['kind'] not in ('nan','positive_infinity','negative_infinity'):
                raise ValueError('Bad nonfinite tag')
            p = normalize(tag['path']); annotated.append(p)
            reviews.append({'section':section,'path':list(p),'kind':tag['kind']})
        literals = [normalize(p) for p in audit['literalNullPaths']]
        combined = annotated+literals
        if len(set(combined)) != len(combined) or set(combined) != set(nulls):
            raise ValueError('Null/tag identity or duplication error')
        for tag in audit['nonfinite']:
            if placeholders[tuple(tag['path'])]!='__project_scalar_v6__:'+tag['kind']:
                raise ValueError('Wrong placeholder kind')
        if any(placeholders[p]!='__project_scalar_v6__:literal_null' for p in literals):
            raise ValueError('Wrong literal-null placeholder')
        if numeric+len(annotated) != audit['numericLeafCount']:raise ValueError('Numeric census mismatch')
        total += audit['numericLeafCount']
    return dict(iteration=row['iter'], numeric_leaves=total, nonfinite_reviews=reviews,
                finite_record=not reviews, scientific_eligibility=False, posterior_qualified=False)


def compare_tsv(directory):
    directory = Path(directory)
    lines = (directory/'C1.log.json').read_text().splitlines();header = load(lines[0])
    if header != {'fields':['iter','prior','likelihood','posterior'],'nested':True,'format':'MCON',
                  'version':'0.2','projectScalarSchema':'native-cjson-explicit-special-values-v6'}:
        raise ValueError('Wrong header')
    records = [load(l) for l in lines[1:]]
    with (directory/'C1.log').open() as handle:tsv = list(csv.DictReader(handle,delimiter='\t'))
    mapping = load((directory/'C1.log.column-map.json').read_text())
    if len(tsv) != len(records):raise ValueError('Row census mismatch')
    compared = 0;reviews = []
    for row, reference in zip(records,tsv):
        audit = read_record(row); reviews.extend(audit['nonfinite_reviews'])
        flat = {'iter':row['iter']}
        def visit(value,prefix):
            for key,item in value.items():
                if isinstance(item,dict):
                    if key.endswith('/'):visit(item,prefix+key)
                    else:
                        for subkey,scalar in item.items():flat[prefix+key+'['+subkey+']']=scalar
                else:flat[prefix+key]=item
        visit({k:v for k,v in row['statistics//'].items() if k!='__project_scalar_v6_quality__'},'')
        visit(row['parameters//'],'')
        if set(flat) != set(mapping):raise ValueError('Mapped field census mismatch')
        if audit['nonfinite_reviews']:raise ValueError('Diagnostic nonfinite row requires separate review')
        for field, column in mapping.items():
            expected = float(reference[column]);value = flat[field]
            if type(value) not in (float,int) or not math.isfinite(expected) or not math.isclose(value,expected,rel_tol=2e-13,abs_tol=0):
                raise ValueError('Mapped scalar discrepancy: '+field)
            compared += 1
    return dict(rows=len(records), mapped_values_compared=compared, nonfinite_reviews=reviews,
                relative_tolerance=2e-13, absolute_tolerance=0.0, scientific_eligibility=False)
