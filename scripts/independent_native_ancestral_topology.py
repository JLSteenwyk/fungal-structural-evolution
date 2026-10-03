"""Manual rooted Newick topology with bitset clades and exact decimal lengths."""
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import math

from independent_native_ancestral_alignment import newick_tokens


def parse_tree(text):
    tokens=list(newick_tokens(text));index=0;stack=[];nodes=[];root=None;expect_child=True
    while index<len(tokens):
        token=tokens[index]
        if expect_child:
            if token=='(':
                stack.append([]);index+=1;continue
            if not isinstance(token,tuple):raise ValueError('Expected named subtree')
            children=[];label=token[1];index+=1
        else:
            if token==',':
                if not stack:raise ValueError('Comma outside clade')
                index+=1;expect_child=True;continue
            if token==')':
                if not stack or not stack[-1]:raise ValueError('Empty clade')
                children=stack.pop();index+=1
                if index>=len(tokens) or not isinstance(tokens[index],tuple):raise ValueError('Named node required')
                label=tokens[index][1];index+=1
            elif token==';':
                if stack or root is None or index!=len(tokens)-1:raise ValueError('Incomplete or extra tree')
                break
            else:raise ValueError('Unexpected tree suffix')
        if not label:raise ValueError('Empty node label')
        length=None
        if index<len(tokens) and tokens[index]==':':
            index+=1
            if index>=len(tokens) or not isinstance(tokens[index],tuple):raise ValueError('Missing edge length')
            try:length=Decimal(tokens[index][1])
            except InvalidOperation as e:raise ValueError('Invalid edge length') from e
            if not length.is_finite() or not math.isfinite(float(length)):raise ValueError('Nonfinite edge length')
            index+=1
        node=dict(label=label,children=children,length=length);nodes.append(node)
        if stack:stack[-1].append(node)
        else:
            if root is not None:raise ValueError('Multiple roots')
            root=node
        expect_child=False
    else:raise ValueError('Tree terminator required')
    if len({n['label'] for n in nodes})!=len(nodes):raise ValueError('Repeated node label')
    return dict(root=root,nodes=nodes,tips=sorted(n['label'] for n in nodes if not n['children']))


def clade_index(tree, tips):
    if tree['tips']!=tips:raise ValueError('Different tip sets')
    bits={tip:1<<i for i,tip in enumerate(tips)};index={};by_label={}
    # Parser appends children before parents, so this is a postorder traversal.
    for node in tree['nodes']:
        if node['children']:
            mask=0
            for child in node['children']:
                if mask & child['mask']:raise ValueError('Overlapping child clades')
                mask|=child['mask']
        else:mask=bits[node['label']]
        if mask in index:raise ValueError('Ambiguous duplicate descendant set, including unary clades')
        node['mask']=mask;index[mask]=node;by_label[node['label']]=node
    if tree['root']['mask']!=(1<<len(tips))-1:raise ValueError('Root does not span tips')
    return index,by_label


def match_trees(source_text,runtime_text):
    source,runtime=map(parse_tree,[source_text,runtime_text]);tips=source['tips']
    left,source_labels=clade_index(source,tips);right,runtime_labels=clade_index(runtime,tips)
    if set(left)!=set(right):raise ValueError('Rooted clade sets differ')
    rows=[];maximum=Fraction(0);negative=0
    for mask in sorted(left):
        a,b=left[mask],right[mask];is_root=a is source['root']
        if is_root!=(b is runtime['root']):raise ValueError('Root correspondence differs')
        delta=None
        if not is_root:
            if a['length'] is None or b['length'] is None:raise ValueError('Nonroot edge length required')
            delta=abs(Fraction(a['length'])-Fraction(b['length']))
            if delta>=Fraction('1e-10'):raise ValueError('Branch length differs at unchanged strict tolerance')
            maximum=max(maximum,delta);negative+=(a['length']<0 or b['length']<0)
        rows.append(dict(source_node=a['label'],runtime_node=b['label'],descendant_mask_hex=hex(mask),
            retained_descendants=mask.bit_count(),is_root=is_root,
            source_length=None if a['length'] is None else str(a['length']),
            runtime_length=None if b['length'] is None else str(b['length']),
            absolute_length_difference=None if delta is None else str(delta),scientific_eligibility=False))
    return dict(tips=tips,rows=rows,source_labels=source_labels,runtime_labels=runtime_labels,
        source_index=left,runtime_index=right,maximum_absolute_length_difference=str(maximum),
        negative_branch_pairs_require_review=int(negative),source_root=source['root']['label'],runtime_root=runtime['root']['label'])
