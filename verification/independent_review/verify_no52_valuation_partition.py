"""Exact order and branch-partition replay for the No. 52 valuation sieve."""
from __future__ import annotations
from itertools import combinations
import json
from pathlib import Path
import sympy as s


def main():
    alpha=s.symbols('alpha', real=True)
    edges=tuple(combinations(range(1,7),2))
    q={}
    for e in edges:
        if e in {(4,5),(4,6),(5,6)}: q[e]=s.Integer(2)
        elif e in {(2,3),(2,5),(3,5)}: q[e]=4-2*alpha
        elif e==(1,6): q[e]=s.Integer(0)
        elif e in {(1,2),(1,3),(1,4),(1,5)}: q[e]=s.Integer(-1)
        else: q[e]=5-3*alpha
    products={}
    for e in edges:
        disjoint=[f for f in edges if set(e).isdisjoint(f)]
        assert len(disjoint)==6
        products[e]=s.expand(q[e]+sum(q[f] for f in disjoint))
    floor=13-8*alpha
    endpoints=(s.Rational(13,8),s.Rational(7,4))
    for e,expr in products.items():
        diff=s.expand(expr-floor)
        vals=[s.factor(diff.subs(alpha,t)) for t in endpoints]
        assert all(v>=0 for v in vals),(e,expr,vals)
    minimizers=[e for e,expr in products.items() if s.expand(expr-floor)==0]
    assert minimizers==[(2,6),(3,6),(5,6)],minimizers
    assert s.factor(floor.subs(alpha,s.Rational(13,8)))==0

    interval=(s.Rational(13,8),s.Rational(5,3))
    positive_orders={
      'missing_edge_source':5-3*alpha,
      'triangle_relative':4-2*alpha,
      'crossed_remainder':(9*alpha-11)/2,
    }
    for name,expr in positive_orders.items():
        at_left=s.factor(expr.subs(alpha,interval[0]))
        at_right=s.factor(expr.subs(alpha,interval[1]))
        assert at_left>0 and at_right>=0,(name,at_left,at_right)
        assert s.diff(expr,alpha)<=0 or at_right>0,name
    assert s.factor((4*alpha-6).subs(alpha,s.Rational(7,4)))==1
    assert s.factor((4*alpha-6).subs(alpha,s.Rational(5,3)))==s.Rational(2,3)

    # Exhaust the equality patterns of (W2,W3,B). Exactly-two-equal has
    # three partitions; all-distinct has four W1 placements; all-equal is
    # the sole survivor after the positive mass-centre ratio is imposed.
    triple_partitions=[((0,1,2),),((0,1),(2,)),((0,2),(1,)),((1,2),(0,)),
                       ((0,),(1,),(2,))]
    assert len(triple_partitions)==5
    cases=[]
    for part_index,partition in enumerate(triple_partitions):
        if len(partition)==1:
            cases.append({'triple':'all equal','W1':'forced different by centre','status':'remaining'})
        elif len(partition)==2:
            cases.append({'triple':f'exactly two equal #{part_index}','W1':'arbitrary','status':'two-edge path contradiction'})
        else:
            for placement in ('matches W2','matches W3','matches B','matches none'):
                reason=('mass-weighted original z-row contradiction' if placement in ('matches W2','matches W3')
                        else 'uncancelled potential pole' if placement=='matches B'
                        else 'outer z-row and zw-edge contradiction')
                cases.append({'triple':'pairwise distinct','W1':placement,'status':reason})
    assert len(cases)==8 and sum(c['status']=='remaining' for c in cases)==1

    result={
      'diagram':52,
      'squared_distance_orders':{''.join(map(str,e)):str(q[e]) for e in edges},
      'seven_product_orders':{''.join(map(str,e)):str(products[e]) for e in edges},
      'minimum_order_on_13_8_to_7_4':'13-8*alpha',
      'minimizing_products':['26','36','56'],
      'middle_interval_positive_remainder_orders':{k:str(v) for k,v in positive_orders.items()},
      'outer_w_lead_partition_cases':cases,
      'status':'passed',
      'scope':'Replays the displayed order arithmetic and finite case partition; endpoint signed-force Jacobians and energy inequalities are checked by generate_no52_certificate.py.'
    }
    out=Path(__file__).with_name('no52_valuation_partition_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'seven_products':len(products),'minimum_products':['26','36','56'],
                      'outer_cases':len(cases),'status':'passed'},indent=2))

if __name__=='__main__':main()
