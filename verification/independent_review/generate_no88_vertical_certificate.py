"""Exact nine-sheet No. 88 vertical-divisor check over Z[x,y,z,omega].

At t=1 the two recovery cubics give v=omega^a,w=omega^b. Differentiate
the CM numerator along each sheet, clear only (1-x^3)(1-y^3), then
take the quadratic-field norms. A constant gcd rules out a common divisor
of the value and derivative. Five conjugacy representatives cover all nine.
"""
from __future__ import annotations
import json
from pathlib import Path
import sympy as s
from exact_geometry import cm_inverse_root_numerator

x,y,z,t,v,w,o=s.symbols('x y z t v w omega')


def pair(expr):
    p=s.Poly(expr,o)
    r=p.rem(s.Poly(o**2+o+1,o))
    return s.Poly(r.nth(0),x,y,z,domain=s.ZZ),s.Poly(r.nth(1),x,y,z,domain=s.ZZ)


def norm(a,b):return a*a-a*b+b*b


def terms(poly):
    return [[*m,int(c)] for m,c in poly.terms()]


def main():
    C=cm_inverse_root_numerator((x,y,z,t,v,w)).as_expr()
    px=1-x**3;py=1-y**3;pz=1-z**3
    # v'(1)=v*pz/py and w'(1)=w*pz/px, since v^3=w^3=1.
    derivative=px*py*s.diff(C,t)+px*pz*v*s.diff(C,v)+py*pz*w*s.diff(C,w)
    records=[]
    reps=((0,0),(0,1),(1,0),(1,1),(1,2))
    for a,b in reps:
        sub={t:1,v:o**a,w:o**b}
        A=pair(C.subs(sub));B=pair(derivative.subs(sub))
        na=norm(*A);nb=norm(*B)
        g=s.gcd(na,nb)
        record={'phase_exponents':[a,b],
                'conjugate_phase_exponents':[(-a)%3,(-b)%3],
                'value_norm_terms_xyz':terms(na),
                'derivative_norm_terms_xyz':terms(nb),
                'gcd':str(g.as_expr()),'gcd_degree':int(g.total_degree()),
                'value_norm_total_degree':int(na.total_degree()),
                'derivative_norm_total_degree':int(nb.total_degree())}
        records.append(record)
        print(f'phases {(a,b)}: norm degrees {na.total_degree()},{nb.total_degree()}; gcd={g.as_expr()}',flush=True)
    result={'diagram':88,'source':'CM determinant plus differentiated monic cubic recoveries',
            'coefficient_ring':'ZZ[x,y,z,omega]/(omega^2+omega+1)',
            'norm_definition':'N(a+b*omega)=a^2-a*b+b^2',
            'derivative_clearing_factor':'(1-x^3)*(1-y^3)',
            'phase_representatives':5,'all_nine_sheets_covered':True,
            'all_norm_gcds_constant':all(r['gcd_degree']==0 for r in records),
            'records':records}
    covered=set()
    for r in records:covered.add(tuple(r['phase_exponents']));covered.add(tuple(r['conjugate_phase_exponents']))
    assert len(covered)==9
    result['all_nine_sheets_covered']=len(covered)==9
    out=Path(__file__).with_name('no88_vertical_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f'Saved {out}; constant norm gcds: {result["all_norm_gcds_constant"]}',flush=True)


if __name__=='__main__':main()
