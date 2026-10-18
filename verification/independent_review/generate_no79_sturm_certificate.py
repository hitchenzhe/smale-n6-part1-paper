"""One-time exact-data generator for the No.79 Sturm certificate.

Requires SymPy 1.14 for polynomial construction only. The delivered
verify_paper_claims.py will replay the resulting coefficient lists using
fractions and integer polynomial arithmetic without SymPy.
"""
from __future__ import annotations
import json
from pathlib import Path
from math import lcm,gcd
import sympy as sp

th=sp.symbols('theta')
q=(1+th)/2
q0=sp.Rational(24,25)
C0=sp.Rational(3125,108)
bound=sp.Rational(1024,1)/(q0**2*C0**2)

def N1(t):
    return sp.expand((1+q*th*t)**2+q**2*(3*th**2*t**2+6*t))

def primitive(expr, scale):
    p=sp.Poly(sp.cancel(scale*expr),th,domain=sp.QQ)
    _,prim=p.clear_denoms(convert=True)
    content,prim=prim.primitive()
    return prim

cases={
    'D_1_5': sp.Rational(125000)*
             (N1(sp.Rational(1,5))**3*(1-th)**4*(1+th)**2-bound),
    'D_ell': sp.Rational(466560000000000000)*
             (N1(th/3+sp.Rational(17,200))**3*(1-th)**4*(1+th)**2-bound),
}
out={}
for name,expr in cases.items():
    P=primitive(expr,1)
    if name=='D_1_5': points=[sp.Rational(0),sp.Rational(3,8),sp.Rational(2,5)]
    else: points=[sp.Rational(3,8),sp.Rational(9,20)]
    chain=sp.sturm(P.as_expr(),th)
    def to_integer(poly):
        coeffs=sp.Poly(poly,th,domain=sp.QQ).all_coeffs()
        den=lcm(*[x.q for x in coeffs])
        vals=[int(x*den) for x in coeffs]
        g=0
        for x in vals: g=gcd(g,abs(x))
        return [x//g for x in reversed(vals)]
    chain_coeffs=[to_integer(c) for c in chain]
    signs=[]; variations=[]
    for x in points:
        sg=[]
        for c in chain:
            val=sp.sign(c.subs(th,x))
            if val: sg.append(int(val))
        signs.append(int(sp.sign(P.as_expr().subs(th,x))))
        variations.append(sum(a!=b for a,b in zip(sg,sg[1:])))
    out[name]={
        'primitive_polynomial_ascending':to_integer(P.as_expr()),
        'sturm_chain_ascending':chain_coeffs,
        'endpoint_rationals':[[int(x.p),int(x.q)] for x in points],
        'polynomial_endpoint_signs':signs,
        'variations':variations,
    }
    print(name,'degree',P.degree(),'endpoint signs',signs,'variations',variations)

path=Path(__file__).with_name('no79_sturm_certificate.json')
path.write_text(json.dumps(out,indent=2),encoding='utf-8')
print(path)
