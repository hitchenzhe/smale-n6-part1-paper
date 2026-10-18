"""Rebuild the No.79 exact Hurwitz certificate from the printed sextic.

Requires SymPy 1.14. All algebra is over QQ/ZZ; no floating point is used.
"""
from __future__ import annotations
import json
from pathlib import Path
from math import comb
import sympy as sp

A,X,z,v=sp.symbols('A X z v')
U=(A+X)*z+1
V=X*z+A+1
W=A*(1-z)
Q=A**2*U*W-A**2*X**2*V*W-X**2*U*V
F=Q**3+27*A**6*X**4*U*V**2*W**3-27*A**6*X**2*U**2*V*W**3+27*A**4*X**2*Q*U*V*W**2
P=sp.Poly(-F.subs(z,1+z),z)
assert P.degree()==6
c=P.all_coeffs() # c6,...,c0, in ZZ[A,X]
c0,c1,c2,c3,c4,c5,c6=sp.symbols('c0:7')
Hc=sp.Matrix([
    [c5,c3,c1,0,0,0],
    [c6,c4,c2,c0,0,0],
    [0,c5,c3,c1,0,0],
    [0,c6,c4,c2,c0,0],
    [0,0,c5,c3,c1,0],
    [0,0,c6,c4,c2,c0],
])
subst={c6:c[0],c5:c[1],c4:c[2],c3:c[3],c2:c[4],c1:c[5],c0:c[6]}
S=A+X+1
data=[]
expect=[(11,11,-1,78),(21,21,-1,253),(30,30,1,496),
        (36,36,1,703),(39,39,-1,820)]
for k,(da,dx,sgn,count) in enumerate(expect,1):
    print('forming Hurwitz minor',k,flush=True)
    delta=sp.Poly(sp.expand(Hc[:k,:k].det(method='domain-ge').subs(subst)),A,X,domain=sp.ZZ)
    print('  expanded',k,'terms',len(delta.terms()),flush=True)
    if k==1: denom=sp.Poly(S,A,X,domain=sp.ZZ)
    elif k==2: denom=sp.Poly(S**3,A,X,domain=sp.ZZ)
    elif k==3: denom=sp.Poly(S**6,A,X,domain=sp.ZZ)
    elif k==4: denom=sp.Poly(X**2*S**10,A,X,domain=sp.ZZ)
    else: denom=sp.Poly(X**6*S**15,A,X,domain=sp.ZZ)
    quotient,rem=sp.div(delta,denom)
    print('  divided',k,flush=True)
    assert rem.is_zero, (k,'factor division remainder')
    if k in (1,2,5): quotient=-quotient
    poly=quotient.as_expr()
    got=(sp.Poly(poly,A,X).degree(A),sp.Poly(poly,A,X).degree(X))
    assert got==(da,dx),(k,got,(da,dx))
    # Expand the positive-denominator-cleared substitution directly with
    # binomial coefficients instead of expanding a rational expression.
    transformed={}
    for (degA,degX),coef in quotient.terms():
        for degV in range(degX,dx+1):
            key=(degA,degV)
            transformed[key]=transformed.get(key,0)+int(coef)*5**degX*8**(dx-degX)*comb(dx-degX,degV-degX)
    coeffs=[x for x in transformed.values() if x]
    print('  substituted',k,'terms',len(coeffs),flush=True)
    signs={sp.sign(q) for q in coeffs}
    assert len(coeffs)==count,(k,len(coeffs),count)
    assert signs=={sgn},(k,signs,sgn)
    data.append({
        'hurwitz_minor':k,
        'factor_removed':{1:'S',2:'S^3',3:'S^6',4:'X^2*S^10',5:'X^6*S^15'}[k],
        'residual_degree_A_X':[da,dx],
        'residual_coefficient_count_after_X_substitution':len(coeffs),
        'residual_coefficient_sign':sgn,
        'transformed_v_degree':max(b for (a,b),coef in transformed.items() if coef),
    })

c0=P.nth(0)
assert sp.expand(c0-X**6*(A+X+1)**6)==0
assert len(data)==5
result={'sextic_degree':6,'constant_coefficient':'X^6*(A+X+1)^6',
        'hurwitz_residual_checks':data,
        'delta6_sign':'positive because c0>0 and delta5>0'}
path=Path(__file__).with_name('no79_hurwitz_certificate.json')
path.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
