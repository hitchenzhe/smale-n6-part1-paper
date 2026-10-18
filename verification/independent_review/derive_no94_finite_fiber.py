"""Derive No. 94 wing cubes from original w rows and one z-row difference.

First reconstruct the finite-field slice x=6 from the original rows, using
the spoke roots (1,x,y) and their common cube scale S. No precomputed quotient
or norm polynomial is used.
"""
from itertools import combinations
from pathlib import Path
import json
import sympy as s
from sympy.polys.fields import field


def det3(A):
    return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
            -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
            +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))


def solve3(A,b,det):
    out=[]
    for j in range(3):
        B=[r[:] for r in A]
        for i in range(3):B[i][j]=b[i]
        out.append(det3(B)/det)
    return out


def derive(K,x,y,masses):
    a,b,c,d=[K(m) for m in masses]
    roots=[K(1),x,y];q=[r**3 for r in roots]
    W=[-(q[2]-q[1])/a,(q[2]-q[0])/b,-(q[1]-q[0])/c]
    P=[1/(r*r*w) for r,w in zip(roots,W)]
    # Coefficients of (p12,p13,p23) from two wing w rows and z1-z2.
    A=[[b*(W[0]-W[1]),c*(W[0]-W[2]),K(0)],
       [a*(W[1]-W[0]),K(0),c*(W[1]-W[2])],
       [(a+b)*(P[0]-P[1]),c*(P[0]-P[2]),-c*(P[1]-P[2])]]
    b0=[W[0],W[1],P[0]-P[1]]
    b1=[-d*q[0]*W[0],-d*q[1]*W[1],-d*(q[0]*P[0]-q[1]*P[1])]
    Delta=det3(A)
    c0=solve3(A,b0,Delta);c1=solve3(A,b1,Delta)
    ties=[]
    for k,(i,j) in enumerate(combinations(range(3),2)):
        factor=((P[i]-P[j])*(W[i]-W[j]))**3
        # In ascending S order: factor*(c0+c1*S)^2-S^2.
        ties.append((factor*c0[k]**2,2*factor*c0[k]*c1[k],factor*c1[k]**2-1))
    guards=[*roots,*W,*[W[i]-W[j] for i,j in combinations(range(3),2)],
            *[q[i]-q[j] for i,j in combinations(range(3),2)]]
    return Delta,ties,guards


def quadratic_resultant(f,g):
    c,b,a=f;F,e,d=g
    return (a*F-c*d)**2-(a*e-b*d)*(b*F-c*e)


def main():
    K,y=field('y',s.GF(17));x=K(6)
    delta,ties,guards=derive(K,x,y,(1,2,4,1))
    raw=[quadratic_resultant(ties[0],ties[j]) for j in (1,2)]
    yy=s.Symbol('y')
    polys=[s.Poly(r.numer.as_expr(),yy,modulus=17) for r in raw]
    print('numerator degrees',[p.degree() for p in polys],flush=True)
    g=s.gcd(*polys)
    print('raw gcd degree',g.degree(),flush=True)
    guard=s.Poly(1,yy,modulus=17)
    for h in guards+[delta]:
        guard*=s.Poly(h.numer.as_expr(),yy,modulus=17)*s.Poly(h.denom.as_expr(),yy,modulus=17)
    removed=[]
    while g.degree()>0:
        h=s.gcd(g,guard)
        if h.degree()==0:break
        removed.append(h.degree());g=g.exquo(h)
    print('remaining gcd degree after guards',g.degree(),flush=True)
    result={'diagram':94,'method':'original-row derivation of three affine wing cubes and squared signed ties',
            'mass_sample_mod17':[1,2,4,1],'fixed_x':6,
            'raw_resultant_numerator_degrees':[p.degree() for p in polys],
            'raw_common_factor_degree':s.gcd(*polys).degree(),
            'remaining_guarded_gcd_degree':g.degree(),
            'removed_guard_gcd_degrees':removed,
            'resultant_numerators_ascending_mod17':[[int(p.nth(i))%17 for i in range(p.degree()+1)] for p in polys],
            'geometric_determinant_numerator':str(delta.numer.as_expr()),
            'geometric_determinant_denominator':str(delta.denom.as_expr())}
    out=Path(__file__).with_name('no94_finite_fiber_slice_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f'Saved {out}')


if __name__=='__main__':main()
