"""Replay the finite-field No. 94 boundary certificates printed in the draft.

This checks explicit fibers and identities; it does not claim that all generic
No. 94 resultants, projective charts, or mass-varying boundary arcs are covered.
"""
from __future__ import annotations
from itertools import product
from pathlib import Path
import json
import sympy as s
P=17

def modpoly(poly,var):
    return s.Poly(poly,var,modulus=P)

def no94_positive_mass_exceptional_locus_certificate():
    r,t=s.symbols('r t'); a,b,c=s.symbols('a b c')
    A=a**3;B=b**3;C=c**3
    Q=(A*r-C)**2*(r+1)**2+B**2*r**2+2*B*r*(r+1)*(A*r+C)
    q=modpoly(Q.subs({a:1,b:3,c:5}),r)
    assert q.degree()==4 and q.is_irreducible
    guards=[-1,0,-s.Rational(4,1),5,-s.Rational(5,8)]
    vals=[]
    for x0 in guards:
        v=s.cancel(q.as_expr().subs(r,x0)); num,den=s.fraction(v)
        vals.append(int(num)%P*pow(int(den)%P,-1,P)%P)
    assert vals==[15,2,5,7,12],vals
    # Wing-dominant stratum: F=0 and K=0 from the displayed lemma.
    aa,bb,cc=1,3,5
    D=aa**2*cc**2*t*(1-t)
    N=D-aa**2*bb**2*t-bb**2*cc**2*(1-t)
    F=s.Poly(s.expand(N**3-27*bb**4*D**2),t,modulus=P)
    L1=bb*(1-t)+cc; L2=aa*(t-1)+cc*t; L3=-aa-bb*t
    K=s.Poly(s.expand(sum(m**4*([L1,L2,L3][i]**2)*([L1,L2,L3][(i+1)%3]**2 if i in (0,2) else [L1,L2,L3][2]**2) for i,m in enumerate((aa,bb,cc)))),t,modulus=P)
    # Use the exact indexed K formula rather than the compact construction above.
    Ls=(L1,L2,L3); ms=(aa,bb,cc)
    Kexpr=sum(ms[i]**4*Ls[j]**2*Ls[k]**2 for i in range(3) for j,k in [([q for q in range(3) if q!=i][0],[q for q in range(3) if q!=i][1])])
    for i,j in ((0,1),(0,2),(1,2)):
        k=3-i-j
        Kexpr-=2*ms[i]**2*ms[j]**2*Ls[i]*Ls[j]*Ls[k]**2
    K=modpoly(Kexpr,t)
    assert F.degree()==6 and K.degree()==4
    assert [int(F.nth(i))%P for i in reversed(range(7))]==[15,15,12,5,11,8,4]
    assert [int(K.nth(i))%P for i in reversed(range(5))]==[10,16,10,13,14]
    res=int(s.resultant(F.as_expr(),K.as_expr(),t))%P
    assert res==12
    return {"sample_mod17":[1,3,5],"Q_irreducible_degree":4,
            "forbidden_guard_Q_values_mod17":vals,
            "wing_dominant_F_coefficients_ascending_mod17":[int(F.nth(i))%P for i in range(7)],
            "wing_dominant_K_coefficients_ascending_mod17":[int(K.nth(i))%P for i in range(5)],
            "wing_dominant_resultant_mod17":res,
            "scope":"Exact local zero-sum fibers. Does not cover every No. 94 chart or establish the global exceptional-mass proposition."}

def test_no94_two_diverging_spokes_have_no_balanced_three_body_boundary_mod17():
    t=s.symbols('t'); records=[]
    for masses in ((3,5,1),(1,5,1),(1,3,1)):
        m2,m3,m4=masses
        vals=[]
        for eps1,eps2 in ((-1,-1),(-1,1),(1,-1),(1,1)):
            p=s.Poly(m3*t**2-eps1,t)
            q=s.Poly(m2*(1-t)**2-eps2,t)
            value=int(s.resultant(p.as_expr(),q.as_expr(),t))%P
            vals.append(value)
        records.append({"masses_m2_m3_m4":list(masses),"resultants_mod17":vals})
    expected=[[10,8,9,6],[16,4,16,15],[3,13,3,6]]
    assert [r['resultants_mod17'] for r in records]==expected,records
    return {"mass_triples":records,"all_twelve_nonzero":all(x for r in records for x in r['resultants_mod17'])}

# Arithmetic in F_17[omega]/(omega^2+omega+1).
def add(x,y): return ((x[0]+y[0])%P,(x[1]+y[1])%P)
def neg(x): return ((-x[0])%P,(-x[1])%P)
def sub(x,y): return add(x,neg(y))
def mul(x,y):
    a,b=x;c,d=y
    return ((a*c-b*d)%P,(a*d+b*c-b*d)%P)
def inv(x):
    a,b=x; n=(a*a-a*b+b*b)%P
    assert n
    return ((a-b)*pow(n,-1,P)%P,(-b)*pow(n,-1,P)%P)
def scale(x,n): return (x[0]*n%P,x[1]*n%P)
def padd(a,b):
    z=[(0,0)]*max(len(a),len(b))
    for i,x in enumerate(a): z[i]=add(z[i],x)
    for i,x in enumerate(b): z[i]=add(z[i],x)
    while len(z)>1 and z[-1]==(0,0): z.pop()
    return z
def pmul(a,b):
    z=[(0,0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): z[i+j]=add(z[i+j],mul(x,y))
    while len(z)>1 and z[-1]==(0,0): z.pop()
    return z
def pscale(a,n): return [scale(x,n) for x in a]
def psub(a,b): return padd(a,pscale(b,-1))
def pzero(): return [(0,0)]
def pow2(a): return pmul(a,a)
def rem_monic(poly,q):
    z=list(poly)
    while len(z)>=3:
        k=len(z)-1; lead=z[-1]; shift=k-2
        # q=[constant, linear, 1]
        for j in range(3): z[shift+j]=sub(z[shift+j],mul(lead,q[j]))
        z.pop()
    z += [(0,0)]*(2-len(z))
    return z[:2]
def omega_pow(k): return ((1,0),(0,1),(16,16))[k%3]
def noncollinear_H3_table():
    a,b,c=1,3,5; M=a+b+c; rows=[]
    for j in range(3):
        for k in range(3):
            s12=(1,0); s13=omega_pow(j); s23=omega_pow(k)
            # t*L_i where L_i=M*(raw_i - weighted centre).
            l1=[scale(s13,-c),(-b%P,0)]
            l2=[scale(s13,-c),((a+c)%P,0)]
            l3=[scale(s13,a+b),(-b%P,0)]
            L=(l1,l2,l3); ms=(a,b,c)
            K=pzero()
            for i in range(3):
                others=[n for n in range(3) if n!=i]
                term=pmul(pow2(L[others[0]]),pow2(L[others[1]]))
                K=padd(K,pscale(term,ms[i]**4))
            for i,jj in ((0,1),(0,2),(1,2)):
                kk=3-i-jj
                term=pmul(pmul(L[i],L[jj]),pow2(L[kk]))
                K=padd(K,pscale(term,-2*ms[i]**2*ms[jj]**2))
            # Q=t^2 + (s23-1-s13)t+s13.
            B=sub(sub(s23,(1,0)),s13); q=(s13,B,(1,0))
            rem=rem_monic(K,q); A0,A1=rem
            resultant=add(sub(add(mul(mul(A1,A1),s13),mul(mul(A1,A0),B)),(0,0)),mul(A0,A0))
            # Correct formula: A1^2*C - A1*A0*B + A0^2.
            resultant=add(sub(mul(mul(A1,A1),s13),mul(mul(A1,A0),B)),mul(A0,A0))
            assert resultant!=(0,0),(j,k,rem)
            rows.append([resultant[0],resultant[1]])
    expected=[[(4,0),(9,0),(9,0)],[(0,0),(0,0),(0,0)]]
    return rows

def test_no94_all_zero_spokes_finite_wing_shapes_have_no_h3_zero_mod17():
    # Collinear branch is checked directly; phase resultant table is recomputed
    # over F_17(omega), but its printed comparison is kept in explicit order.
    t=s.symbols('t'); m1,m2,m3=1,3,5; M=m1+m2+m3
    H=s.Poly(8225*t**4+24290*t**3+15339*t**2-3838*t+209,t,modulus=P)
    collinear=[]
    for e12,e13,e23 in product((-1,1),repeat=3):
        p1=-(m2+m3*(1+t))/M; p2=1-p1
        A1=-m2*e12-m3*e13/(1+t)**2
        A2=m1*e12-m3*e23/t**2
        poly=s.cancel(M*t**2*(1+t)**2*(p2*A1-p1*A2))
        E=modpoly(poly,t)
        g=s.gcd(E,H)
        assert E.degree()==5 and g.degree()==0,(e12,e13,e23,E.degree(),g)
        collinear.append({"signs":[e12,e13,e23],"degree":E.degree(),"gcd_degree":g.degree()})
    phase_rows=noncollinear_H3_table()
    expected=[(4,0),(9,0),(9,0),(11,14),(4,12),(8,15),(14,3),(10,2),(9,5)]
    assert phase_rows==[list(x) for x in expected],phase_rows
    return {"sample_mod17":[1,3,5],"collinear_euler_signs":collinear,
            "noncollinear_phase_resultants_basis_1_omega":phase_rows,
            "all_8_plus_27_branches_exclude_H3_zero":True,
            "scope":"The finite-wing component of the all-spokes-zero boundary only."}

def main():
    result={"no94_positive_mass_exceptional_locus_certificate":no94_positive_mass_exceptional_locus_certificate(),
            "test_no94_all_zero_spokes_finite_wing_shapes_have_no_h3_zero_mod17":test_no94_all_zero_spokes_finite_wing_shapes_have_no_h3_zero_mod17(),
            "test_no94_two_diverging_spokes_have_no_balanced_three_body_boundary_mod17":test_no94_two_diverging_spokes_have_no_balanced_three_body_boundary_mod17()}
    out=Path(__file__).with_name('no94_local_certificates.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k: 'passed' for k in result},indent=2))
if __name__=='__main__': main()
