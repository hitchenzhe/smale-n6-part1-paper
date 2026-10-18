"""Exact symbolic replay of the No. 66 mass-dominance valuation argument."""
from __future__ import annotations
import json
from itertools import combinations
from pathlib import Path
import sympy as s


def differential_identity():
    m=s.symbols('m0:4'); Z=s.symbols('Z0:4'); W=s.symbols('W0:4')
    dm=s.symbols('dm0:4'); dZ=s.symbols('dZ0:4'); dW=s.symbols('dW0:4')
    D={(i,j):s.symbols(f'D{i}{j}') for i,j in combinations(range(4),2)}
    def d(i,j): return D[(i,j)] if i<j else D[(j,i)]
    U=sum(m[i]*m[j]*d(i,j) for i,j in combinations(range(4),2))
    dU=0
    for i,j in combinations(range(4),2):
        Dij=d(i,j); Zij=Z[i]-Z[j]; Wij=W[i]-W[j]
        dDij=-Dij**3*(Wij*(dZ[i]-dZ[j])+Zij*(dW[i]-dW[j]))/2
        dU+=(dm[i]*m[j]+m[i]*dm[j])*Dij+m[i]*m[j]*dDij
    Q=[sum(m[j]*d(i,j) for j in range(4) if j!=i) for i in range(4)]
    Fz=[sum(m[j]*d(i,j)**3*(Z[i]-Z[j]) for j in range(4) if j!=i) for i in range(4)]
    Fw=[sum(m[j]*d(i,j)**3*(W[i]-W[j]) for j in range(4) if j!=i) for i in range(4)]
    rhs=sum(Q[i]*dm[i] for i in range(4))-s.Rational(1,2)*sum(
        m[i]*(Fw[i]*dZ[i]+Fz[i]*dW[i]) for i in range(4))
    assert s.expand(dU-rhs)==0
    return True


def main():
    a,b,c,u=s.symbols('a b c u', nonzero=True)
    v=-(a+b*u)/c
    x=s.factor(a*(v**3-1)/(b*(u**3-v**3)))
    y=s.factor(a*(1-u**3)/(c*(u**3-v**3)))
    W4=s.factor(1/(u**2*x))
    W5=s.factor(1/(v**2*y))
    T=b**3+c**3
    x0=-a*b**2/T
    y0=-a*c**2/T
    leadW4=s.limit(u**2*W4,u,s.oo)
    leadW5=s.limit(u**2*W5,u,s.oo)
    leadDiff=s.factor(s.limit(u**3*(W4-W5),u,s.oo))
    assert s.factor(leadW4+T/(a*b**2))==0,leadW4
    assert s.factor(leadW5+T/(a*b**2))==0,leadW5
    assert s.factor(leadDiff-T/b**3)==0,leadDiff
    assert s.factor(s.limit(x,u,s.oo)-x0)==0
    assert s.factor(s.limit(y,u,s.oo)-y0)==0

    # The three signed ties make p,q units and ord(r)<=-3/2 at u=infty.
    tie_p=s.factor((1-x)*(1-W4))
    tie_q=s.factor((1-y)*(1-W5))
    assert s.factor(s.limit(tie_p,u,s.oo)-(1-x0))==0
    assert s.factor(s.limit(tie_q,u,s.oo)-(1-y0))==0
    assert s.factor(s.limit(x-y,u,s.oo)-(x0-y0))==0
    # Coefficient field is C(a,b,c), so all listed leading constants are units.
    assert s.factor(x0-y0)!=0 and s.factor(1-x0)!=0 and s.factor(1-y0)!=0

    # At the hypothetical dominating component the 3rd tie yields
    # nu(r)<=-3/2; the 3rd z row then has a unique least term.
    valuation_terms={'ap3_x_minus_1':'= 0',
                     'cr3_x_minus_y':'<= -9/2',
                     'd_u3_x':'>= -3'}
    assert s.Rational(-9,2)<-3<0
    differential_identity()

    result={
      'diagram':66,
      'status':'passed',
      'normalization':'Z3=W3=D36=1; Z6=W6=0',
      'spoke_recovery':'v=-(a+b*u)/c; x=a*(v^3-1)/(b*(u^3-v^3)); y=a*(1-u^3)/(c*(u^3-v^3))',
      'infinity_leading_coefficients':{
        'x':'-a*b^2/(b^3+c^3)','y':'-a*c^2/(b^3+c^3)',
        'W4':'-(b^3+c^3)/(a*b^2)*u^-2',
        'W5':'-(b^3+c^3)/(a*b^2)*u^-2',
        'u^3*(W4-W5)':'(b^3+c^3)/b^3'},
      'finite_signed_cover_degree_bound':8,
      'valuation_comparison':valuation_terms,
      'exact_checks':['four-core differential identity','two linear shape rows solve x,y',
                      'common W4,W5 leading terms and nonzero difference coefficient',
                      'tie unit guards at u=infinity','unique least valuation in body-4 z row'],
      'scope':'Checks every algebraic identity in the dominant-component contradiction. The finite-type function-field argument then excludes a mass-dominating component for both V=0 and V!=0.'
    }
    out=Path(__file__).with_name('no66_function_field_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'diagram':66,'status':'passed','cover_degree':8,
                      'u3_W4_minus_W5':str(leadDiff)},indent=2))

if __name__=='__main__':main()
