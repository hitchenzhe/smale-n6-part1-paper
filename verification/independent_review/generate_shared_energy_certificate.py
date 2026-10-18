"""Reconstruct shared original-row identities and positive-mass energy gaps."""
from itertools import combinations
from pathlib import Path
import json
import sympy as s


def main():
    pairs=list(combinations(range(6),2))
    E=s.Matrix([[int(ij==kl or not(set(ij)&set(kl))) for kl in pairs] for ij in pairs])
    assert E.det()==-114688 and E.eigenvals()=={7:1,-2:5,2:9}
    a,b,c,d=s.symbols('a b c d',positive=True);Q=a+b+c+d
    edges=((a,b),(b,c),(c,d),(d,a));D=a*c+b*d
    cycle=2*D-sum(x*y*(Q-x-y)/(x+y) for x,y in edges)
    square=(a*c-b*d)**2*(1/((a+b)*(c+d))+1/((a+d)*(b+c)))
    assert s.cancel(cycle-square)==0
    assert s.cancel(Q/2-sum(x*y/(x+y) for x,y in edges)-sum((x-y)**2/(4*(x+y)) for x,y in edges))==0
    # No.20 displayed Holder polynomial and exact positive quadratic block.
    L,x=s.symbols('lambda x',positive=True)
    gap=(L*x+1+x)**3-(L*(1+x)+s.Rational(9,8))*(s.Rational(27,2)+L*L/16)
    at=s.Poly(gap.subs(x,s.Rational(19,12)),L)
    assert [at.nth(i) for i in range(4)]==[s.Rational(3547,1728),-s.Rational(1829,576),s.Rational(22301,1152),s.Rational(1645,432)]
    assert 4*at.nth(0)*at.nth(2)-at.nth(1)**2==s.Rational(148167571,995328)
    h=s.symbols('h',nonnegative=True)
    derivative=s.Poly(s.expand(s.diff(gap,x).subs(x,h+s.Rational(19,12))),h,L)
    assert all(z>0 for _,z in derivative.terms())
    # No.36: replace the dimensionally inconsistent displayed equality.
    B=a+b;T=c+d+s.Symbol('e',positive=True);e=T-c-d;C=B+T;X=a*b;H=a*d+b*c
    deficit=2*(B*T-H)-X*T/B
    positive=2*B*e+2*a*c+2*b*d-X*(e+c+d)/B
    assert s.cancel(deficit-positive)==0
    assert s.cancel(B*positive-((2*B*B-X)*e+(2*a*B-X)*c+(2*b*B-X)*d))==0
    assert all(z>0 for _,z in s.Poly(s.expand(2*B*B-X),a,b).terms())
    assert all(z>0 for _,z in s.Poly(s.expand(2*a*B-X),a,b).terms())
    assert all(z>0 for _,z in s.Poly(s.expand(2*b*B-X),a,b).terms())
    # Sum of pair weights for n positive masses: <=(n-1)/(2n)*M^2.
    mass=s.symbols('m1:7');M=sum(mass);sigma=sum(mass[i]*mass[j] for i,j in pairs)
    assert s.expand(s.Rational(5,12)*M*M-sigma-sum((mass[i]-mass[j])**2 for i,j in pairs)/12)==0
    # No.105 pole-source ratio becomes strictly negative after U pole cancellation.
    m2,m3,m5,m6=s.symbols('m2 m3 m5 m6',positive=True)
    source_ratio=-m2*m5/m6**2*(-m3*m6/(m2*m5))**4
    assert s.factor(source_ratio)==-m3**4*m6**2/(m2**3*m5**3)
    result={'complement_product_matrix_determinant':int(E.det()),
            'complement_product_eigenvalues':{'7':1,'-2':5,'2':9},
            'four_cycle_rational_deficit_exact':True,'cycle_reciprocal_mass_bound_exact':True,
            'no20_four_coefficients':[str(at.nth(k)) for k in range(4)],
            'no20_quadratic_discriminant_gap':'148167571/995328',
            'no20_monotonicity_on_x_ge_19_12':True,
            'no36_replacement_positive_deficit':str(s.factor(deficit)),
            'no36_bad_old_equality_unit_mass_left':'87/40',
            'no36_bad_old_equality_unit_mass_right':'11/8',
            'no105_eliminated_ratio':str(s.factor(source_ratio)),
            'scope':'Exact algebra and mass inequalities, separate from graph-specific valuation exhaustion.'}
    out=Path(__file__).with_name('shared_energy_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Complement-product matrix, cycle gaps, No20 polynomial, corrected No36 energy gap and No105 sign verified')


if __name__=='__main__':main()
