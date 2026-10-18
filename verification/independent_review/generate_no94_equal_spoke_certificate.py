"""Check every equal-spoke cube phase over the complete H3 quartic field.

With W=(-1/a,(r+1)/b,-r/c), H3 fixes the irreducible quartic r.
For equal spoke cubes D_i4=s*e_i, e_i^3=1, the original wing rows
give the three normalized wing cubes as (s^-3-m4)*c_i.
The three signed ties are C_i*t^2=1. Unequal C_i exclude every
finite t and their nonzero differences exclude a common projective infinity.
Only the phase e_1=1 is fixed; the other two have nine choices.
"""
from itertools import combinations
from pathlib import Path
import json
from flint import fmpz_mod_poly_ctx,fq_default_ctx
from derive_no94_finite_fiber import det3,solve3


def main():
    a,b,c,d=1,3,5,1
    import sympy as sp
    r=sp.Symbol('r');aa,bb,cc=a**3,b**3,c**3
    Q=sp.Poly((aa*r-cc)**2*(r+1)**2+bb**2*r**2+2*bb*r*(r+1)*(aa*r+cc),r)
    P=fmpz_mod_poly_ctx(17);q=P([int(Q.nth(k))%17 for k in range(5)])
    assert q.degree()==4 and q.is_irreducible()
    F=fq_default_ctx(modulus=q,var='r',fq_type='FQ_NMOD');r=F.gen();one=F.one();zero=F.zero()
    omega=None
    for offset in range(17):
        candidate=(r+offset)**((17**4-1)//3)
        if candidate!=one and candidate**3==one:
            omega=candidate;break
    assert omega is not None
    phases=(one,omega,omega**2);W=(-one/a,(r+1)/b,-r/c)
    records=[]
    for j,ej in enumerate(phases):
        for k,ek in enumerate(phases):
            roots=(one,ej,ek);Z=[1/(ei*ei*wi) for ei,wi in zip(roots,W)]
            A=[[b*(W[0]-W[1]),c*(W[0]-W[2]),zero],
               [a*(W[1]-W[0]),zero,c*(W[1]-W[2])],
               [(a+b)*(Z[0]-Z[1]),c*(Z[0]-Z[2]),-c*(Z[1]-Z[2])]]
            determinant=det3(A)
            assert determinant!=zero,('singular phase',j,k)
            coef=solve3(A,[W[0],W[1],Z[0]-Z[1]],determinant)
            C=[((Z[i]-Z[jj])*(W[i]-W[jj]))**3*coef[n]**2
               for n,(i,jj) in enumerate(combinations(range(3),2))]
            assert C[0]!=C[1] or C[0]!=C[2]
            assert any(z!=zero for z in C)
            records.append({'cube_phase_exponents':[j,k],
                            'row_determinant':str(determinant),
                            'quadratic_tie_coefficients_C':[str(z) for z in C],
                            'C1_minus_C2':str(C[0]-C[1]),'C1_minus_C3':str(C[0]-C[2]),
                            'no_common_finite_or_infinite_scale':True})
    result={'diagram':94,'mass_sample_mod17':[a,b,c,d],
            'H3_quartic_coefficients_ascending_mod17':[int(q[i]) for i in range(5)],
            'H3_quartic_irreducible':True,'omega_as_polynomial_in_r':str(omega),
            'four_Frobenius_roots_covered':True,'nine_phase_records':records,
            'all_equal_spoke_projective_scale_fibers_empty':True,
            'scope':'Equal-spoke boundary chart at the H3 sample; other boundary charts require their own arguments.'}
    out=Path(__file__).with_name('no94_equal_spoke_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('All nine equal-spoke phases, all four H3 ratio conjugates, finite and infinite scale excluded')
    print(f'Saved {out}')


if __name__=='__main__':main()
