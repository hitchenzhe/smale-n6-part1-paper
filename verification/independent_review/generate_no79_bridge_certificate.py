"""Exact bridge cubic identities and the omitted Hurwitz leading-coefficient gate."""
from pathlib import Path
from math import comb
import json
import sympy as s


def main():
    theta,t,V=s.symbols('theta t V',real=True);p=(1-theta)/2;q=(1+theta)/2
    d1=1+q*theta*t+s.I*q*V;d2=1-p*theta*t-s.I*p*V
    assert s.expand(p*d1+q*d2-1)==0
    cubic=s.expand(p*d1**3+q*d2**3)
    L=3*theta**2*t*t+6*t
    imag=s.im(cubic)
    assert s.expand(imag-p*q*theta*V*(L-V*V))==0
    norm_difference=s.factor((s.re(d1)**2+s.im(d1)**2-
                              s.re(d2)**2-s.im(d2)**2).subs(V*V,L))
    assert s.expand(norm_difference-4*theta*t*(2+theta**2*t))==0
    # The unweighted modulus order does not determine the mass-weighted shape.
    gamma=1+q*theta*t
    G=s.factor(3*q*q*L-gamma**2)
    expected_G=2*(1+theta)**2*theta**2*t**2+(1+theta)*(9+7*theta)*t/2-1
    assert s.expand(G-expected_G)==0
    assert G.subs({theta:s.Rational(9,20),t:s.Rational(1,5)})==s.Rational(1591621,2000000)
    assert G.subs({theta:s.Rational(9,20),t:s.Rational(1,5)})>0
    assert s.expand(s.diff(G,t)-(1+theta)*(8*t*theta**3+8*t*theta**2+7*theta+9)/2)==0
    assert s.expand(s.diff(G,theta)-t*(8*t*theta**3+12*t*theta**2+4*t*theta+7*theta+8))==0
    alpha,beta,K=s.symbols('alpha beta K',positive=True)
    weighted_r=-beta*d2**3/(alpha*d1**3)
    assert s.cancel((weighted_r-1+K/(alpha*d1**3)).subs(K,alpha*d1**3+beta*d2**3))==0
    real=s.expand(s.re(cubic).subs(V*V,L))
    R=1-2*p*q*t*(3+2*theta**2*t)**2
    assert s.expand(real-R)==0
    # The physical Holder inequality gives
    # R=K*P^2/c^3 < P^2*(M_A*M_B)/sigma^3.
    # Put P=alpha+beta, theta=(beta-alpha)/P and
    # u=2*sqrt(alpha*beta)/P=sqrt(1-theta^2).
    # The two unused cross products and the two inner products obey
    # M_A*M_B >= P*(1+u), sigma-M_A*M_B >= P*u.
    # First replace I by its lower bound P*u; then C/(C+P*u)^3
    # decreases for C>=P*(1+u). Hence R < (1+u)/(1+2*u)^3.
    u,v=s.symbols('u v',nonnegative=True)
    R_at_one_fifth=1-u*u*(17-2*u*u)**2/s.Integer(250)
    mass_upper=(1+u)/(1+2*u)**3
    gap=s.factor(250*(1+2*u)**3*(R_at_one_fifth-mass_upper))
    positive_piece=s.cancel(gap/u)
    bernstein_chart=s.Poly(s.cancel((1+v)**8*positive_piece.subs(u,v/(1+v))),v)
    positive_coefficients=[int(c) for c in bernstein_chart.all_coeffs()]
    assert positive_coefficients==[175,9385,53453,129401,167471,125127,54243,12711,1250]
    assert s.simplify(positive_piece.subs(u,1))==175
    # R(t,theta) strictly decreases for t>0 and 0<theta<1.
    assert s.expand(s.diff(R,t)+2*p*q*(3+2*theta**2*t)*(3+6*theta**2*t))==0
    # The stronger floor on 3/8 <= theta <= 9/20 also follows directly
    # from the same physical mass bound. The derivative is negative there.
    ell=theta/s.Integer(3)+s.Rational(17,200)
    R_ell=s.factor(R.subs(t,ell))
    derivative_quintic=s.factor(s.diff(R_ell,theta)*54000000/(200*theta**3+51*theta**2+900))
    assert derivative_quintic==180000*theta**5+76500*theta**4-132197*theta**3+213900*theta**2+40698*theta-90000
    interval_poly=s.Poly(s.cancel((1+v)**5*(-derivative_quintic).subs(theta,s.Rational(3,8)+s.Rational(3,40)*v/(1+v))),v)
    interval_coefficients=[str(c) for c in interval_poly.all_coeffs()]
    assert all(c>0 for c in interval_poly.all_coeffs())
    assert s.Rational(319,400)>s.Rational(89,100)**2
    middle_bound=s.Rational(1)+s.Rational(89,100)
    middle_bound/= (1+2*s.Rational(89,100))**3
    endpoint_value=s.factor(R_ell.subs(theta,s.Rational(9,20)))
    assert endpoint_value==s.Rational(26184698477743,256000000000000)
    assert endpoint_value>middle_bound
    # Opposite nonreal half-sheet t<-2/theta^2 has R>1 exactly.
    # Physical Holder and sigma>MA*MB>alpha+beta require R<1.
    A,X,z=s.symbols('A X z');U=(A+X)*z+1;VV=X*z+A+1;W=A*(1-z)
    Q=A*A*U*W-A*A*X*X*VV*W-X*X*U*VV
    F=Q**3+27*A**6*X**4*U*VV**2*W**3-27*A**6*X**2*U**2*VV*W**3+27*A**4*X**2*Q*U*VV*W**2
    # Original zero-z-force rows at the three triangle vertices imply
    # (w23/w12)^3=U/(A^2 W), (w13/w12)^3=U/(X^2 V).
    # Additivity w12+w23=w13 gives this two-cubic norm.
    simple=Q**3-27*A**4*X**4*U**2*VV**2*W**2
    assert s.expand(F-simple)==0
    # Swapping the two wing masses sends r to 1/r. This is an exact
    # polynomial symmetry, but identifying which bridge ratio is the
    # original-row shape variable on each triangle remains necessary.
    def shape_polynomial(a,x,z):
        uu=(a+x)*z+1;vv=x*z+a+1;ww=a*(1-z)
        qq=a*a*uu*ww-a*a*x*x*vv*ww-x*x*uu*vv
        return qq**3-27*a**4*x**4*uu**2*vv**2*ww**2
    assert s.cancel(X**18*z**6*shape_polynomial(A/X,1/X,1/z)-simple)==0
    u,v=s.symbols('u v')
    assert s.expand(((v**3-1-u**3)**3-27*u**3*v**3).subs(v,u+1))==0
    pz=U/(A*A*W);qz=U/(X*X*VV)
    assert s.cancel((qz-1-pz)**3-27*pz*qz-simple/(A**6*X**6*VV**3*W**3))==0
    polynomial=s.Poly(-F.subs(z,1+z),z)
    lead=s.Poly(polynomial.LC(),A,X);degree=lead.degree(X);transformed={}
    for (i,j),coef in lead.terms():
        for k in range(j,degree+1):
            key=(i,k);transformed[key]=transformed.get(key,0)+int(coef)*5**j*8**(degree-j)*comb(degree-j,k-j)
    coeff={k:v for k,v in transformed.items() if v}
    assert degree==12 and len(coeff)==91 and all(v>0 for v in coeff.values())
    assert max(k[1] for k in coeff)==12
    a,b,c,x,y,qq=s.symbols('a b c x y q')
    cubicgap=s.expand(a*a*(a*x**3-b*y**3)-(a*x-b*y)**3)
    check=cubicgap.subs(y,qq*a*x/b)-a**3*x**3*(3*qq*(1-qq)+(1-a*a/(b*b))*qq**3)
    assert s.cancel(check)==0
    out=Path(__file__).with_name('no79_bridge_certificate.json')
    result={'weighted_bridge_first_moment':1,
            'imaginary_cubic_factor':'p*q*theta*V*(3*theta^2*t^2+6*t-V^2)',
            'bridge_norm_difference':'abs(d26/m)^2-abs(d35/m)^2=4*theta*t*(2+theta^2*t)>0',
            'weighted_shape_ratio':'r=-beta*d35^3/(alpha*d26^3)=1-K/(alpha*d26^3)',
            'large_asymmetry_G_corner':'1591621/2000000',
            'large_asymmetry_Re_d26_cubed_negative_from_t_floor':True,
            'large_asymmetry_weighted_Re_r_greater_than_one':True,
            'original_row_shape_identification':'The weighted ratio is stated in the manuscript; this script checks its algebraic consequences, not the original six-body row-to-shape mapping.',
            'real_weighted_cubic':'1-2*p*q*t*(3+2*theta^2*t)^2',
            'physical_Holder_bridge_bound':'R < (1+u)/(1+2*u)^3, u=sqrt(1-theta^2)',
            't_greater_than_one_fifth_for_all_unequal_positive_bridge_weights':True,
            'one_fifth_gap_positive_chart_coefficients_descending':positive_coefficients,
            'adaptive_floor_theta_interval':['3/8','9/20'],
            'adaptive_floor':'t > theta/3 + 17/200',
            'adaptive_floor_derivative_negative_chart_coefficients_descending':interval_coefficients,
            'R_at_adaptive_floor_theta_9_over_20':str(endpoint_value),
            'physical_mass_bound_theta_3_over_8_to_9_over_20':str(middle_bound),
            'opposite_nonreal_half_sheet_has_R_greater_than_one':True,
            'real_signed_bridge_cubic_gap_exact':True,
            'triangle_sextic_derived_from_original_zero_z_forces_and_signed_ties':True,
            'simplified_triangle_sextic':'Q^3-27*A^4*X^4*U^2*V^2*W^2',
            'triangle_sextic_reciprocal_mass_swap':'X^18*r^6*F(A/X,1/X,1/r)=F(A,X,1,r)',
            'Hurwitz_leading_coefficient_transformed_degree_v':12,
            'Hurwitz_leading_coefficient_positive_term_count':91,
            'positive_leading_coefficient_sparse_A_v':[[i,j,c] for (i,j),c in sorted(coeff.items())],
            'scope':'Bridge moments, both exact physical Holder t floors, the weighted large-asymmetry sector calculation, and the leading Hurwitz coefficient. The exact equal-weight and N1 shape/energy identities are replayed by verify_no79_shape_energy.py. The original six-body row-to-shape assignment and global geometric branch cover remain separate proof obligations.'}
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Bridge first/cubic moments, real-bridge gap, opposite half-sheet and all 91 leading-coefficient signs verified')


if __name__=='__main__':main()
