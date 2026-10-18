"""Exact No. 52 algebra: ties, inner shape, endpoint Jacobians and positivity.

The outer valuation exhaustion and full six-body fractional-jet assertions
are proof obligations separate from these independently reconstructed rows.
"""
from itertools import combinations
from math import comb
from fractions import Fraction as Q
from pathlib import Path
import json
import sympy as s


def sparse(poly):return [[*m,str(c)] for m,c in poly.terms()]


def triangle_rows(a,b,d):
    u,v,p,q,r,t,k=s.symbols('u v p q r t k')
    variables=(u,v,p,q,r,t,k)
    rows=[b*r**3*(u-v)+d*t**3*u,
          a*r**3*(v-u)+d*k**3*v,
          b*r**3*(p-q)+d*t**3*p,
          a*r**3*(q-p)+d*k**3*q,
          r*r*(u-v)*(p-q)-1,t*t*u*p-1,k*k*v*q-1]
    return variables,rows


def main():
    delta,z,wi,wj,w5=s.symbols('delta z wi wj w5')
    expr=delta-delta**3*z*(wi-w5)-delta**3*z*(w5-wj)
    assert s.expand(expr-delta*(1-delta**2*z*(wi-wj)))==0
    # General six-body mass-centred pair-square identity.
    m=s.symbols('m1:7');Z=s.symbols('z1:7');W=s.symbols('w1:7')
    pair_sum=sum(m[i]*m[j]*(Z[i]-Z[j])*(W[i]-W[j]) for i,j in combinations(range(6),2))
    rhs=sum(m)*sum(a*b*c for a,b,c in zip(m,Z,W))-sum(a*b for a,b in zip(m,Z))*sum(a*b for a,b in zip(m,W))
    assert s.expand(pair_sum-rhs)==0
    A,C,t=s.symbols('A C t');u=-A/t**3;v=-t/A;cube=-C*t**3/(t**3+A)
    Fz4=cube*(1-u)+C;Fz5=A*cube*(u-1)+C*t**3*u
    Fw4=cube*(1-v)+C;Fw5=A*cube*(v-1)+C*t**3*v
    assert s.cancel(Fz4)==s.cancel(Fz5)==0
    assert s.cancel(Fw5-t**3*Fw4)==0
    assert s.cancel(Fw4-C*(A**2-t**4)/(A*(t**3+A)))==0
    squared=A*t**3/((t**3+A)*(A+t))
    shape=C**2*(A+t)**3*(t**3+A)-A**3*t**3
    assert s.cancel((cube**2-squared**3)/shape-t**6/((t**3+A)**3*(A+t)**3))==0
    positive_gap=s.Poly(s.expand((A+t)**3*(t**3+A)-A**3*t**3),A,t)
    assert all(c>0 for _,c in positive_gap.terms())
    print('Tie cancellation, six-body pair-square identity and inner shape verified',flush=True)

    x,y=s.symbols('x y',positive=True);a=x*x;b=y*y;d=x**3*y**3/(a+b)**2
    variables,rows=triangle_rows(a,b,d);u0,v0,p0,q0,r0,t0,k0=variables
    point={u0:1,v0:-a/b,p0:1,q0:-b/a,r0:-x*y/(a+b),t0:1,k0:1}
    J=s.Matrix(rows).jacobian(variables).subs(point)
    left=s.Matrix([-x**4/y**4,1,x**4/(3*y**4),-x**4/(3*y**4),
                   -x**5/(y*(a+b)),x**7/(y*(a+b)**2),x**5*y/(a+b)**2])
    assert all(s.cancel(q)==0 for q in J.T*left)
    gauge_kernel=s.Matrix([1,-a/b,-3,3*b/a,-x*y/(a+b),1,1])
    assert all(s.cancel(z)==0 for z in J*gauge_kernel)
    transverse=s.factor(J.extract([0,2,3,4,5,6],[1,2,3,4,5,6]).det(method='domain-ge'))
    expected_transverse=9*x**4*y**10*(x-y)**2*(x+y)**2/(x*x+y*y)**6
    assert s.cancel(transverse-expected_transverse)==0
    print('alpha=5/3 left kernel verified from the seven original triangle rows',flush=True)
    # Original crossed z sources at the first relative order 1/3.
    # The w-source at this order is common to both wings and projects to zero.
    Z5,p,T2,T3,X,Y=s.symbols('Z5 p T2 T3 X Y',nonzero=True)
    P=a+b;q=-b*p/a
    E2=T2/p;E3=T3/q
    projected=left[0]*(Z5-E2)+left[1]*(Z5-E3)
    Ucross=a*T2+b*T3
    assert s.cancel(projected-(a*Ucross-(a*a-b*b)*Z5*p)/(b*b*p))==0
    assert s.cancel(left[2]+left[3])==0
    reference=Z5*(a*p+b*q)-Ucross
    assert s.cancel(reference+(b*b*p/a)*projected)==0
    # z-centre: X*Z4+(P+d)*Z5=0, while Z6=0.
    Z4=-(P+d)*Z5/X
    D24sq=1/((Z5-Z4)*p);D26sq=1/(Z5*p)
    D34sq=1/((Z5-Z4)*q);D36sq=1/(Z5*q)
    assert s.cancel(D24sq/D26sq-X/(X+P+d))==0
    assert s.cancel(D34sq/D36sq-X/(X+P+d))==0
    assert s.cancel(D34sq/D24sq+a/b)==0
    assert s.cancel(D36sq/D26sq+a/b)==0
    for sign in (-1,1):
        coefficient=a+sign*s.I*x*y
        assert s.expand(coefficient*s.conjugate(coefficient)-a*P)==0
    # |Ucross|^3 = |(a^2-b^2)/a|*|cross coefficient|^2.
    assert s.cancel(((a*a-b*b)/a)*a*P-(a-b)*P**2)==0
    print('alpha=5/3 crossed-source solvability, reference cancellation, distance ratios and norm coefficient verified',flush=True)

    a,b,d=s.symbols('a b d',positive=True)
    variables,rows=triangle_rows(a,b,d);u0,v0,p0,q0,r0,t0,k0=variables
    # The sole equal-mass endpoint candidate has zero internal w-forces,
    # whereas the original two w rows require the same nonzero coefficient.
    equal_branch={b:a,d:a/4,u0:1,v0:-1,p0:1,q0:-1,
                  r0:-s.Rational(1,2),t0:1,k0:1}
    assert all(s.cancel(rows[i].subs(equal_branch))==0 for i in (2,3))
    print('alpha=7/4 equal-mass signed candidate has zero internal w-forces',flush=True)
    Js=s.Matrix(rows).jacobian(variables)
    branches=[(-1,1,a*a*b*b/((a+b)*(a*a+b*b)),
               -32*b**8*(a-b)**4*(a*a+a*b+b*b)**4/(a*(a+b)**7*(a*a+b*b)**4)),
              (-1,-1,a*a*b*b/((a+b)*(a*a-b*b)),
               32*b**8*(a*a-a*b+b*b)**4/(a*(a-b)**4*(a+b)**7)),
              (1,-1,-a*a*b*b/((a+b)*(a*a-b*b)),
               -32*b**8*(a*a-a*b+b*b)**4/(a*(a-b)**4*(a+b)**7))]
    branch_records=[]
    sources=s.Matrix([[1,0],[1,0],[0,1],[0,1],[0,0],[0,0],[0,0]])
    weighted=s.Matrix([[0,0,a,b,0,0,0]])
    for e23,e35,mass,expected in branches:
        point={u0:1,v0:-a/b,p0:1,q0:-a/b,r0:e23*b/(a+b),t0:1,k0:e35*b/a,d:mass}
        jb=Js.subs(point).applyfunc(s.cancel)
        determinant=s.factor(jb.det(method='domain-ge'))
        assert s.cancel(determinant-expected)==0,(e23,e35,determinant)
        response=jb.inv(method='DM')*sources
        projection=(weighted*response).applyfunc(s.cancel)
        assert projection==s.zeros(1,2)
        branch_records.append({'signs':[e23,e35],'d_mass':str(mass),
                               'Jacobian_determinant':str(determinant),
                               'weighted_w_source_responses':['0','0']})
        print(f'alpha=7/4 signs {(e23,e35)}: determinant and both source projections verified',flush=True)

    l,u,q,w=s.symbols('lambda u q w');L=1+l**3
    X=L*u**2/(1+2*u);v=u/(1+u);Y=l**3*q
    assert s.cancel(v**2-X/(X+L))==0
    H=l*l+l**3+L*(X+Y)+X*Y;C5=L+X+Y
    cleared=s.Poly(s.cancel((1+u)**2*(1+2*u)**3*(H**3-C5*(Y+X*v)**2)),q)
    assert cleared.degree()==3
    controls=[]
    for j in range(4):
        control=s.Poly(sum(s.Rational(comb(j,k),comb(3,k))*cleared.nth(k) for k in range(j+1)),l,u)
        e=control.degree(l)
        transformed=s.Poly(sum(c*w**aa*(1+2*w)**(e-aa)*u**bb for (aa,bb),c in control.terms()),w,u)
        assert all(c>=0 for (aa,bb),c in transformed.terms() if bb!=3)
        upoly=s.Poly(transformed.as_expr(),u)
        c2,c3,c4=[upoly.nth(k) for k in (2,3,4)]
        disc=s.Poly(4*c2*c4-c3*c3,w)
        assert all(c>=0 for _,c in disc.terms()) and not disc.is_zero
        constant=s.Poly(upoly.nth(0),w)
        assert all(c>=0 for _,c in constant.terms()) and not constant.is_zero
        controls.append({'Bernstein_index':j,'lambda_clearing_degree':e,
                         'transformed_control_sparse_w_u':sparse(transformed),
                         'four_c2_c4_minus_c3_squared_coefficients_w':sparse(disc),
                         'strict_constant_block_coefficients_w':sparse(constant)})
        print(f'Bernstein control {j}: nonnegative outside u^3 block; quadratic block and strict constant checked',flush=True)
    assert s.det(s.Matrix([[-s.Rational(1,2),-s.Rational(3,2)],[-s.Rational(3,2),-s.Rational(1,2)]]))==-2
    alpha=s.symbols('alpha')
    pairs=tuple(combinations(range(1,7),2));orders={ij:s.Integer(-1) for ij in pairs}
    for ij in ((4,5),(4,6),(5,6)):orders[ij]=s.Integer(2)
    for ij in ((2,3),(2,5),(3,5)):orders[ij]=4-2*alpha
    for ij in ((2,4),(2,6),(3,4),(3,6)):orders[ij]=5-3*alpha
    orders[(1,6)]=s.Integer(0)
    product_orders={ij:s.expand(orders[ij]+sum(o for kl,o in orders.items() if not(set(ij)&set(kl)))) for ij in pairs}
    for expr in product_orders.values():
        diff=s.Poly(expr-(13-8*alpha),alpha)
        assert diff.degree()<=1 or diff.is_zero
        assert diff.eval(1)>=0 and diff.eval(s.Rational(7,4))>=0
    assert {ij for ij,o in product_orders.items() if o==13-8*alpha}=={(2,6),(3,6),(5,6)}
    endpoint=Q(7,4);beta=4-3*endpoint;dcross=-(1+beta)/2
    source_orders={'crossed_z_relative':1+3*dcross-1,
                   'crossed_w_relative':beta+3*dcross-(-2),
                   'body1_z_relative':1+3*Q(1,2)-1,
                   'body1_w_relative':-2+3*Q(1,2)-(-2),
                   'coordinate_relative_z':endpoint-1,
                   'coordinate_relative_w':beta-(-2)}
    assert min(source_orders.values())==Q(3,8)
    potential_orders={'common_source_weighted_shape':1+beta+Q(3,8),
                      'relative_inertia':endpoint+beta,
                      'crossed_source_contraction':beta+1+3*dcross}
    assert min(potential_orders.values())==Q(1,8)
    result={'diagram':52,'tie_identity_corrected_index':'w5-wj, not w5-w3',
            'pair_square_identity':'sum_{i<j} mi*mj*(zi-zj)*(wi-wj)=M*I-(sum mi*zi)*(sum mi*wi)',
            'inner_shape':str(shape),'inner_force_ratio':'Fw5=t^3*Fw4',
            'inner_C_less_than_one_positive_gap_terms':sparse(positive_gap),
            'alpha_5_3_left_kernel_exact':True,
            'alpha_5_3_crossed_projection':'a*Ucross=(a^2-b^2)*Z5*p',
            'alpha_5_3_transverse_six_row_minor':str(transverse),
            'alpha_5_3_only_kernel_is_exact_scaling_family':True,
            'alpha_5_3_reference_potential_cancels':True,
            'alpha_5_3_crossed_ratio_correction':'D24^2/D26^2=X/(X+a+b+d), not D24^2/D25^2',
            'alpha_5_3_crossed_pair_squared_modulus':'a*(a+b)',
            'alpha_7_4_equal_mass_candidate_internal_w_forces_zero':True,
            'alpha_7_4_branches':branch_records,
            'source_scaling_determinant':'-2',
            'seven_product_orders':{''.join(map(str,ij)):str(e) for ij,e in product_orders.items()},
            'minimum_seven_product_order':'13-8*alpha on 1<=alpha<=7/4',
            'alpha_7_4_noncommon_source_relative_orders':{k:str(v) for k,v in source_orders.items()},
            'alpha_7_4_reference_potential_orders':{k:str(v) for k,v in potential_orders.items()},
            'alpha_7_4_uniform_IFT_error_order':'3/8',
            'alpha_7_4_UQ_order_lower_bound':'1/8',
            'correct_parameterization':'X=(1+lambda^3)*u^2/(1+2*u), v=u/(1+u)',
            'Bernstein_controls':controls,
            'strict_H5_cubed_gap_certified':True,
            'm1_derivative_bound':'C5-(Y^2*C5)^(1/3)-(1/3)*((Y+X*v)/(m1+C5))^(2/3)>3/4-1/3=5/12',
            'derivative_guard_facts':['C5>=1','0<Y/C5<1/8','0<(Y+X*v)/(m1+C5)<1'],
            'scope':'Exact endpoint algebra, positivity and order arithmetic; the manuscript proves the graph-to-order and implicit-function reductions used by these certificates.'}
    out=Path(__file__).with_name('no52_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f'Saved {out}')


if __name__=='__main__':main()
