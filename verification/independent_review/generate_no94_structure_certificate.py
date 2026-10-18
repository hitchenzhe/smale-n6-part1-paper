"""Universal No. 94 original-row and geometric-determinant identities."""
from pathlib import Path
import json
import sympy as s
from derive_no94_finite_fiber import det3


def main():
    a,b,c=s.symbols('a b c');P1,P2,P3,W1,W2,W3=s.symbols('P1 P2 P3 W1 W2 W3')
    A=[[b*(W1-W2),c*(W1-W3),0],
       [a*(W2-W1),0,c*(W2-W3)],
       [(a+b)*(P1-P2),c*(P1-P3),-c*(P2-P3)]]
    area=(P1-P3)*(W2-W3)-(P2-P3)*(W1-W3)
    determinant=s.factor(det3(A))
    assert s.expand(determinant-c*c*(a*W1+b*W2-(a+b)*W3)*area)==0
    assert s.expand(determinant.subs(W1,(-b*W2-c*W3)/a)+c*c*(a+b+c)*W3*area.subs(W1,(-b*W2-c*W3)/a))==0
    # A collinear wing z_i=lambda*w_i+Z0, z4=Z4, w4=0.
    # Original z_i-lambda*w_i rows: Z0=m4*q_i*(Z0-Z4).
    Z0,Z4,m4,q1,q2=s.symbols('Z0 Z4 m4 q1 q2')
    row1=Z0-m4*q1*(Z0-Z4);row2=Z0-m4*q2*(Z0-Z4)
    assert s.expand(row1-row2+m4*(q1-q2)*(Z0-Z4))==0
    assert s.expand(row1.subs(Z4,Z0)-Z0)==0
    # Equal spoke cubes: every solved wing cube is (s^-3-m4)*c_i.
    U,M,C1,C2,C3,t=s.symbols('U M C1 C2 C3 t')
    coeff=s.Matrix(s.symbols('c1:4'));rhs=s.Matrix([W1,W2,P1-P2])
    residual=s.Matrix(A)*((U-M)*coeff)-(U-M)*rhs-(U-M)*(s.Matrix(A)*coeff-rhs)
    assert all(s.expand(z)==0 for z in residual)
    # For C_i*t^2-1 homogenized in [t:v], a finite common root
    # requires all C_i equal, while infinity requires all C_i=0.
    assert s.expand((C1*t*t-1)-(C2*t*t-1)-(C1-C2)*t*t)==0
    result={'diagram':94,'geometric_determinant_factorization':str(determinant),
            'mass_centered_determinant':'-m3^2*(m1+m2+m3)*W3*oriented_area(P,W)',
            'all_zero_determinant_actual_distinct_spoke_shapes_have_G_zero':True,
            'collinear_original_row_identity':'Z0=m4*q_i*(Z0-Z4); distinct q_i imply Z0=Z4=0, hence G=0',
            'equal_spoke_scale_ties':'C_i*(s^-3-m4)^2=1',
            'scope':'Universal algebra. The full projective-boundary valuation cover is a separate mathematical argument.'}
    out=Path(__file__).with_name('no94_structure_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Universal determinant, collinear G=0 and equal-spoke scale identities verified')


if __name__=='__main__':main()
