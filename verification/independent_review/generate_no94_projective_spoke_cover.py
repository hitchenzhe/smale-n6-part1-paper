"""Exact projective spoke-ratio partition on the No. 94 H3=0 fiber.

This certifies the list of relative spoke-ratio charts. It does not
certify the force-equation exclusions within each chart.
"""
from pathlib import Path
import json
import sympy as s


def main():
    r,L,a,b,c=s.symbols('r L a b c')
    wing_roots=(a,-b/(r+1),c/r)
    H=sum(wing_roots)**2-4*sum(wing_roots[i]*wing_roots[j]
                                     for i in range(3) for j in range(i+1,3))
    Q=(a*r-c)**2*(r+1)**2+b*b*r*r+2*b*r*(r+1)*(a*r+c)
    assert s.cancel(r*r*(r+1)**2*H-Q)==0
    assert s.expand(Q.subs(r,0)-c*c)==0
    assert s.expand(Q.subs(r,-1)-b*b)==0
    assert s.expand(s.Poly(Q,r).LC()-a*a)==0
    q=(L,L+r,L+r+1)
    assert [s.expand(q[1]-q[0]),s.expand(q[2]-q[1]),s.expand(q[2]-q[0])]==[r,1,r+1]
    zero_locations=(s.Integer(0),-r,-r-1)
    assert [s.expand(zero_locations[1]-zero_locations[0]),
            s.expand(zero_locations[2]-zero_locations[1]),
            s.expand(zero_locations[2]-zero_locations[0])]==[-r,-1,-r-1]
    q17=s.Poly(Q.subs({a:1,b:27,c:125}),r,modulus=17)
    assert q17.degree()==4 and q17.is_irreducible
    assert [int(q17.eval(0))%17,int(q17.eval(-1))%17,int(q17.LC())%17]==[2,15,1]
    data={
        'diagram':94,
        'exact_H3_quartic_from_original_weighted_roots':True,
        'Q_at_r_zero':'c^2', 'Q_at_r_minus_one':'b^2', 'Q_at_infinity_leading':'a^2',
        'mass_sample_mod17':[1,3,5,1],
        'sample_Q_zero_minus_one_infinity_mod17':[2,15,1],
        'projective_spoke_ratio':'[q1:q2:q3]=[L:L+r:L+r+1]',
        'chart_partition':[
            'finite L outside 0,-r,-r-1: three nonzero relative spokes',
            'finite L in 0,-r,-r-1: exactly one zero relative spoke',
            'L=infinity: all three relative spoke cubes equal'],
        'scope':'Complete projective ratio partition on H3=0; original force-equation exclusion within each chart is a separate obligation.'
    }
    out=Path(__file__).with_name('no94_projective_spoke_cover.json')
    out.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print('No. 94 H3 quartic and three-chart projective spoke partition verified')


if __name__=='__main__':
    main()
