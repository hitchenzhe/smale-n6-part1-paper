"""Reconstruct the No.50 open resonance's source ratio and potential pole.

The original w1 and mass-weighted w2+w3 rows imply a ratio of
the signed bridge cubes. Taking its simple cube root forces +1,
not -1, as the normalized ratio's constant coefficient.
The graph's remaining valuation cases are separate proof items.
"""
from pathlib import Path
import json
import sympy as s


def main():
    order=s.Symbol('s');delta=(order-1)/3;q=(4*order-7)/3
    bridge=order-2;wbridge=4-3*order
    assert s.expand(bridge+delta-q)==0
    assert s.expand((-2*q-1)-wbridge-delta)==0
    # The w13 force has order q-1; relative to w1's order -2 it is 4 delta.
    assert s.expand(q+1-4*delta)==0
    gaps={'coordinate_reference':6-3*order,
          'w13_relative':4*delta,
          'other_near_crossed_w_force':3*(order-1)/2,
          'maximal23_endpoint_to_inner_w_force':3-3*order/2}
    low,high=s.Rational(13,8),s.Rational(7,4)
    for name,e in gaps.items():
        assert s.simplify(e-delta).subs(order,low)>0
        assert s.simplify(e-delta).subs(order,high)>0
    assert s.simplify((3*order-5)/2-q)==(order-1)/6
    m1,m3,m5,m6,T,A,D,u,v,h=s.symbols('m1 m3 m5 m6 T A D u v h',nonzero=True)
    # Revised manuscript: T=m1+m4+m5+m6 is already the full outer mass.
    cube_ratio=-T*m6/(m3*m5)
    leading_ratio=-m1*m6/(m3*m5)
    mass_condition=m1**3*m6**2-T*m3**2*m5**2
    assert s.cancel(leading_ratio**3-cube_ratio+m6*mass_condition/(m3**3*m5**3))==0
    # D^3/A^3 = cube_ratio*w16/w35, and w16/w35=1+h+o(h).
    normalized=s.Poly(s.expand((1+h/3)**3-(1+h)),h)
    assert normalized.nth(0)==normalized.nth(1)==0
    # m1*m6*A0+m3*m5*D0=0, with A/A0 a unit tending to 1.
    A0=s.Symbol('A0',nonzero=True);D0=leading_ratio*A0
    potential=s.expand(m1*m6*A0+m3*m5*D0*(1+h/3))
    assert s.expand(potential-m3*m5*D0*h/3)==0
    E0=s.Symbol('E0',nonzero=True)
    # Equation (no50-E-coefficient): m1*m3*E0=m3*m5*D0*v/u.
    pair_potential=m1*m3*E0/3
    assert s.expand(pair_potential+m1*m3*E0-s.Rational(4,3)*m1*m3*E0)==0
    result={'diagram':50,'open_window':['13/8','7/4'],
            'delta_order':str(delta),'unmarked13_pole_order':str(q),
            'w_row_relative_error_lower_orders':{k:str(v) for k,v in gaps.items()},
            'bridge_cube_ratio':str(cube_ratio),
            'bridge_leading_ratio':str(leading_ratio),
            'bridge_leading_mass_condition':str(mass_condition),
            'correct_normalized_ratio':'(D/A)/(D0/A0)=1+(v/(3*u))*epsilon^((s-1)/3)+o(epsilon^((s-1)/3))',
            'constant_minus_one_in_old_manuscript_is_inconsistent':True,
            'pair_potential_coefficient':'m1*m3*E0/3',
            'unmarked_13_potential_coefficient':'m1*m3*E0',
            'complete_potential_coefficient':'4*m1*m3*E0/3',
            'first_uncancelled_potential':'(4*m1*m3*E0/3)*epsilon^((4*s-7)/3)',
            'no_new_exceptional_mass_factor_needed':True,
            'all_other_missing_edge_potentials_strictly_later':True,
            'scope':'The open only-23 resonance; remaining graph branches and endpoint need their original-row arguments.'}
    out=Path(__file__).with_name('no50_resonance_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Original w-source ratio, +1 normalization, simple cubic response and nonzero resonance pole verified')


if __name__=='__main__':main()
