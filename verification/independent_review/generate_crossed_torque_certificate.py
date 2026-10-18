"""Original six-body cluster identities and No. 93/98 valuation transitions."""
from itertools import combinations
from pathlib import Path
import json
import sympy as s


def main():
    masses=s.symbols('m1:7');z=s.symbols('z1:7');w=s.symbols('w1:7')
    edges=tuple(combinations(range(6),2));roots=dict(zip(edges,s.symbols('d0:15')))
    Fz=[];Fw=[]
    for i in range(6):
        Fz.append(sum(masses[j]*roots[tuple(sorted((i,j)))]**3*(z[i]-z[j]) for j in range(6) if j!=i))
        Fw.append(sum(masses[j]*roots[tuple(sorted((i,j)))]**3*(w[i]-w[j]) for j in range(6) if j!=i))
    upper=(3,4,5);lower=(0,1,2)
    torque=sum(masses[i]*(z[i]*Fw[i]-w[i]*Fz[i]) for i in upper)
    crossed=sum(masses[i]*masses[j]*roots[(j,i)]**3*(w[i]*z[j]-z[i]*w[j]) for i in upper for j in lower)
    assert s.expand(torque-crossed)==0
    center_z=sum(masses[i]*Fz[i] for i in upper)
    cross_z=sum(masses[i]*masses[j]*roots[(j,i)]**3*(z[i]-z[j]) for i in upper for j in lower)
    assert s.expand(center_z-cross_z)==0
    center_w=sum(masses[i]*Fw[i] for i in upper)
    cross_w=sum(masses[i]*masses[j]*roots[(j,i)]**3*(w[i]-w[j]) for i in upper for j in lower)
    assert s.expand(center_w-cross_w)==0
    K,Z,W=s.symbols('K Z W')
    assert s.expand(-K**3*Z*W-K+K*(K*K*Z*W+1))==0
    alpha=s.symbols('alpha');order=1-alpha/2
    no93={
        'body4_z_position':alpha,'body4_crossed_inverse_distance':order,
        'body4_crossed_z_force':3*order+alpha,
        'body4_crossed_w_force':3*order-2,
        'body4_crossed_torque':3*order+alpha-2,
        'body56_crossed_torque':s.Integer(2),
        'relative_tie_corrections':[2-alpha,s.Integer(4),6-alpha],
    }
    assert s.simplify(no93['body4_crossed_torque']-order)==0
    assert s.simplify(no93['body4_crossed_z_force']-2)==order
    assert all(s.simplify(c)>0 if not c.free_symbols else (c.subs(alpha,-2)>0 and c.subs(alpha,2)>=0) for c in no93['relative_tie_corrections'])
    a,b,c,x=s.symbols('a b c x')
    assert s.expand((a+b+c)**2-3*(a*b+a*c+b*c)-((a-b)**2+(a-c)**2+(b-c)**2)/2)==0
    assert s.expand((1+x/3)**3-(1+x)-x*x/3-x**3/27)==0
    result={'diagrams':[93,98],
            'exact_original_cluster_torque_identity':True,
            'exact_mass_weighted_cluster_row_identities':True,
            'crossed_tie_leading_constraint':'K^2*Z*W=-1',
            'crossed_torque_leading_term':'-K^3*Z*W=K',
            'no93_domain':'-2<alpha<2; s=1-alpha/2 in (0,2)',
            'no93_orders':{k:([str(a) for a in v] if isinstance(v,list) else str(v)) for k,v in no93.items()},
            'no93_row_closure_argument':['w5,w6=O(epsilon^2) by their original rows',
                       'upper internal w gaps have order 2, hence w4=O(epsilon^2)',
                       'weighted lower z center is o(epsilon^2), lower z gaps order 2, hence all lower z=O(epsilon^2)',
                       'alpha>=2 contradicts D4j=o(1) via its original signed tie',
                       '4j torque order s is strictly earlier than 5j/6j torque order 2'],
            'no98_internal_orders_z_w_D':{'upper':[-2,2,0],'lower':[2,-2,0]},
            'no98_crossed_orders_z_w_D':[-2,-2,2],
            'no98_block_center_argument':'Mass-weighted original rows cancel internal forces; crossed forces have order 4, so each leading internal boundary is mass centered.',
            'no98_Lagrange_energy_identities_exact':True,
            'diagram_input':'The stated strokes/circles and maximal-stroke rules of the cited Chang-Chen classification. This file verifies the original-row algebra after that explicit input.'}
    out=Path(__file__).with_name('crossed_torque_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'cluster_torque_identity':True,'cluster_centers':True,
                      'no93_torque_order':str(order),'no98_torque_order':2,
                      'Lagrange_energy_identities':True},indent=2))
    print(f'Saved {out}')


if __name__=='__main__':main()
