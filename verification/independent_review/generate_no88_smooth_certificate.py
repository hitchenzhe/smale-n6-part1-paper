"""Reconstruct both affine orientations, masses and H at No. 88's F1012 point."""
from itertools import combinations
from pathlib import Path
import json
import sympy as s
from exact_geometry import EDGES,cm_inverse_root_numerator
from finite_field import CyclotomicField


def main():
    f=CyclotomicField(101)
    roots=tuple(map(f.c,(2,3,4,75,(0,37),(0,18))))
    ds={e:f.pow(r,-2) for e,r in zip(EDGES,roots)}
    q={e:f.sub(1,f.pow(r,3)) for e,r in zip(EDGES,roots)}
    assert all(a!=f.zero for a in q.values())
    A,B,C,D,E,F=(ds[e] for e in EDGES)
    b=f.sub(f.add(A,B),D)
    discriminant=f.sub(f.mul(b,b),f.mul(4,f.mul(A,B)))
    sq=f.sqrt(discriminant)
    assert sq!=f.zero
    records=[]
    for sign in (1,-1):
        z3=f.div(f.add(b,f.mul(sign,sq)),f.mul(2,A));w3=f.div(B,z3)
        rhs1=f.sub(f.add(A,C),E);rhs2=f.sub(f.add(B,C),F)
        det=f.sub(f.mul(A,z3),w3)
        z4=f.div(f.sub(f.mul(rhs1,z3),rhs2),det)
        w4=f.div(f.sub(f.mul(A,rhs2),f.mul(w3,rhs1)),det)
        z=(f.zero,f.one,z3,z4);w=(f.zero,A,w3,w4)
        for ij in EDGES:
            i,j=ij
            assert f.mul(f.sub(z[i],z[j]),f.sub(w[i],w[j]))==ds[ij]
        areas=[]
        for omitted in range(4):
            keep=[i for i in range(4) if i!=omitted]
            areas.append(f.mul((-1)**omitted,f.det([[1,z[i],w[i]] for i in keep])))
        assert all(a!=f.zero for a in areas)
        m=[f.one]
        for i in range(1,4):
            k=next(k for k in range(4) if k not in (i,0))
            q0k=q[tuple(sorted((0,k)))];qik=q[tuple(sorted((i,k)))]
            m.append(f.div(f.mul(areas[i],q0k),f.mul(areas[0],qik)))
        M=f.sum(m);m=[f.div(a,M) for a in m]
        Zc=f.sum(f.mul(a,c) for a,c in zip(m,z));Wc=f.sum(f.mul(a,c) for a,c in zip(m,w))
        Z=[f.sub(a,Zc) for a in z];W=[f.sub(a,Wc) for a in w]
        assert all(a!=f.zero for a in W)
        residual=[]
        for i in range(4):
            for pos in (Z,W):
                r=f.sum(f.mul(m[j],f.mul(f.pow(roots[EDGES.index(tuple(sorted((i,j))))],3),f.sub(pos[i],pos[j]))) for j in range(4) if j!=i)
                residual.append(f.sub(pos[i],r))
        assert all(r==f.zero for r in residual),residual
        ui=[f.div(f.mul(a,a),b) for a,b in zip(m,W)]
        e1=f.sum(ui);e2=f.sum(f.mul(ui[i],ui[j]) for i,j in combinations(range(4),2))
        e4=f.one
        for a in ui:e4=f.mul(e4,a)
        H=f.sub(f.pow(f.sub(f.pow(e1,2),f.mul(4,e2)),2),f.mul(64,e4))
        assert H!=f.zero
        records.append({'orientation_sign':sign,'affine_z':[list(a) for a in z],
                        'affine_w':[list(a) for a in w],'normalized_masses':[list(a) for a in m],
                        'centered_W':[list(a) for a in W],'H':list(H),
                        'eight_original_force_row_residuals':[list(a) for a in residual]})
    x,y,z,t,v,w=s.symbols('x y z t v w');variables=(x,y,z,t,v,w)
    cm=cm_inverse_root_numerator(variables)
    e1=(1-y**3)*(v**3-1)+(1-z**3)*(1-t**3)
    e2=(1-x**3)*(w**3-1)+(1-z**3)*(1-t**3)
    jac=[[f.eval(s.Poly(s.diff(expr,a),*variables),roots) for a in (t,v,w)] for expr in (cm.as_expr(),e1,e2)]
    jacdet=f.det(jac)
    assert jacdet!=f.zero
    # The manuscript uses the primitive numerator after fixing x,y,z=2,3,4.
    # Its integer content is 32, hence that Jacobian is the generic one /32.
    fixed=s.Poly(cm.as_expr().subs({x:2,y:3,z:4}),t,v,w)
    fixed_content=int(fixed.content())
    primitive_jac=f.div(jacdet,fixed_content)
    assert primitive_jac==(97,63),(fixed_content,primitive_jac)
    assert [r['H'] for r in records]==[[93,99],[99,57]]
    result={'diagram':88,'field':'F101[omega]/(omega^2+omega+1)',
            'affine_normalization':'z1=w1=0, z2=1, w2=D12^-2, sum(m_i)=1',
            'distance_ties_verified':True,'both_force_systems_verified':True,
            'quadratic_orientation_discriminant':list(discriminant),
            'CM_recovery_jacobian_rows':['CM','opposite_y','opposite_x'],
            'CM_recovery_jacobian_columns':['t','v','w'],
            'generic_cleared_numerator_jacobian':list(jacdet),
            'specialized_numerator_integer_content':fixed_content,
            'primitive_specialized_numerator_jacobian':list(primitive_jac),
            'orientations':records}
    out=Path(__file__).with_name('no88_smooth_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'H_values':[r['H'] for r in records],
                      'generic_jacobian':list(jacdet),'primitive_specialized_jacobian':list(primitive_jac)},indent=2))
    print(f'Saved {out}')


if __name__=='__main__':main()
