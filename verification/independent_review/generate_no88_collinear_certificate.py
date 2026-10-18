"""Rebuild No. 88's 32 collinear sign classes over ZZ[u,v] and QQ.

This supersedes the incorrectly numbered old No. 66 modular script. Raw
minor gcds are computed over the integer polynomial ring, before any modular
reduction. Sample kernel, root norm and force slope are exact rationals.
"""
from __future__ import annotations

from itertools import combinations, product
import json
from pathlib import Path
import sympy as s

u,v=s.symbols('u v')
positions=(s.Integer(0),s.Integer(1),u,v)
edges=tuple(combinations(range(4),2))


def build(signs):
    signs=dict(zip(edges,(1,*signs)))
    gaps={ij:s.Poly((positions[ij[0]]-positions[ij[1]])**2,u,v) for ij in edges}
    L=s.prod(gaps.values())
    force=[]
    for i in range(4):
        row=[]
        for j in range(4):
            if i==j:row.append(s.Poly(0,u,v));continue
            ij=tuple(sorted((i,j)))
            row.append(L.exquo(gaps[ij])*signs[ij]*(1 if i<j else -1))
        force.append(row)
    rows=[[force[2][j]+force[0][j]*(u-1)-force[1][j]*u for j in range(4)],
          [force[3][j]+force[0][j]*(v-1)-force[1][j]*v for j in range(4)]]
    return rows,force,L


def root_H(masses,W):
    a=[s.cancel(m*m/w) for m,w in zip(masses,W)]
    e1=sum(a);e2=sum(a[i]*a[j] for i,j in combinations(range(4),2));e4=s.prod(a)
    return s.factor((e1**2-4*e2)**2-64*e4)


def run():
    records=[]
    for index,tail in enumerate(product((-1,1),repeat=5),1):
        rows,force,L=build(tail)
        minors=[rows[0][i]*rows[1][j]-rows[0][j]*rows[1][i] for i,j in edges]
        g=minors[0]
        for minor in minors[1:]:g=s.gcd(g,minor)
        content,primitive=g.primitive()
        sub={u:2,v:3}
        C=s.Matrix([[f.eval(sub) for f in row] for row in rows])
        reduced,pivots=C.rref()
        free=[j for j in range(4) if j not in pivots]
        basis=[]
        for j in free:
            b=s.zeros(4,1);b[j]=1
            for r,k in enumerate(pivots):b[k]=-reduced[r,j]
            basis.append(b)
        assert len(basis)==2,(tail,C)
        m=basis[0]-10*basis[1];M=sum(m)
        xi=(0,1,2,3)
        assert M!=0,(tail,'M=0')
        center=sum(m[i]*xi[i] for i in range(4))/M
        W=[s.factor(q-center) for q in xi]
        assert all(W),(tail,'a centered coordinate is zero')
        H=root_H(m,W)
        L0=L.eval(sub)
        F=s.Matrix([[f.eval(sub)/L0 for f in row] for row in force])*m
        slope=s.factor(F[1]-F[0])
        assert H!=0 and slope!=0,(tail,H,slope)
        record={
            'edge_order':['12','13','14','23','24','34'],
            'signs':[1,*tail],
            'raw_maximal_minor_gcd':str(g.as_expr()),
            'primitive_gcd':str(primitive.as_expr()),
            'gcd_is_constant':primitive.total_degree()==0,
            'gcd_content':int(content),
            'minor_total_degrees':[int(f.total_degree()) for f in minors],
            'sample_rank':len(pivots),
            'rref_pivots_zero_based':list(pivots),
            'kernel_basis':[[str(c) for c in b] for b in basis],
            'mass_vector_v0_minus_10v1':[str(c) for c in m],
            'total_mass':str(M),
            'centered_coordinates':[str(c) for c in W],
            'root_norm_H':str(H),
            'force_slope_F2_minus_F1':str(slope),
            'all_sample_mass_entries_nonzero':all(c!=0 for c in m),
        }
        records.append(record)
        print(f'class {index}/32: gcd={g.as_expr()}, rank={len(pivots)}, H and slope nonzero',flush=True)
    result={
        'diagram':88,'field':'ZZ[u,v] for minors; QQ for sample witnesses',
        'source':'No. 88 collinear force-slope equations in derived_relations.tex',
        'sympy_version':s.__version__,
        'classes':32,'full_signed_choices_with_global_reversal':64,
        'all_raw_minor_gcds_constant':all(r['gcd_is_constant'] for r in records),
        'all_sample_ranks_two':all(r['sample_rank']==2 for r in records),
        'records':records,
    }
    out=Path(__file__).with_name('no88_collinear_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f'Saved {out}; all raw minor gcds constant: {result["all_raw_minor_gcds_constant"]}',flush=True)
    return result


if __name__=='__main__':run()
