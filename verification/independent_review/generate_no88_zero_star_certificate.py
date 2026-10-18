"""Universal No. 88 zero-star polynomial and mass-scale checks.

The same symbolic identities cover all four star vertices and all 27 cube
phases on their three spokes, provided the original distance/circle guards
hold. No finite-field sample substitutes for the universal identities.
"""
from itertools import combinations,product
from pathlib import Path
import json
import sympy as s
from exact_geometry import EDGES,cm_squared_distance_polynomial


def main():
    C=cm_squared_distance_polynomial()
    records=[]
    # For every star and every one of its three free nonstar edges, the
    # coefficient of that squared-distance variable squared is a nonzero
    # opposite spoke. This is phase independent.
    for star in range(4):
        nonstar=[j for j in range(4) if j!=star]
        for free in combinations(nonstar,2):
            k=EDGES.index(tuple(sorted(free)))
            opposite=tuple(sorted(set(range(4))-set(free)))
            coefficient=sum(c*s.prod(C.gens[j]**n for j,n in enumerate(m) if j!=k)
                            for m,c in C.terms() if m[k]==2)
            expected=-2*C.gens[EDGES.index(opposite)]
            assert s.expand(coefficient-expected)==0
            records.append({'star':star+1,'free_nonstar_edge':[j+1 for j in free],
                            'opposite_spoke':[j+1 for j in opposite],
                            'leading_coefficient':str(expected)})
    scale=s.symbols('scale');a,b,c=s.symbols('a b c');wa,wb,wc=s.symbols('wa wb wc')
    Cmass=a+b+c;K=a*wa+b*wb+c*wc
    total=1+scale*Cmass
    masses=[1/total,*[scale*q/total for q in (a,b,c)]]
    center=scale*K/total
    centered=[-center,wa-center,wb-center,wc-center]
    scaled_u=[s.cancel(scale*m*m/w) for m,w in zip(masses,centered)]
    limits=[s.cancel(q).subs(scale,0) for q in scaled_u]
    assert s.cancel(limits[0]+1/K)==0 and limits[1:]==[0,0,0]
    e1=sum(limits);e2=sum(limits[i]*limits[j] for i,j in combinations(range(4),2));e4=s.prod(limits)
    Hlimit=s.factor((e1**2-4*e2)**2-64*e4)
    assert s.cancel(Hlimit-1/K**4)==0
    phases=[list(p) for p in product(range(3),repeat=3)]
    result={'diagram':88,'star_vertices':4,'spoke_phase_triples':phases,
            'star_phase_labels':4*len(phases),'CM_quadratic_identities':records,
            'phase_independent_CM_nonidentity':True,
            'mass_normalization':'m_star=1/(1+scale*(a+b+c)); m_others=scale*(a,b,c)/(1+scale*(a+b+c))',
            'uncentered_star_W':'0','weighted_nonstar_W':'K=a*wa+b*wb+c*wc',
            'scaled_ui_limits':[str(q) for q in limits],
            'limit_scale4_H':str(Hlimit),
            'guards':['all three spoke distances nonzero',
                      'wa,wb,wc nonzero by spoke distance ties',
                      'K nonzero: otherwise the centered star W is zero, forbidden by its circle',
                      'total normalized mass denominator nonzero'],
            'scope':'Checks the universal nonidentity and free mass-scale cut; the rank-one stress supplies the stated fixed nonstar mass ratios.'}
    out=Path(__file__).with_name('no88_zero_star_certificate.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'star_phase_labels':result['star_phase_labels'],
                      'CM_identities':len(records),'limit_scale4_H':str(Hlimit)},indent=2))
    print(f'Saved {out}')


if __name__=='__main__':main()
