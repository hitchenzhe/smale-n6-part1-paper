"""Small exact geometric constructions shared by independent generators."""
from itertools import combinations,permutations
import sympy as s

EDGES=tuple(combinations(range(4),2))


def cm_squared_distance_polynomial():
    """Permutation expansion of the 5x5 four-point CM determinant."""
    d=s.symbols('s12 s13 s14 s23 s24 s34')
    A=s.zeros(5)
    for i in range(1,5):A[0,i]=A[i,0]=1
    for (i,j),q in zip(EDGES,d):A[i+1,j+1]=A[j+1,i+1]=q
    terms=[]
    for p in permutations(range(5)):
        n=sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))
        terms.append((-1)**n*s.prod(A[i,p[i]] for i in range(5)))
    return s.Poly(sum(terms),*d,domain=s.ZZ)


def cm_inverse_root_numerator(roots):
    """Clear product(D_ij^4) from CM(D_ij^-2), without rational det expansion."""
    squared=cm_squared_distance_polynomial()
    return s.Poly.from_dict({tuple(4-2*k for k in m):c for m,c in squared.terms()},
                           roots,domain=s.ZZ)
