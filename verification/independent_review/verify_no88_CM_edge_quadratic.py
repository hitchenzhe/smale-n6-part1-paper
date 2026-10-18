"""Exact determinant expansion for the No. 88 zero-star Cayley--Menger step."""
from __future__ import annotations

from itertools import combinations, permutations
import json
from pathlib import Path

N=6
EDGES=tuple(combinations(range(4),2))
INDEX={edge:i for i,edge in enumerate(EDGES)}


def add(a,b):
    out=dict(a)
    for m,c in b.items():
        out[m]=out.get(m,0)+c
        if not out[m]:del out[m]
    return out


def mul(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():
            z=tuple(i+j for i,j in zip(x,y));out[z]=out.get(z,0)+c*d
    return {m:c for m,c in out.items() if c}


def det(matrix):
    result={}
    for perm in permutations(range(len(matrix))):
        inversions=sum(perm[i]>perm[j] for i in range(len(perm)) for j in range(i+1,len(perm)))
        term={(0,)*N:(-1 if inversions%2 else 1)}
        for i,j in enumerate(perm):
            term=mul(term,matrix[i][j])
            if not term:break
        result=add(result,term)
    return result


def main():
    zero={(0,)*N:0};one={(0,)*N:1}
    cm=[[dict(zero) for _ in range(5)] for _ in range(5)]
    for i in range(4):cm[0][i+1]=one;cm[i+1][0]=one
    for (i,j),k in INDEX.items():
        mon=[0]*N;mon[k]=1
        cm[i+1][j+1]=cm[j+1][i+1]={tuple(mon):1}
    polynomial=det(cm)
    records=[]
    for edge,k in INDEX.items():
        # The coefficient of x_edge^2 is a polynomial in the remaining
        # squared distances; the Cayley--Menger identity predicts -2 times
        # the variable on the opposite edge.
        coeff={}
        for exponents,c in polynomial.items():
            if exponents[k]==2:
                reduced=list(exponents);reduced[k]=0
                coeff[tuple(reduced)]=c
        other=tuple(i for i in range(4) if i not in edge)
        opposite=(other[0],other[1]);expected=[0]*N;expected[INDEX[opposite]]=1
        expected_poly={tuple(expected):-2}
        if coeff!=expected_poly:
            raise AssertionError((edge,coeff,expected_poly))
        records.append({"free_edge_zero_based":list(edge),
                        "opposite_edge_zero_based":list(opposite),
                        "coefficient_of_free_squared_distance_squared":"-2 times opposite squared distance"})
    record={"source":"direct permutation expansion of the four-point Cayley--Menger determinant",
            "ring":"Z[s_12,s_13,s_14,s_23,s_24,s_34] where s_ij is a squared distance",
            "determinant_sparse_term_count":len(polynomial),
            "each_edge_degree_two":all(max(m[k] for m in polynomial)==2 for k in range(N)),
            "six_opposite_edge_leading_coefficients_verified":len(records)==6,
            "records":records}
    out=Path(__file__).with_name("no88_CM_edge_quadratic_certificate.json")
    out.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:record[k] for k in ("determinant_sparse_term_count","each_edge_degree_two","six_opposite_edge_leading_coefficients_verified")},indent=2))
    print(f"saved {out}")


if __name__=="__main__":main()
