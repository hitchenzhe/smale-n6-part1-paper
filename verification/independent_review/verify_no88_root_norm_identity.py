"""Exact integer-polynomial proof of the No. 88 sixteen-sign root norm."""
from __future__ import annotations

from itertools import combinations, product
import hashlib
import json
from pathlib import Path


def add(a,b):
    out=dict(a)
    for m,c in b.items():
        out[m]=out.get(m,0)+c
        if out[m]==0:del out[m]
    return out


def scale(a,c):return {m:c*v for m,v in a.items() if c*v}


def mul(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():
            z=tuple(i+j for i,j in zip(x,y));out[z]=out.get(z,0)+c*d
    return {m:c for m,c in out.items() if c}


def power(a,n):
    out={(0,0,0,0):1}
    while n:
        if n&1:out=mul(out,a)
        a=mul(a,a);n>>=1
    return out


def canonical(poly):
    return sorted([[*m,c] for m,c in poly.items()])


def main():
    roots=[{tuple(1 if i==j else 0 for i in range(4)):1} for j in range(4)]
    signed_product={(0,0,0,0):1}
    for signs in product((-1,1),repeat=4):
        linear={}
        for sign,r in zip(signs,roots):linear=add(linear,scale(r,sign))
        signed_product=mul(signed_product,linear)

    # u_i = r_i^2; construct e_1,e_2,e_4 exactly in the r variables.
    u=[power(r,2) for r in roots]
    e1={};e2={};e4={(0,0,0,0):1}
    for item in u:e1=add(e1,item);e4=mul(e4,item)
    for i,j in combinations(range(4),2):e2=add(e2,mul(u[i],u[j]))
    h=add(power(add(power(e1,2),scale(e2,-4)),2),scale(e4,-64))
    h_squared=power(h,2)
    if signed_product!=h_squared:
        raise AssertionError("the displayed H formula does not square to the 16-sign product")
    encoded=canonical(signed_product)
    digest=hashlib.sha256(json.dumps(encoded,separators=(",",":")).encode()).hexdigest()
    record={"source":"direct expansion of the No. 88 displayed root norm",
            "ring":"Z[r1,r2,r3,r4]",
            "sign_choices":16,
            "identity":"product_{epsilon_i=+-1}(epsilon_1*r1+...+epsilon_4*r4) = ((e1^2-4e2)^2-64e4)^2",
            "H_sparse_term_count":len(h),"product_sparse_term_count":len(signed_product),
            "identity_verified_exactly":True,
            "expanded_product_sha256":digest,
            "expanded_product_sparse_terms":encoded}
    out=Path(__file__).with_name("no88_root_norm_certificate.json")
    out.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:record[k] for k in ("sign_choices","H_sparse_term_count","product_sparse_term_count","identity_verified_exactly","expanded_product_sha256")},indent=2))
    print(f"saved {out}")


if __name__=="__main__":main()
