"""Replay the five No. 8 seven-product groups in main_merged_revised.tex.

The revised proof uses the actual order k of the three 1-to-T separations;
it only proves k <= h.  The older linear certificate substituted h for k
and therefore was not a certificate for this revised argument.
"""

from itertools import combinations
import json
from pathlib import Path

import sympy as s


def main():
    A, B, k = s.symbols("A B k", real=True)
    rho = A+B
    eta_q = (4+2*B)/3
    eta_t = (4+2*A)/3
    eta_1 = A+k
    edges = tuple(combinations(range(1, 7), 2))
    order = {e: rho for e in edges}
    for e in ((2,3), (2,6), (3,6)):
        order[e] = eta_q
    for e in ((4,5), (4,6), (5,6)):
        order[e] = eta_t
    for e in ((1,4), (1,5), (1,6)):
        order[e] = eta_1
    product = {
        e: s.expand(order[e] + sum(order[f] for f in edges if not set(e)&set(f)))
        for e in edges
    }
    groups = (
        (((1,2),(1,3)), 3*rho+eta_q+3*eta_t, (8-rho)/3,
         s.Rational(4,3)*(2+4*A+3*B)),
        (((1,4),(1,5),(4,5)), 2*rho+3*eta_q+eta_t+eta_1,
         3*eta_q+eta_t, 3*A+2*B+k),
        (((1,6),(4,6),(5,6)), 4*rho+eta_q+eta_t+eta_1,
         8*(1+rho)/3, 3*A+2*B+k),
        (((2,3),), eta_q+3*eta_t+3*eta_1,
         eta_q+3*eta_t-6*rho, 3*(3*A+2*B+k)),
        (((2,4),(2,5),(2,6),(3,4),(3,5),(3,6)),
         3*rho+eta_q+eta_t+2*eta_1,
         eta_q+eta_t-rho, 2*(3*A+2*B+k)),
    )
    seen = set()
    output_groups = []
    for members, formula, lower, difference in groups:
        assert all(e not in seen for e in members)
        seen.update(members)
        assert all(s.expand(product[e]-formula)==0 for e in members)
        assert s.expand(formula-lower-difference)==0
        output_groups.append({
            "pairs": ["".join(map(str,e)) for e in members],
            "order": str(s.expand(formula)),
            "lower_bound": str(s.expand(lower)),
            "nonnegative_slack": str(s.expand(difference)),
        })
    assert seen == set(edges)
    # The first slack is nonnegative by 4A+3B>=-2; the other
    # four by k>=-3A-2B.  rho>-1 follows from A<1 and
    # A+3rho=4A+3B>=-2.  eta_q,eta_t>0 follows from A,B>-2.
    out = {
        "diagram": 8,
        "manuscript": "main_merged_revised.tex",
        "position_supports": {"z":[1,4,5], "w":[1,2,3,4,5,6]},
        "assumptions_from_original_rows": [
            "rho=A+B<0", "A<1", "B<1", "4A+3B>=-2",
            "k>=-3A-2B", "k<=(4-A)/3",
        ],
        "deduced_positive_quantities": [
            "rho>-1", "eta_Q=(4+2B)/3>0", "eta_T=(4+2A)/3>0",
            "3A+2B+k>=0", "2+4A+3B>=0",
        ],
        "groups": output_groups,
        "seven_products_covered": 15,
        "all_group_lower_bounds_positive": True,
        "scope": "Exact complement-product algebra after the original-row order bounds in Lemma 5.9; those bounds are a separate analytic proof step.",
    }
    path = Path(__file__).with_name("no8_revised_order_certificate.json")
    path.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(f"Verified {len(groups)} groups and all {len(edges)} seven-products; saved {path}")


if __name__ == "__main__":
    main()
