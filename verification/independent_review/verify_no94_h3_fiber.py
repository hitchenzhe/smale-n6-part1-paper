"""Replay the No. 94 raw H3=0 fiber elimination over F_17.

This independently derives all three raw signed-tie resultants from the
original row model, eliminates y against the H3 numerator, and checks the
degree-1608 eliminations and their degree-1120 common factor. The script does
not yet saturate the resulting fibers by all geometric/circle/collision
guards or certify the projective pre-specialization closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from itertools import combinations

from verify_no94_distinct_spoke_fiber import (
    CTX, RF, derive, quad_resultant, multidegree,
)
from flint import fq_default_ctx, fq_default_poly_ctx, fmpz_mod_poly_ctx


def eval_x_to_extension(poly, alpha, field_ctx, y_ctx):
    coeffs = {}
    for (ix, iy), coefficient in poly.to_dict().items():
        coeffs[iy] = coeffs.get(iy, field_ctx(0)) + field_ctx(int(coefficient)) * (alpha ** ix)
    degree = max(coeffs, default=0)
    return y_ctx([coeffs.get(j, field_ctx(0)) for j in range(degree + 1)])


def main():
    _, ties, guards, W = derive((1, 3, 5, 1))
    resultants = [quad_resultant(ties[i], ties[j])
                  for i, j in combinations(range(3), 2)]
    R = [r.n for r in resultants]
    degrees = [multidegree(f) for f in R]
    if degrees != [[102, 102], [96, 102], [102, 96]]:
        raise AssertionError(("pair resultant degrees", degrees))

    masses = [1, 3, 5]
    u = [-masses[i]**2 / W[i] for i in range(3)]
    H3 = (u[0] + u[1] + u[2])**2 - 4 * (
        u[0] * u[1] + u[0] * u[2] + u[1] * u[2]
    )
    Q = H3.n
    qdeg = multidegree(Q)
    if qdeg != [12, 12]:
        raise AssertionError(("H3 degree", qdeg))

    eliminants = [f.resultant(Q, "y") for f in R]
    edeg = [int(f.degrees()[0]) for f in eliminants]
    if edeg != [1608, 1608, 1608]:
        raise AssertionError(("elimination degrees", edeg))
    common = eliminants[0].gcd(eliminants[1]).gcd(eliminants[2])
    if common.degrees()[0] != 1120:
        raise AssertionError(("common x degree", common.degrees()))
    _, factors = common.factor()
    factor_records = sorted(
        [{"degree": int(f.degrees()[0]), "multiplicity": int(e),
          "coefficients_ascending": [int(f.to_dict().get((k, 0), 0)) % 17
                                      for k in range(int(f.degrees()[0]) + 1)]}
         for f, e in factors],
        key=lambda z: z["degree"],
    )
    expected = [(1, 360), (2, 360), (4, 10)]
    if [(f["degree"], f["multiplicity"]) for f in factor_records] != expected:
        raise AssertionError(factor_records)

    # On each irreducible x factor, compute the exact common y-fiber over its
    # finite extension and remove every explicitly inverted guard.
    fiber_records = []
    guard_names = ["spoke_1", "spoke_2", "spoke_3",
                   "circle_1", "circle_2", "circle_3",
                   "circle_difference_12", "circle_difference_13", "circle_difference_23",
                   "spoke_cube_difference_12", "spoke_cube_difference_13",
                   "spoke_cube_difference_23", "geometric_determinant"]
    for record in factor_records:
        degree = record["degree"]
        coeffs = record["coefficients_ascending"]
        if degree == 1:
            field_ctx = fq_default_ctx(17, 1, "z")
            alpha = field_ctx(-coeffs[0] * pow(coeffs[1], -1, 17) % 17)
        else:
            modulus = fmpz_mod_poly_ctx(17)(coeffs)
            field_ctx = fq_default_ctx(modulus=modulus, fq_type="FQ_NMOD")
            alpha = field_ctx.gen()
        y_ctx = fq_default_poly_ctx(field_ctx)
        fiber = eval_x_to_extension(R[0], alpha, field_ctx, y_ctx)
        for poly in [R[1], R[2], Q]:
            fiber = fiber.gcd(eval_x_to_extension(poly, alpha, field_ctx, y_ctx))
        fiber_degree = int(fiber.degree())
        if fiber_degree != {1: 12, 2: 12, 4: 1}[degree]:
            raise AssertionError((degree, fiber_degree))
        guard_product = y_ctx([field_ctx(1)])
        guard_hits = []
        for name, guard in zip(guard_names, guards):
            guard_gcd_degree = 0
            for guard_poly in (guard.n, guard.d):
                evaluated_guard = eval_x_to_extension(guard_poly, alpha, field_ctx, y_ctx)
                guard_product *= evaluated_guard
                if evaluated_guard != 0:
                    guard_gcd_degree = max(
                        guard_gcd_degree,
                        int(fiber.gcd(evaluated_guard).degree()),
                    )
            if guard_gcd_degree:
                guard_hits.append({"guard": name, "gcd_degree": guard_gcd_degree})
        guarded_part = fiber.gcd(guard_product)
        surviving_degree = fiber_degree - int(guarded_part.degree())
        if surviving_degree != 0:
            raise AssertionError((degree, fiber_degree, surviving_degree))
        fiber_records.append({"x_factor_degree": degree,
                              "common_y_fiber_degree": fiber_degree,
                              "guarded_saturated_degree": surviving_degree,
                              "guards_meeting_fiber": guard_hits})

    result = {
        "diagram": 94,
        "field": "F_17[x,y] via python-flint nmod_mpoly",
        "mass_witness_mod17": [1, 3, 5, 1],
        "pair_resultant_multidegrees_xy": degrees,
        "H3_numerator_multidegree_xy": qdeg,
        "three_y_elimination_degrees": edeg,
        "common_x_factor_degree": int(common.degrees()[0]),
        "common_x_factors": factor_records,
        "guard_saturated_y_fibers": fiber_records,
        "pairwise_guard_count": len(guards),
        "scope": "Replays the raw 102/96-degree H3 affine fiber elimination and exact guard saturation at the stated finite-field mass sample. The remaining projective boundary charts are separate."
    }
    out = Path(__file__).with_name("no94_h3_raw_fiber_certificate.json")
    out.write_text(json.dumps(result, indent=2, default=int) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, default=int))


if __name__ == "__main__":
    main()
