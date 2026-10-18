"""Rebuild the No. 88 finite-field degree-90 norm from the displayed paper data.

Requires SymPy 1.14.  The companion standard-library verifier checks the saved
coefficient vector and its irreducibility over F_101 independently.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp
from math import gcd
from exact_geometry import cm_inverse_root_numerator


P = 101
X0, Y0, Z0 = 52, 26, 44
t, v, w = sp.symbols("t v w")
x, y, z = sp.symbols("x y z")


def build_cm_numerator():
    # Expand the 22-term squared-distance determinant, then substitute
    # reciprocal squared roots. This avoids a symbolic rational determinant.
    generic=build_general_cm_numerator()
    raw=sp.Poly(generic.as_expr().subs({x:X0,y:Y0,z:Z0}),t,v,w)
    constant_denominator=X0**4*Y0**4*Z0**4
    common=gcd(int(raw.content()),constant_denominator)
    poly=raw.exquo_ground(common)
    denominator=(constant_denominator//common)*t**4*v**4*w**4
    if poly.total_degree() != 10:
        raise AssertionError((poly.total_degree(), denominator))
    return poly.as_expr(), str(denominator)


def build_general_cm_numerator():
    return cm_inverse_root_numerator((x,y,z,t,v,w))


def f2_add(a, b):
    return ((a[0] + b[0]) % P, (a[1] + b[1]) % P)


def f2_mul(a, b):
    # omega^2 = -omega - 1.
    return ((a[0] * b[0] - a[1] * b[1]) % P,
            (a[0] * b[1] + a[1] * b[0] - a[1] * b[1]) % P)


def f2_pow(a, n):
    out = (1, 0)
    while n:
        if n & 1:
            out = f2_mul(out, a)
        a = f2_mul(a, a)
        n >>= 1
    return out


def eval_extension(poly, values):
    result = (0, 0)
    for monom, coeff in poly.terms():
        term = (int(coeff) % P, 0)
        for value, exponent in zip(values, monom):
            term = f2_mul(term, f2_pow(value, exponent))
        result = f2_add(result, term)
    return result


def main():
    cm, denominator = build_cm_numerator()
    cm_mod_poly = sp.Poly(cm, t, v, w, modulus=P)
    general_cm = build_general_cm_numerator()
    c_v = ((1 - Z0**3) * pow((1 - Y0**3) % P, -1, P)) % P
    c_w = ((1 - Z0**3) * pow((1 - X0**3) % P, -1, P)) % P
    # Opposite-edge equations give v^3=1-c_v(1-t^3),
    # w^3=1-c_w(1-t^3).
    f_v = v**3 - (1 - c_v * (1 - t**3))
    f_w = w**3 - (1 - c_w * (1 - t**3))
    cm_f = sp.Poly(cm, t, v, w, modulus=P).as_expr()
    first = sp.Poly(sp.resultant(cm_f, f_v, v),t,w,modulus=P).as_expr()
    norm = sp.resultant(first, f_w, w)
    norm_poly = sp.Poly(norm, t, modulus=P).monic()
    if norm_poly.degree() != 90:
        raise AssertionError((norm_poly.degree(), norm_poly.LC()))
    coeffs = [int(norm_poly.nth(i)) % P for i in range(91)]
    cm_terms = [
        {"x": int(monom[0]), "y": int(monom[1]), "z": int(monom[2]),
         "t": int(monom[3]), "v": int(monom[4]), "w": int(monom[5]),
         "coefficient": int(coeff) % P}
        for monom, coeff in sp.Poly(general_cm, x, y, z, t, v, w, modulus=P).terms()
    ]
    omega = (0, 1)
    shape_point = ((2, 0), (3, 0), (4, 0), (75, 0), f2_mul((37, 0), omega),
                   f2_mul((18, 0), omega))
    X, Y, Z, T, V, W = shape_point
    eq_y = f2_add(f2_mul(((1 - 3**3) % P, 0), f2_add(f2_pow(V, 3), (P - 1, 0))),
                  f2_mul(((1 - 4**3) % P, 0), ((1 - pow(75, 3, P)) % P, 0)))
    eq_x = f2_add(f2_mul(((1 - 2**3) % P, 0), f2_add(f2_pow(W, 3), (P - 1, 0))),
                  f2_mul(((1 - 4**3) % P, 0), ((1 - pow(75, 3, P)) % P, 0)))
    cm_at_point = eval_extension(general_cm, shape_point)
    if eq_y != (0, 0) or eq_x != (0, 0) or cm_at_point != (0, 0):
        raise AssertionError((eq_y, eq_x, cm_at_point))
    record = {
        "source": "independent reconstruction from derived_relations.tex equations (No. 88)",
        "field": "F_101",
        "base_specialization": {"x": X0, "y": Y0, "z": Z0},
        "cm_numerator_total_degree_t_v_w": 10,
        "cm_cleared_denominator": denominator,
        "cm_numerator_terms_mod_101": cm_terms,
        "smooth_shape_point_mod_1012": {
            "coordinates_in_1_omega_basis": [[2, 0], [3, 0], [4, 0],
                                                [75, 0], [0, 37], [0, 18]],
            "opposite_equation_residuals": [list(eq_y), list(eq_x)],
            "cayley_menger_residual": list(cm_at_point),
        },
        "v_cubic_recovery_coefficient": c_v,
        "w_cubic_recovery_coefficient": c_w,
        "resultant_order": ["Res_v(CM numerator, v recovery)", "Res_w(previous, w recovery)"],
        "norm_degree": norm_poly.degree(),
        "norm_leading_coefficient_before_monic_normalization": int(sp.Poly(norm, t, modulus=P).LC()),
        "norm_monic_coefficients_ascending_mod_101": coeffs,
        "monic_norm_irreducible": bool(norm_poly.is_irreducible),
    }
    if not record["monic_norm_irreducible"]:
        raise AssertionError("degree-90 fiber failed irreducibility over F_101")
    out = Path(__file__).with_name("no88_fiber_certificate.json")
    out.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: record[k] for k in (
        "field", "base_specialization", "cm_numerator_total_degree_t_v_w",
        "v_cubic_recovery_coefficient", "w_cubic_recovery_coefficient",
        "norm_degree", "norm_leading_coefficient_before_monic_normalization",
        "monic_norm_irreducible")}, indent=2))
    print(f"saved {out}")


if __name__ == "__main__":
    main()
