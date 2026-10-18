"""Exact symbolic checks for the No. 79 shape-to-energy bounds in the revised proof.

This checks the displayed sextic normal form, its wing-swap symmetry, the
equal-weight and unequal-weight bridge identities, and the rational positivity
certificate used for the universal t > 1/5 floor. It does not infer the
six-body valuation pattern or replace the geometric branch-coverage proof.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as s


def require_zero(expr: s.Expr, name: str) -> None:
    if s.cancel(s.factor(expr)) != 0:
        raise AssertionError(f"exact identity failed: {name}")


def main() -> None:
    a, b, d, r = s.symbols("a b d r", nonzero=True)
    u = (a + b) * r + d
    v = b * r + a + d
    w = a * (1 - r)
    tau, nu = w / u, w / v
    R = -1 + a**2 * nu / b**2 - a**2 * tau / d**2
    qhat = d**2 * a**2 * u * w - a**2 * b**2 * v * w - b**2 * d**2 * u * v
    f_hom = s.expand(qhat**3 - 27 * a**4 * b**4 * d**4 * u**2 * v**2 * w**2)
    require_zero(qhat - b**2 * d**2 * u * v * R, "Qhat=b^2 d^2 u v R")
    require_zero(
        f_hom / (b**6 * d**6 * u**3 * v**3)
        - (R**3 - 27 * a**4 * tau * nu / (b**2 * d**2)),
        "No. 79 normalized signed shape sextic (Eq. 30)",
    )
    # The positive real sector 0<r<1 forces the third triangle mass d>b.
    rr = s.symbols("rr", positive=True)
    Utri = (a + b) * rr + d
    Vtri = b * rr + a + d
    Wtri = a * (1 - rr)
    Qtri = a**2 * Wtri * (d**2 * Utri - b**2 * Vtri) - b**2 * d**2 * Utri * Vtri
    Ftri = Qtri**3 - 27 * a**4 * b**4 * d**4 * Utri**2 * Vtri**2 * Wtri**2
    require_zero(Vtri - Utri - Wtri, "triangle separation identity V-U=W")
    # Exact decomposition: if d<=b, both terms on the right are nonpositive,
    # and the second is strictly negative because W>0.
    require_zero(
        d**2 * Utri - b**2 * Vtri
        - ((d**2 - b**2) * Utri - b**2 * Wtri),
        "triangle sign decomposition",
    )
    assert s.factor((d**2 * Utri - b**2 * Vtri).subs(d, b) + b**2 * Wtri) == 0
    # Bridge weights alpha=lambda*b^2, beta=lambda*d^2 therefore obey p<n.
    lam = s.symbols("lam", positive=True)
    p_bridge, n_bridge = lam * b**2, lam * d**2
    require_zero(n_bridge - p_bridge - lam * (d**2 - b**2),
                 "positive-real bridge weight ordering")

    # Same-sign real bridge roots (r<0) are excluded by weighted Holder.
    px, ny, dx, dy = s.symbols("px ny dx dy", positive=True)
    holder_gap = s.expand((px + ny)**2 * (px * dx**3 + ny * dy**3)
                          - (px * dx + ny * dy)**3)
    require_zero(
        holder_gap - px * ny * (dx - dy)**2 * (
            (2 * px + ny) * dx + (px + 2 * ny) * dy
        ),
        "weighted Holder gap factorization",
    )
    alpha_w, beta_w, d26, d35 = s.symbols(
        "alpha_w beta_w d26 d35", nonzero=True
    )
    weighted_r = -beta_w * d35**3 / (alpha_w * d26**3)
    Kbridge = alpha_w * d26**3 + beta_w * d35**3
    require_zero(weighted_r - (1 - Kbridge / (alpha_w * d26**3)),
                 "weighted bridge ratio and K identity")
    qbridge = s.symbols("qbridge", positive=True)
    p0, n0, x0, y0 = s.symbols("p0 n0 x0 y0", positive=True)
    K0 = p0 * x0**3 - n0 * y0**3
    c0 = p0 * x0 - n0 * y0
    q0bridge = n0 * y0 / (p0 * x0)
    require_zero(
        p0**2 * K0 - c0**3
        - p0**3 * x0**3 * (
            3 * q0bridge * (1 - q0bridge)
            + (1 - p0**2 / n0**2) * q0bridge**3
        ),
        "unequal-sign bridge cubic gap",
    )

    # The final asymmetric sector: its bridge cubic has negative real part,
    # while r=1-K/(alpha*d26^3) then forces Re(r)>1.
    th, tt = s.symbols("th tt", positive=True)
    pp, qq = (1 - th) / 2, (1 + th) / 2
    LL = 3 * th**2 * tt**2 + 6 * tt
    gamma = 1 + qq * th * tt
    GG = 2 * (1 + th)**2 * th**2 * tt**2 + (1 + th) * (9 + 7 * th) * tt / 2 - 1
    require_zero(gamma**2 - 3 * qq**2 * LL + GG, "final-sector bridge real-part factor")
    if s.factor(s.diff(GG, th)).is_positive is not True or s.factor(s.diff(GG, tt)).is_positive is not True:
        raise AssertionError("final-sector G is not strictly increasing in both variables")
    if s.factor(GG.subs({th: s.Rational(9, 20), tt: s.Rational(1, 5)})) != s.Rational(1591621, 2000000):
        raise AssertionError("final-sector G boundary value changed")
    alpha_pos, beta_pos = s.symbols("alpha_pos beta_pos", positive=True)
    theta_from_weights = (beta_pos - alpha_pos) / (alpha_pos + beta_pos)
    require_zero(
        alpha_pos / beta_pos
        - (1 - theta_from_weights) / (1 + theta_from_weights),
        "triangle wing-product equals bridge-weight ratio",
    )

    A, B, z = s.symbols("A B z", nonzero=True)

    def normalized_sextic(x: s.Expr, y: s.Expr, root: s.Expr) -> s.Expr:
        U = (x + y) * root + 1
        V = y * root + x + 1
        W = x * (1 - root)
        Q = x**2 * U * W - x**2 * y**2 * V * W - y**2 * U * V
        return s.expand(Q**3 - 27 * x**4 * y**4 * U**2 * V**2 * W**2)

    require_zero(
        f_hom - d**18 * normalized_sextic(a / d, b / d, r),
        "homogeneous sextic agrees with the original triangle sextic",
    )
    require_zero(
        normalized_sextic(A, B, z)
        - z**6 * normalized_sextic(A / B, 1 / B, 1 / z) * B**18,
        "wing-swap/reciprocal shape identity",
    )

    # Equal bridge weights: exact Re(r), its monotonicity, and the endpoint.
    s0 = s.symbols("s0", real=True)
    re_r = -(1 - s0) * (1 - 14 * s0 + s0**2) / (1 + s0) ** 3
    require_zero(
        s.diff(re_r, s0) - 6 * (s0 - 3) * (3 * s0 - 1) / (s0 + 1) ** 4,
        "equal-weight Re(r) derivative",
    )
    if s.factor(re_r.subs(s0, s.Rational(1, 4))) != s.Rational(117, 125):
        raise AssertionError("equal-weight Re(r) endpoint changed")
    x = s.symbols("x", real=True)
    r_equal = -((1 - s.I * x) / (1 + s.I * x)) ** 3
    re_equal = s.expand_complex(r_equal).as_real_imag()[0]
    require_zero(re_equal - re_r.subs(s0, x**2), "equal-weight parametrized Re(r)")
    delta_equal = s.expand_complex((1 - r_equal) * (1 - s.conjugate(r_equal)))
    require_zero(
        delta_equal - 4 * (1 - 3 * x**2) ** 2 / (1 + x**2) ** 3,
        "equal-weight delta identity",
    )
    rb = 1 / r
    delta2 = (1 - r) * (1 - rb)
    H = a + b + d
    require_zero(
        u * ((a + b) * rb + d) - (H**2 - (a + b) * d * delta2),
        "equal-weight |u| norm identity",
    )
    require_zero(
        v * (b * rb + a + d) - (H**2 - b * (a + d) * delta2),
        "equal-weight |v| norm identity",
    )

    q0 = s.Rational(24, 25)
    if not s.sqrt(s.Rational(121, 125)) > q0:
        raise AssertionError("equal-weight denominator lower bound is too weak")
    if not s.Rational(64, 1) / (5 * s.sqrt(5)) < 6:
        raise AssertionError("equal-weight delta constant is not below 6")

    # The mass identity gives the sharp one-triangle lower bound used twice.
    t = s.symbols("t", nonnegative=True)
    mass_identity = (2 * t - 3) ** 2 * (27 * t**3 + 216 * t**2 + 76 * t + 12)
    require_zero(
        108 * (1 + t) ** 5 - 3125 * t**3 - mass_identity,
        "mass identity (Eq. 32)",
    )
    C0 = s.Rational(3125, 108)
    if not C0**2 > 625:
        raise AssertionError("mass lower bound does not dominate 625")

    # Unequal bridge weights: moments, squared moduli and weighted |r|.
    theta, t = s.symbols("theta t", positive=True)
    p, q = (1 - theta) / 2, (1 + theta) / 2
    L = 3 * theta**2 * t**2 + 6 * t
    N1 = (1 + q * theta * t) ** 2 + q**2 * L
    N2 = (1 - p * theta * t) ** 2 + p**2 * L
    require_zero(N1 - N2 - 4 * theta * t * (2 + theta**2 * t), "N1-N2")
    f0 = 1 - (1 - theta**2) * t * (3 + 2 * theta**2 * t) ** 2 / 2
    bridge_bracket = 4 + (1 - theta**2) * (
        4 * theta**2 * (1 + theta**2) * t**3
        + (15 * theta**2 + 9) * t**2
        + 12 * t
    )
    require_zero(
        q**2 * N2**3 - p**2 * N1**3 - theta * f0 * bridge_bracket / 4,
        "weighted bridge modulus comparison",
    )
    V_im = s.symbols("V_im", real=True)
    d26 = 1 + q * theta * t + s.I * q * V_im
    d35 = 1 - p * theta * t - s.I * p * V_im
    bridge_cubic = s.expand(p * d26**3 + q * d35**3)
    cubic_real = s.re(bridge_cubic)
    cubic_imag = s.im(bridge_cubic)
    require_zero(cubic_real.subs(V_im**2, L) - f0, "real weighted cubic moment f0")
    require_zero(
        cubic_imag - p * q * theta * V_im * (L - V_im**2),
        "imaginary weighted bridge cubic factor",
    )
    # Exact positive-coefficient chart proves f0(1/5)>(1+h)/(1+2h)^3.
    h, chart_v = s.symbols("h chart_v", positive=True)
    f0_at_floor = 1 - h**2 * (17 - 2 * h**2) ** 2 / 250
    mass_upper = (1 + h) / (1 + 2 * h) ** 3
    positive_gap = s.factor(250 * (1 + 2 * h) ** 3 * (f0_at_floor - mass_upper) / h)
    chart = s.Poly(
        s.cancel((1 + chart_v) ** 8 * positive_gap.subs(h, chart_v / (1 + chart_v))),
        chart_v,
    )
    gap_coefficients_descending = [int(coef) for coef in chart.all_coeffs()]
    expected_gap_coefficients_descending = [
        175, 9385, 53453, 129401, 167471, 125127, 54243, 12711, 1250
    ]
    if gap_coefficients_descending != expected_gap_coefficients_descending:
        raise AssertionError("the rational positivity chart changed")
    if not all(coef > 0 for coef in gap_coefficients_descending):
        raise AssertionError("the transformed t-floor gap is not positive")

    # Convert the two shape-root lower bounds and the Holder upper bound to Eq. (33).
    require_zero(
        16 / (p**4 * q**2)
        - 1024 / ((1 - theta) ** 4 * (1 + theta) ** 2),
        "N1-bound mass-ratio constant",
    )
    N1_rhs = s.factor(1024 / (q0**2 * C0**2))
    if N1_rhs != s.Rational(20736, 15625):
        raise AssertionError("Eq. (33) constant changed")

    result = {
        "source": "sections/矫正/1.png through 6.png",
        "status": "passed",
        "exact_checks": [
            "homogeneous signed triangle sextic and Eq. (30) normal form",
            "positive-real triangle forces d>b and hence p<n",
            "bridge mass-weighted ratio, K identity and signed cubic gap",
            "final sector Re(d26^3)<0 and strict Re(r)>1 bridge",
            "wing-swap reciprocal identity",
            "equal-weight Re(r) monotonicity and 117/125 endpoint",
            "equal-weight q0 and delta constants",
            "mass identity (32)",
            "unequal bridge N1-N2 identity",
            "unequal weighted bridge modulus comparison and cubic moment",
            "positive-coefficient rational chart for t>1/5",
            "Eq. (33) mass-ratio constant",
        ],
        "equal_weight_constants": {
            "q0": "24/25",
            "Re_r_at_s_1_over_4": "117/125",
            "delta_squared_upper": "16/125",
            "mass_factor_C0": "3125/108",
        },
        "t_floor_gap_coefficients_descending": gap_coefficients_descending,
        "unequal_weight_bound_constant": str(N1_rhs),
        "scope": (
            "Checks exact algebraic identities and rational constants in the new No. 79 "
            "shape-energy derivation. It does not independently regenerate the diagram, "
            "valuation branches, original-row-to-shape assignment, or global geometric cover."
        ),
    }
    out = Path(__file__).with_name("no79_shape_energy_certificate.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("No. 79 exact shape-sextic, equal-weight, and unequal-weight identities passed")


def test_no79_shape_energy_constants_replay_exactly() -> None:
    main()


if __name__ == "__main__":
    test_no79_shape_energy_constants_replay_exactly()
