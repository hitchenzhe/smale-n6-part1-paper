"""Small exact checks tied to main_merged_revised.tex (2026-10-16).

This is an independent algebra check, not a proof of diagram branch coverage.
The JSON deliberately records contradictions as findings rather than hiding them
behind a successful process exit.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as s


def zero(expr):
    return s.cancel(s.together(expr)) == 0


def main():
    base = Path(__file__).parent
    root = base.parent.parent
    manuscript = root / "main_merged_revised.tex"
    source = manuscript.read_text(encoding="utf-8")
    results = {}

    # No. 66: every identity used in the new one-parameter reduction.
    a, b, c, d, u = s.symbols("a b c d u", nonzero=True)
    v = -(a + b * u) / c
    x = a * (v**3 - 1) / (b * (u**3 - v**3))
    y = a * (1 - u**3) / (c * (u**3 - v**3))
    w4, w5 = 1 / (u**2 * x), 1 / (v**2 * y)
    results["no66_spoke_and_centre"] = all((
        zero(a + b*u + c*v),
        zero(a + b*x + c*y),
        zero(a + b*u**3*x + c*v**3*y),
        zero(u**2*x*w4-1),
        zero(v**2*y*w5-1),
    ))
    T = b**3 + c**3
    results["no66_infinity_limits"] = all((
        zero(s.limit(x, u, s.oo) + a*b**2/T),
        zero(s.limit(y, u, s.oo) + a*c**2/T),
        zero(s.limit(u**2*w4, u, s.oo) + T/(a*b**2)),
        zero(s.limit(u**2*w5, u, s.oo) + T/(a*b**2)),
        zero(s.limit(u**2*(w4-w5), u, s.oo)),
    ))
    # The third tie gives r^2(x-y)(W4-W5)=1.  If W4-W5 vanishes
    # at order strictly beyond u^-2, the middle term of row 4 has
    # a pole stronger than u^3, while the other two terms do not.
    results["no66_difference_nonzero_generic"] = not zero(w4-w5)
    results["no66_difference_order_three"] = zero(
        s.limit(u**3 * (w4-w5), u, s.oo) - T / b**3
    )

    # Formal differential of U_C, with dD supplied by each signed tie.
    ms = s.symbols("m0:4")
    zs = s.symbols("z0:4")
    ws = s.symbols("w0:4")
    dms = s.symbols("dm0:4")
    dzs = s.symbols("dz0:4")
    dws = s.symbols("dw0:4")
    pairs = [(i, j) for i in range(4) for j in range(i+1, 4)]
    ds = {ij: s.Symbol(f"D{ij[0]}{ij[1]}") for ij in pairs}
    dU = 0
    Qs = [0] * 4
    fz = [0] * 4
    fw = [0] * 4
    for (i, j), Dij in ds.items():
        dz, dw = zs[i]-zs[j], ws[i]-ws[j]
        dD = -Dij**3 * (dw*(dzs[i]-dzs[j])+dz*(dws[i]-dws[j]))/2
        dU += (dms[i]*ms[j]+ms[i]*dms[j])*Dij+ms[i]*ms[j]*dD
        Qs[i] += ms[j]*Dij
        Qs[j] += ms[i]*Dij
        fz[i] += ms[j]*Dij**3*dz
        fz[j] -= ms[i]*Dij**3*dz
        fw[i] += ms[j]*Dij**3*dw
        fw[j] -= ms[i]*Dij**3*dw
    rhs = sum(Qs[i]*dms[i]-ms[i]*(fw[i]*dzs[i]+fz[i]*dws[i])/2 for i in range(4))
    results["no66_differential_identity"] = zero(dU-rhs)

    # No. 79: evaluate the modulus ordering of the *printed* bridges.
    theta, t = s.symbols("theta t", positive=True)
    p, q = (1-theta)/2, (1+theta)/2
    L = 3*theta**2*t**2 + 6*t
    d26_abs2 = (1+q*theta*t)**2 + q**2*L
    d35_abs2 = (1-p*theta*t)**2 + p**2*L
    gap = s.factor(d26_abs2-d35_abs2)
    results["no79_printed_modulus_gap"] = str(gap)
    results["no79_norm_gap_identity"] = zero(gap-4*theta*t*(2+theta**2*t))
    gamma = 1+q*theta*t
    G = s.factor(3*q**2*L-gamma**2)
    results["no79_large_asymmetry_G"] = str(G)
    results["no79_large_asymmetry_G_identity"] = zero(
        G-2*(1+theta)**2*theta**2*t**2-(1+theta)*(9+7*theta)*t/2+1
    )
    results["no79_large_asymmetry_G_corner"] = str(
        G.subs({theta:s.Rational(9,20),t:s.Rational(1,5)})
    )
    results["no79_large_asymmetry_G_corner_positive"] = (
        bool(G.subs({theta:s.Rational(9,20),t:s.Rational(1,5)}) > 0)
    )
    results["no79_large_asymmetry_G_corner_positive"] = (
        bool(G.subs({theta:s.Rational(9,20),t:s.Rational(1,5)}) > 0)
    )
    results["no79_weighted_ratio_present"] = (
        "r=-\\frac{\\beta d_{35}^{3}}{\\alpha d_{26}^{3}}" in source
    )
    results["no79_weighted_sector_identity_present"] = (
        "r=1-\\frac{K}{\\alpha d_{26}^3}" in source
    )

    # No. 52: verify the exact parameterization and two key positive gaps.
    lam, X, Y = s.symbols("lam X Y", positive=True)
    LL = 1+lam**3
    H5 = lam**2+lam**3+LL*(X+Y)+X*Y
    Q = LL+Y
    J = X+lam**2
    results["no52_H5_QJ_gap"] = str(s.factor(H5-Q*J))
    results["no52_gap_identity"] = zero(H5-Q*J-(1-lam**2)*(lam**3+Y)-lam**3*Y)
    results["no52_Cauchy_gap"] = str(s.factor(J**2-(X**2+lam**4*(X+1))))
    results["no52_Cauchy_identity"] = zero(J**2-X**2-lam**4*(X+1)-lam**2*(2-lam**2)*X)

    status = {k: v for k, v in results.items() if isinstance(v, bool)}
    assert all(status.values()), status
    output = {
        "manuscript": manuscript.name,
        "manuscript_sha256": hashlib.sha256(manuscript.read_bytes()).hexdigest(),
        "checks": results,
        "mathematical_status": "No.79 weighted sector algebra verified; original-row mapping and shape-energy bounds remain separate proof obligations",
        "scope": "Exact local algebra only; no branch-exhaustion or original-row coverage proof.",
    }
    out = base / "revised_core_claims_audit.json"
    out.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
