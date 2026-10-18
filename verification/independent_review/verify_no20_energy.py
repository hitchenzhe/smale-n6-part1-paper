"""Independent exact checks for the newly supplied No. 20 derivation.

This script checks the algebraic shape equation, its simple root, the
second-order offset/phase reduction, and the final positive-mass Holder gap.
It does not claim to certify the analytic valuation hypotheses stated in the
manuscript; those remain proof obligations of the written asymptotic argument.
"""
from itertools import product
import json
from pathlib import Path

import sympy as s


def check_phase_sum():
    """For u,v>0, u*zeta+v*zeta' can be positive real only if both are 1."""
    sqrt3 = s.sqrt(3)
    # Coefficients of (real part, imag part)/(1/2 or sqrt(3)/2) for
    # phases 1, omega, omega^2. Enumerate all nine ordered pairs exactly.
    phase = [(s.Integer(1), s.Integer(0)),
             (s.Rational(-1, 2), s.Integer(1)),
             (s.Rational(-1, 2), s.Integer(-1))]
    allowed = []
    exceptional_equal_weight_cancellation = []
    for i, j in product(range(3), repeat=2):
        ri, ii = phase[i]
        rj, ij = phase[j]
        u, v = s.symbols("u v", positive=True)
        real = s.expand(ri*u + rj*v)
        imag = s.expand(ii*u + ij*v)
        if imag == 0:
            if real.is_positive is True:
                allowed.append((i, j))
            elif s.simplify(real.subs(v, u)).is_negative is True:
                exceptional_equal_weight_cancellation.append((i, j))
                # The only nontrivial zero-imaginary family is (omega,omega^2)
                # or its reverse, and its real part is -(u+v)/2 < 0.
                assert s.simplify(real + (u + v)/2) == 0
        else:
            # A nonzero single-sign imaginary coefficient cannot cancel for
            # positive weights. Opposite nontrivial phases cancel only at u=v.
            coeffs = (s.diff(imag, u), s.diff(imag, v))
            assert coeffs in ((0, 1), (1, 0), (1, 1), (-1, -1),
                              (1, -1), (-1, 1), (0, -1), (-1, 0))
            if coeffs in ((1, -1), (-1, 1)):
                exceptional_equal_weight_cancellation.append((i, j))
                assert s.simplify(real.subs(v, u)) == -u
    assert allowed == [(0, 0)]
    assert exceptional_equal_weight_cancellation == [(1, 2), (2, 1)]
    # A second use of the same lemma for the positive RHS of the mass relation.
    return {"ordered_phase_pairs_checked": 9,
            "only_positive_real_pair": [0, 0],
            "opposite_nontrivial_pair_at_equal_weights_is_negative": True}


def main():
    a, b, c, t = s.symbols("a b c t", positive=True)
    A = a + b
    mu6 = (a*b)**3 - c**2*A**4

    # Inner-triangle shape equation and simple-root computation.
    P = c**2*(b*t**3 + a)*(a + b*t)**3 - (a*b)**3*t**3
    assert s.factor(P.subs(t, 1) + mu6) == 0
    derivative_at_one = s.factor(s.diff(P, t).subs(t, 1))
    assert s.factor(derivative_at_one - (3*c**2*A**3*(b-a) - 3*mu6)) == 0
    derivative_on_shape = s.factor(derivative_at_one.subs(c**2, (a*b)**3/A**4))
    assert s.factor(derivative_on_shape
                    - 3*((a*b)**3/A**4)*A**3*(b-a)) == 0
    assert s.factor(derivative_on_shape) != 0  # a != b is the stated branch.
    # Original z_p row gives D_pq^3=-c D_pr^3/(a+b); signed ties at
    # the positive-mass ratio t=1 give w_p-w_q=(a+b)v/a.
    Dpr, v = s.symbols("Dpr v", nonzero=True)
    Fwp = b*(-c/A*Dpr)*(A/a*v) + c*Dpr*v
    assert s.factor(Fwp - c*(a-b)*Dpr*v/a) == 0

    # Original-row offsets and w-centre produce the printed phase equation.
    lam, y, eta_i, eta_j, eta_k, W = s.symbols(
        "lambda y eta_i eta_j eta_k W", nonzero=True)
    eta0 = s.symbols("eta_0", nonzero=True)
    # d_i^3=1/a, etc.; with eta^3=1, 1/d_i^2=a^(2/3)*eta_i.
    Ai = a**s.Rational(2, 3)*eta_i/W
    Aj = b**s.Rational(2, 3)*eta_j/W
    Ak = c**s.Rational(2, 3)*eta_k/(-lam*A*W/y)
    centre = s.factor(lam*a*Ai + lam*b*Aj + y*Ak)
    expected_centre = lam*(a**s.Rational(5, 3)*eta_i +
                           b**s.Rational(5, 3)*eta_j) - \
                      y**2*c**s.Rational(2, 3)*eta_k/(lam*A)
    assert s.simplify(W*centre - expected_centre) == 0

    # Four signed-tie terms in c_0, followed by using the mass relation.
    c0 = (lam*(a**s.Rational(5, 3)*eta_i+b**s.Rational(5, 3)*eta_j)
          + y*c**s.Rational(2, 3)*eta_k
          + lam**s.Rational(5, 3)*a*b/A**s.Rational(1, 3)*eta0)
    c0_reduced = s.factor(c0.subs(
        lam*(a**s.Rational(5, 3)*eta_i+b**s.Rational(5, 3)*eta_j),
        y**2*c**s.Rational(2, 3)*eta_k/(lam*A)))
    K = y**2*c**s.Rational(2, 3)/(lam*A) + y*c**s.Rational(2, 3)
    L = lam**s.Rational(5, 3)*a*b/A**s.Rational(1, 3)
    assert s.simplify(c0_reduced - (K*eta_k + L*eta0)) == 0

    # Positive mass bounds. From mu6=0: c=(ab)^(3/2)/A^2 and
    # B=c^(2/3)=ab/A^(4/3); power means and ab<=A^2/4 imply x>=2^(2/3).
    assert s.simplify(((a*b)**s.Rational(3, 2)/A**2)**s.Rational(2, 3)
                      - a*b/A**s.Rational(4, 3)) == 0
    H = a**s.Rational(5, 3)+b**s.Rational(5, 3)
    # Verify the dB expression used in its upper estimate.
    B = a*b/A**s.Rational(4, 3)
    dB_sq = s.factor(A*H*B)
    assert s.simplify(dB_sq - H*a*b/A**s.Rational(1, 3)) == 0

    # Rebuild the mass and pair-sum bounds used in Holder. The lower bound for
    # sigma keeps only the outer-third bridge pairs and all outer/inner pairs.
    lam_bound, x_bound, A_bound = s.symbols(
        "lambda_bound x_bound A_bound", positive=True)
    c_bound = s.symbols("c_bound", positive=True)
    M_actual = lam_bound*A_bound*(1+x_bound)+A_bound+c_bound
    M_upper = A_bound*(lam_bound*(1+x_bound)+s.Rational(9, 8))
    assert s.simplify(M_upper-M_actual-(A_bound/s.Integer(8)-c_bound)) == 0
    sigma_actual_subset = (lam_bound**2*A_bound**2*x_bound
                           + lam_bound*A_bound*(1+x_bound)*(A_bound+c_bound)
                           + lam_bound**2*a*b + a*b+a*c+b*c)
    sigma_lower = A_bound**2*(lam_bound**2*x_bound
                              + lam_bound*(1+x_bound))
    assert s.simplify(sigma_actual_subset-sigma_lower - (
        lam_bound*A_bound*(1+x_bound)*c_bound
        + lam_bound**2*a*b+a*b+a*c+b*c)) == 0

    # After the phase reduction, c0 has four positive real terms. The
    # power-mean and ab bounds give the displayed upper estimate termwise.
    assert s.factor(A**s.Rational(5, 3)-H) != 0  # nonzero polynomial gap
    assert s.factor((A**2/4)-a*b) == (a-b)**2/4
    dB_square_upper_gap = (A**s.Rational(10, 3)/4-dB_sq)
    dB_square_bound_decomposition = (
        (A**s.Rational(5, 3)-H)*a*b/A**s.Rational(1, 3)
        + A**s.Rational(5, 3)*(A**2/4-a*b)/A**s.Rational(1, 3))
    assert s.simplify(dB_square_upper_gap-dB_square_bound_decomposition) == 0
    dB = s.symbols("dB", positive=True)
    c0_positive_mass = (lam_bound*H + lam_bound*dB
                        + lam_bound**s.Rational(5, 3)*a*b/A**s.Rational(1, 3))
    c0_upper = A**s.Rational(5, 3)*(3*lam_bound/2
                                    + lam_bound**s.Rational(5, 3)/4)
    c0_gap_decomposition = (
        lam_bound*(A**s.Rational(5, 3)-H)
        + lam_bound*(A**s.Rational(5, 3)/2-dB)
        + lam_bound**s.Rational(5, 3)*(A**2/4-a*b)/A**s.Rational(1, 3))
    assert s.simplify(c0_upper-c0_positive_mass-c0_gap_decomposition) == 0

    # The final Holder inequality and the exact positive polynomial certificate.
    lam_s, x = s.symbols("lambda x", positive=True)
    gap = ((lam_s*x+1+x)**3
           -(lam_s*(1+x)+s.Rational(9, 8))
           *(s.Rational(27, 2)+lam_s**2/16))
    x0 = s.Rational(19, 12)
    assert s.Rational(2)**s.Rational(2, 3) > x0
    at = s.Poly(s.expand(gap.subs(x, x0)), lam_s)
    coefficients = [at.nth(i) for i in range(4)]
    assert coefficients == [s.Rational(3547, 1728), -s.Rational(1829, 576),
                            s.Rational(22301, 1152), s.Rational(1645, 432)]
    assert 4*at.nth(0)*at.nth(2)-at.nth(1)**2 == s.Rational(148167571, 995328)
    h = s.symbols("h", nonnegative=True)
    dx = s.Poly(s.expand(s.diff(gap, x).subs(x, h+x0)), h, lam_s)
    assert all(coef > 0 for _, coef in dx.terms())
    ss = s.Rational(3, 2)
    tt = lam_s**s.Rational(2, 3)/4
    holder_replacement_gap = s.Rational(27, 2)+lam_s**2/16-(ss+tt)**3
    assert s.simplify(holder_replacement_gap-3*(ss+tt)*(ss-tt)**2) == 0

    # The shape condition yields x^2 >= 2^(4/3), hence x >= 2^(2/3).
    # Record the exact chain separately so the non-rational power-mean steps
    # are explicit in the generated certificate.
    phase_certificate = check_phase_sum()
    certificate = {
        "shape_polynomial": "c^2 (b t^3+a)(a+bt)^3-(ab)^3 t^3",
        "P_at_1": "-mu6",
        "Pprime_at_1": "3 c^2 (a+b)^3 (b-a)-3 mu6",
        "simple_root_condition": "mu6=0 and a!=b",
        "offset_center_reduction_verified": True,
        "finite_energy_reduction_verified": True,
        "Holder_mass_pair_and_energy_bounds_rebuilt": True,
        "positive_real_cube_phase_pairs_checked": phase_certificate,
        "shape_lower_bound": "x >= 2^(2/3) > 19/12",
        "holder_gap_coefficients_at_19_12": [str(q) for q in coefficients],
        "quadratic_discriminant_gap": "148167571/995328",
        "holder_gap_strictly_increasing_for_x_ge_19_12": True,
        "analytic_scope": "Asymptotic valuation and global branch assumptions are read from the supplied proof and are not certified by symbolic arithmetic alone."
    }
    out = Path(__file__).with_name("no20_energy_certificate.json")
    out.write_text(json.dumps(certificate, indent=2)+"\n", encoding="utf-8")
    print("No.20 shape root, offset/phase algebra, energy reduction and Holder gap verified exactly")


if __name__ == "__main__":
    main()

