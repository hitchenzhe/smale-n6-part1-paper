"""Recompute the No. 94 raw signed-tie resultant fiber over F_17.

This is deliberately a finite-field, denominator-explicit implementation
using python-flint multivariate polynomials. It derives the wing-edge cubes
from the original two wing w rows and one z-row difference, forms the three
quadratic signed ties in the common spoke scale, and checks the two bivariate
pair resultants. This is the raw 102-degree model; it does not certify the
manuscript's later degree-120 quotient/norm calculations.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from flint import nmod_mpoly_ctx

P = 17
CTX = nmod_mpoly_ctx.get(("x", "y"), P, "lex")
ONE = CTX.from_dict({(0, 0): 1})
ZERO = CTX.from_dict({})
X, Y = CTX.gens()


def scalar(n: int):
    return CTX.from_dict({(0, 0): n % P})


class RF:
    """Reduced rational function over F_17[x,y] using FLINT gcds."""
    __slots__ = ("n", "d")

    def __init__(self, n, d=ONE):
        if isinstance(n, int):
            n = scalar(n)
        if isinstance(d, int):
            d = scalar(d)
        if n == ZERO:
            self.n, self.d = ZERO, ONE
            return
        g = n.gcd(d)
        n, d = n // g, d // g
        lc = d.to_dict()[max(d.to_dict())]
        inv = pow(int(lc), -1, P)
        self.n, self.d = n * inv, d * inv

    @staticmethod
    def coerce(value):
        return value if isinstance(value, RF) else RF(scalar(int(value)))

    def __add__(self, other):
        other = RF.coerce(other)
        g = self.d.gcd(other.d)
        d1, d2 = self.d // g, other.d // g
        return RF(self.n * d2 + other.n * d1, d1 * other.d)

    __radd__ = __add__

    def __neg__(self):
        return RF(-self.n, self.d)

    def __sub__(self, other):
        return self + (-RF.coerce(other))

    def __rsub__(self, other):
        return RF.coerce(other) - self

    def __mul__(self, other):
        other = RF.coerce(other)
        g1, g2 = self.n.gcd(other.d), other.n.gcd(self.d)
        return RF((self.n // g1) * (other.n // g2),
                  (self.d // g2) * (other.d // g1))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = RF.coerce(other)
        if other.n == ZERO:
            raise ZeroDivisionError
        return self * RF(other.d, other.n)

    def __rtruediv__(self, other):
        return RF.coerce(other) / self

    def __pow__(self, exponent):
        if exponent < 0:
            return RF(self.d ** (-exponent), self.n ** (-exponent))
        return RF(self.n ** exponent, self.d ** exponent)


def det3(A):
    return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
            - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
            + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))


def solve3(A, b, det):
    result = []
    for j in range(3):
        B = [row[:] for row in A]
        for i in range(3):
            B[i][j] = b[i]
        result.append(det3(B) / det)
    return result


def derive(masses):
    a, b, c, d = [RF(scalar(v)) for v in masses]
    roots = [RF(scalar(1)), RF(X), RF(Y)]
    cubes = [z**3 for z in roots]
    W = [-(cubes[2] - cubes[1]) / a,
         (cubes[2] - cubes[0]) / b,
         -(cubes[1] - cubes[0]) / c]
    spoke = [1 / (z**2 * w) for z, w in zip(roots, W)]
    A = [[b * (W[0] - W[1]), c * (W[0] - W[2]), 0],
         [a * (W[1] - W[0]), 0, c * (W[1] - W[2])],
         [(a + b) * (spoke[0] - spoke[1]),
          c * (spoke[0] - spoke[2]), -c * (spoke[1] - spoke[2])]]
    rhs0 = [W[0], W[1], spoke[0] - spoke[1]]
    rhs1 = [-d * cubes[0] * W[0], -d * cubes[1] * W[1],
            -d * (cubes[0] * spoke[0] - cubes[1] * spoke[1])]
    determinant = det3(A)
    c0, c1 = solve3(A, rhs0, determinant), solve3(A, rhs1, determinant)
    ties = []
    for k, (i, j) in enumerate(combinations(range(3), 2)):
        factor = ((spoke[i] - spoke[j]) * (W[i] - W[j]))**3
        ties.append((factor * c0[k]**2,
                     2 * factor * c0[k] * c1[k],
                     factor * c1[k]**2 - 1))
    guards = [*roots, *W,
              *[W[i] - W[j] for i, j in combinations(range(3), 2)],
              *[cubes[i] - cubes[j] for i, j in combinations(range(3), 2)],
              determinant]
    return determinant, ties, guards, W


def quad_resultant(f, g):
    c, b, a = f
    F, e, d = g
    return (a * F - c * d)**2 - (a * e - b * d) * (b * F - c * e)


def multidegree(f):
    terms = f.to_dict()
    return [max((e[0] for e in terms), default=-1),
            max((e[1] for e in terms), default=-1)]


def main():
    _, ties, guards, _ = derive((1, 2, 4, 1))
    print("derived original-row tie coefficients", flush=True)
    resultants = [quad_resultant(ties[i], ties[j])
                  for i, j in combinations(range(3), 2)]
    polys = [r.n for r in resultants]
    degrees = [multidegree(f) for f in polys]
    gcd = polys[0].gcd(polys[1]).gcd(polys[2])
    if degrees != [[102, 102], [96, 102], [102, 96]] or gcd.total_degree() != 0:
        raise AssertionError({"degrees": degrees, "gcd": str(gcd)})

    # One x-slice checks positive y-degree and coprimality; y-slices check the
    # possible pure-x common factor. Full bivariate data is retained below.
    slice_x = [f.subs({"x": 6}) for f in polys[:2]]
    slice_x_degrees = [int(f.degrees()[1]) for f in slice_x]
    slice_x_gcd = slice_x[0].gcd(slice_x[1])
    slice_y_first = [polys[0].subs({"y": y0}) for y0 in (2, 6)]
    slice_y_first_degrees = [int(f.degrees()[0]) for f in slice_y_first]
    slice_y_first_gcd = slice_y_first[0].gcd(slice_y_first[1])
    if slice_x_degrees != [102, 102] or slice_x_gcd.degrees()[1] != 0:
        raise AssertionError((slice_x_degrees, str(slice_x_gcd)))
    if slice_y_first_degrees != [102, 102] or slice_y_first_gcd.degrees()[0] != 0:
        raise AssertionError((slice_y_first_degrees, str(slice_y_first_gcd)))

    data = {
        "diagram": 94,
        "field": "F_17[x,y] via python-flint nmod_mpoly",
        "mass_witness_mod17": [1, 2, 4, 1],
        "model": "raw pair resultants from original-row wing-cube solve and signed ties",
        "pair_order": [[0, 1], [0, 2], [1, 2]],
        "resultant_multidegrees_xy": degrees,
        "bivariate_gcd_degree": int(gcd.total_degree()),
        "x6_slice_y_degrees": slice_x_degrees,
        "x6_slice_gcd_y_degree": int(slice_x_gcd.degrees()[1]),
        "first_resultant_y_slices": {
            "y_values": [2, 6],
            "x_degrees": slice_y_first_degrees,
            "gcd_degree": int(slice_y_first_gcd.degrees()[0])
        },
        "guard_count": len(guards),
        "resultant_sparse_terms": [
            [[int(e[0]), int(e[1]), int(c) % P] for e, c in sorted(f.to_dict().items())]
            for f in polys
        ],
        "scope": "This certifies the raw finite-field 102/96-degree model, not the manuscript's quotient-degree-120 model or its norm/saturation claims."
    }
    out = Path(__file__).with_name("no94_raw_bivariate_fiber_certificate.json")
    out.write_text(json.dumps(data, indent=2, default=int) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in data.items() if k != "resultant_sparse_terms"},
                     indent=2, default=int))


if __name__ == "__main__":
    main()
