"""Independent rational interval replay for the 752 No. 93 torque labels.

The root routines reconstruct the signed Euler quintics and Lagrange phase
quadratics from the displayed three-body equations. Candidate centers are
numeric, but every accepted root count and nonzero torque bound uses exact
Fraction arithmetic and a Rouché inequality.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction as F
from pathlib import Path

from generate_no98_interval_certificate import (
    cscale, cross_root_disk, cinf, disk, euler_roots,
    lagrange_root_disks, positions_euler,
)

PUBLISHED_MARGIN = 3050302563971362176000000
MARGIN_SCALE = 10**24


def main() -> None:
    masses = (1, 2, 4)
    euler = euler_roots(masses)
    lagrange = lagrange_root_disks(masses)
    shapes = [("Euler", row["signs"], positions_euler(row, masses)) for row in euler]
    shapes += [("Lagrange", row["phase_ids"], row["W_disks"]) for row in lagrange]
    assert len(shapes) == 94
    one = disk((F(1), F(0)))
    least = None
    least_label = None
    count = 0
    for kind, label, wing in shapes:
        roots = [cross_root_disk(one, w) for w in wing]
        centers = [cscale(root[0], mass) for root, mass in zip(roots, masses)]
        radius = sum((mass * root[1] for root, mass in zip(roots, masses)), F(0))
        for signs in itertools.product((-1, 1), repeat=3):
            real = sum((sign * center[0] for sign, center in zip(signs, centers)), F(0))
            imag = sum((sign * center[1] for sign, center in zip(signs, centers)), F(0))
            lower = max(abs(real), abs(imag)) - radius
            if lower <= 0:
                raise ArithmeticError(("No. 93 torque disk contains zero", kind, label, signs, lower))
            if least is None or lower < least:
                least, least_label = lower, (kind, label, signs)
            count += 1
    assert count == (40 + 54) * 8 == 752
    computed_floor = int(least * MARGIN_SCALE)
    assert computed_floor == PUBLISHED_MARGIN
    result = {
        "diagram": 93,
        "masses_m1_m2_m3": list(masses),
        "Euler_root_labels": len(euler),
        "Lagrange_root_labels": len(lagrange),
        "signed_torque_labels": count,
        "all_torque_disks_exclude_zero": True,
        "least_exact_lower_bound": str(least),
        "margin_scale": MARGIN_SCALE,
        "least_lower_bound_times_1e24_floor": computed_floor,
        "manuscript_printed_margin": PUBLISHED_MARGIN,
        "manuscript_margin_is_conservative_lower_bound_at_same_scale": True,
        "manuscript_margin_matches_computed_least": computed_floor == PUBLISHED_MARGIN,
        "least_label": least_label,
        "scope": "All root labels and torque intervals at one positive mass witness; original diagram-to-torque transition is a separate analytic step.",
    }
    out = Path(__file__).with_name("no93_interval_certificate.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in (
        "Euler_root_labels", "Lagrange_root_labels", "signed_torque_labels",
        "least_lower_bound_times_1e24_floor", "manuscript_printed_margin",
        "manuscript_margin_is_conservative_lower_bound_at_same_scale",
        "manuscript_margin_matches_computed_least")}, indent=2))


if __name__ == "__main__":
    main()
