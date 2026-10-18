"""Collect the latest curated appendix replay records without certifying gaps."""
from __future__ import annotations
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
MANIFEST = HERE / "appendix_replay_manifest.json"

def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    jobs = [HERE.parent / "run_logs" / "verify_n6_structure.run.json"]
    names = [
        "generate_no8_revised_order_certificate", "generate_shared_energy_certificate",
        "verify_no20_energy",
        "generate_no50_resonance_certificate", "generate_no79_bridge_certificate",
        "verify_no79_shape_energy",
        "generate_no79_sturm_certificate", "verify_no79_hurwitz",
        "generate_no52_certificate", "verify_no52_valuation_partition",
        "verify_revised_core_claims", "verify_no66_function_field_cover",
        "generate_no88_collinear_certificate", "generate_no88_fiber_certificate",
        "generate_no88_vertical_certificate", "generate_no88_smooth_certificate",
        "generate_no88_zero_star_certificate", "verify_no88_root_norm_identity",
        "verify_no88_CM_edge_quadratic", "generate_crossed_torque_certificate",
        "generate_no93_interval_certificate", "generate_no98_interval_certificate",
        "generate_no63_certificate", "derive_no94_finite_fiber",
        "generate_no94_structure_certificate", "generate_no94_equal_spoke_certificate",
        "generate_no94_projective_spoke_cover", "generate_no94_local_certificates",
        "verify_no94_distinct_spoke_fiber", "verify_no94_h3_fiber",
    ]
    jobs += [HERE / "run_logs" / f"{name}.run.json" for name in names]
    records = []
    for path in jobs:
        if not path.exists():
            records.append({"record": path.name, "status": "missing"})
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        records.append({"record": path.name, "script": data.get("script"),
                        "exit_code": data.get("exit_code"),
                        "elapsed_seconds": data.get("elapsed_seconds"),
                        "peak_rss_mb": data.get("peak_rss_mb"),
                        "status": "passed" if data.get("exit_code") == 0 else "failed"})
    cert93 = json.loads((HERE / "no93_interval_certificate.json").read_text(encoding="utf-8"))
    cert98 = json.loads((HERE / "no98_interval_certificate.json").read_text(encoding="utf-8"))
    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "main_merged_revised.tex",
        "jobs": records,
        "all_listed_jobs_passed": all(r.get("status") == "passed" for r in records),
        "no93": {k: cert93[k] for k in ("Euler_root_labels", "Lagrange_root_labels", "signed_torque_labels", "least_exact_lower_bound", "least_lower_bound_times_1e24_floor", "manuscript_printed_margin", "manuscript_margin_is_conservative_lower_bound_at_same_scale", "manuscript_margin_matches_computed_least")},
        "no98": {"masses_upper": cert98["masses_upper"], "masses_lower": cert98["masses_lower"], "upper_euler_root_disks": cert98["upper_euler_root_disks"], "lower_euler_root_disks": cert98["lower_euler_root_disks"], "lower_lagrange_root_disks": cert98["lower_lagrange_root_disks"], "margin_scale": cert98["margin_scale"], "manuscript_printed_margin": cert98["manuscript_printed_margin"], "manuscript_margin_normalization_defined_in_source": cert98["manuscript_margin_normalization_defined_in_source"], **cert98["torque_sweep"]},
        "appendix_coverage": manifest["entries"],
        "known_claim_mismatches": manifest["known_claim_mismatches"],
        "no63": json.loads((HERE / "no63_certificate.json").read_text(encoding="utf-8")),
        "no94_local_certificates": json.loads((HERE / "no94_local_certificates.json").read_text(encoding="utf-8")),
        "no66_function_field": json.loads((HERE / "no66_function_field_certificate.json").read_text(encoding="utf-8")),
        "no52_valuation_partition": json.loads((HERE / "no52_valuation_partition_certificate.json").read_text(encoding="utf-8")),
        "no94_finite_fiber_slice": json.loads((HERE / "no94_finite_fiber_slice_certificate.json").read_text(encoding="utf-8")),
        "no94_raw_bivariate_fiber": json.loads((HERE / "no94_raw_bivariate_fiber_certificate.json").read_text(encoding="utf-8")),
        "no94_h3_raw_fiber": json.loads((HERE / "no94_h3_raw_fiber_certificate.json").read_text(encoding="utf-8")),
        "note": "A passing computational job certifies only its recorded scope; blocked or analytic proof obligations remain outside the computational replay.",
    }
    out = HERE / "appendix_replay_results.json"
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Saved {out}; {sum(r.get('status') == 'passed' for r in records)}/{len(records)} jobs passed.")
    if not result["all_listed_jobs_passed"]:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
