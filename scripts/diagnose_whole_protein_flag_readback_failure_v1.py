#!/usr/bin/env python3
"""Diagnose one saved finite-difference discrepancy without changing any fit."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from cached_matched_ml import CachedMatchedML, profiled_ml
from matched_ml_gradient import evaluate_gradient
from matched_mixed_covariance import MatchedCovariance
from whole_protein_flag_followup import arrays


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--fit-plan", type=Path, required=True)
    parser.add_argument("--input-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    saved = json.loads(args.case.read_text())
    identity = saved["identity"]
    fit_plan = json.loads(args.fit_plan.read_text())
    recipe = None
    with args.input_manifest.open() as handle:
        for line in handle:
            candidate = json.loads(line)
            if candidate["fit_input_id"] == identity["fit_input_id"]:
                recipe = candidate
                break
    if recipe is None or recipe["sha256"] != identity["input_sha256"]:
        raise ValueError("case input is absent or differs from the materialized-input manifest")
    original_item = {
        "fit_input_id": identity["fit_input_id"], "tree": identity["tree"],
        "path": identity["original_fit"], "sha256": identity["original_fit_sha256"],
    }
    original, _, background, family, factor, x, y, scales = arrays(original_item, recipe, fit_plan)
    result = saved["result"]
    theta = np.asarray(result["selected_theta"], dtype=float)
    if theta.shape != (3,) or not np.isfinite(theta).all():
        raise ValueError("saved selected parameters are invalid")
    saved_gradient = np.asarray(result["refinement"]["analytic_gradient"], dtype=float)
    recalculated_gradient = evaluate_gradient(CachedMatchedML(background, family, factor, x, y), np.expm1(theta))["log1p_ratio_gradient"]
    np.testing.assert_allclose(recalculated_gradient, saved_gradient, rtol=1e-9, atol=1e-9)
    upper = float(np.log1p(result["refinement"]["maximum_ratio"]))

    def objective(point: np.ndarray) -> float:
        return float(profiled_ml(MatchedCovariance(background, family, factor, 1.0, *np.expm1(point)), x, y)["negative_profiled_ml"])

    diagnostics = []
    for step in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
        values = []
        schemes = []
        for axis in range(3):
            direction = np.eye(3)[axis] * step
            if theta[axis] < step:
                estimate = (-3 * objective(theta) + 4 * objective(theta + direction) - objective(theta + 2 * direction)) / (2 * step)
                scheme = "forward_three_point"
            elif theta[axis] > upper - step:
                estimate = (3 * objective(theta) - 4 * objective(theta - direction) + objective(theta - 2 * direction)) / (2 * step)
                scheme = "backward_three_point"
            else:
                estimate = (objective(theta + direction) - objective(theta - direction)) / (2 * step)
                scheme = "central_two_point"
            values.append(float(estimate)); schemes.append(scheme)
        values = np.asarray(values)
        diagnostics.append({"step": step, "scheme_by_axis": schemes, "finite_difference_gradient": values.tolist(),
                            "difference_from_saved_analytic_gradient": (values - saved_gradient).tolist(),
                            "maximum_absolute_difference": float(np.max(np.abs(values - saved_gradient)))})
    payload = {
        "schema_version": 1,
        "status": "complete_single_case_finite_difference_diagnostic",
        "case": str(args.case), "case_sha256": sha256(args.case),
        "fit_plan": str(args.fit_plan), "fit_plan_sha256": sha256(args.fit_plan),
        "input_manifest": str(args.input_manifest), "input_manifest_sha256": sha256(args.input_manifest),
        "identity": identity, "selected_theta": theta.tolist(), "upper_log1p_ratio": upper,
        "saved_analytic_gradient": saved_gradient.tolist(), "step_diagnostics": diagnostics,
        "scope": "One saved numerical-review case only. It records finite-difference sensitivity without replacing the saved fit, accepting the follow-up result, changing a tolerance, or making any sequence-structure inference."
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": payload["status"], "maximum_absolute_differences": [x["maximum_absolute_difference"] for x in diagnostics]}))


if __name__ == "__main__":
    main()
