"""Small, guarded residual surrogate for future MusselFlow calibration.

This module intentionally does not invent training data and is not enabled by
the Grasshopper optimizer.  It provides the safe numerical layer needed once
paired baseline/CFD/flume/field observations exist:

    residual = trusted_observation - reduced_order_baseline

An ensemble of ridge-regressed quadratic models estimates the residual and its
between-model uncertainty.  Every prediction also receives an out-of-domain
test.  The ecological grammar and hard constraints remain outside the learned
model, so a new policy limit does not require retraining and cannot be silently
overridden by the surrogate.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from musselflow_ecogrammar_core import resolved_lists


MODEL_SCHEMA_VERSION = 1


LAYOUT_FEATURE_NAMES = [
    "log_obstacle_count",
    "log_domain_area_m2",
    "plan_coverage_fraction",
    "major_mean_m",
    "major_cv",
    "minor_mean_m",
    "minor_cv",
    "aspect_mean",
    "height_mean_m",
    "height_cv",
    "porosity_mean",
    "porosity_std",
    "vertical_blockage_mean",
    "nearest_spacing_over_major",
    "current_speed_mean_m_s",
    "current_speed_max_m_s",
    "current_speed_cv",
    "weighted_normal_orientation",
    "log1p_twine_reynolds_mean",
    "temperature_c",
    "salinity_psu",
    "log1p_chlorophyll_ug_l",
    "log1p_tsm_mg_l",
    "log1p_total_animals",
    "log1p_dry_biomass_kg",
    "log1p_baseline_capture_g_day",
    "baseline_probe_speed_ratio_mean",
    "baseline_obstacle_food_mean",
    "baseline_min_do_mg_l",
    "log1p_baseline_deposition_kg_m2_day",
]


def _as_2d(values, name):
    array = np.asarray(values, dtype=float)
    if array.ndim == 1:
        array = array.reshape((1, -1))
    if array.ndim != 2 or not np.all(np.isfinite(array)):
        raise ValueError("%s must be a finite 2D numeric array." % name)
    return array


def quadratic_basis(standardised_features):
    """Linear plus upper-triangular quadratic terms, without an intercept."""
    z = _as_2d(standardised_features, "standardised_features")
    feature_count = z.shape[1]
    first, second = np.triu_indices(feature_count)
    quadratic = z[:, first]*z[:, second]
    return np.concatenate((z, quadratic), axis=1)


def _coefficient_of_variation(values):
    values = np.asarray(values, dtype=float)
    mean = float(np.mean(values))
    return float(np.std(values)/mean) if mean > 1e-12 else 0.0


def extract_layout_features(obstacles, domain_polygon, flow_vectors,
                            config, baseline_result):
    """Build one explicit residual-learning feature row.

    The feature schema is versioned by ``LAYOUT_FEATURE_NAMES``.  A future
    model package must be retrained if this schema or represented physics
    changes.
    """
    obstacles = _as_2d(obstacles, "obstacles")
    if obstacles.shape[1] != 6 or len(obstacles) == 0:
        raise ValueError("obstacles must have six descriptor columns.")
    polygon = _as_2d(domain_polygon, "domain_polygon")
    if polygon.shape[1] != 2 or len(polygon) < 3:
        raise ValueError("domain_polygon must have at least three XY rows.")
    flows = _as_2d(flow_vectors, "flow_vectors")
    if flows.shape[1] != 2 or len(flows) == 0:
        raise ValueError("flow_vectors must have two plan components.")
    config = resolved_lists(config, len(obstacles), len(flows))
    weights = np.asarray(config["scenario.weights"], dtype=float)

    x = polygon[:, 0]
    y = polygon[:, 1]
    area = 0.5*abs(float(
        np.dot(x, np.roll(y, -1))-np.dot(y, np.roll(x, -1))))
    if area <= 0.0:
        raise ValueError("domain_polygon has zero area.")

    major = obstacles[:, 2]
    minor = obstacles[:, 3]
    height = obstacles[:, 5]
    porosity = np.asarray(config["structure.porosity"], dtype=float)
    centres = obstacles[:, :2]
    if len(centres) > 1:
        delta = centres[:, None, :]-centres[None, :, :]
        distance = np.sqrt(np.sum(delta*delta, axis=2))
        np.fill_diagonal(distance, np.inf)
        nearest_spacing = float(np.mean(np.min(distance, axis=1)))
    else:
        nearest_spacing = math.sqrt(area)

    speeds = np.linalg.norm(flows, axis=1)
    directions = np.divide(
        flows, speeds[:, None],
        out=np.tile(np.array([[1.0, 0.0]]), (len(flows), 1)),
        where=speeds[:, None] > 1e-15)
    major_axes = np.column_stack((
        np.cos(obstacles[:, 4]), np.sin(obstacles[:, 4])))
    # 0 = aligned with flow, 1 = normal/bluff to flow.
    normal_orientation = []
    for direction in directions:
        normal_orientation.append(float(np.mean(
            1.0-np.abs(major_axes@direction))))
    weighted_orientation = float(np.dot(weights, normal_orientation))
    reynolds = (
        speeds*config["structure.twine_diameter_m"] /
        config["hydrodynamics.kinematic_viscosity_m2_s"])

    animal_count = np.asarray(
        config["stocking.mussels_per_obstacle"], dtype=float)
    dry_override = np.asarray(
        config["stocking.dry_tissue_kg_per_obstacle"], dtype=float)
    inferred_dry = (
        animal_count*config["species.mean_dry_tissue_g"]/1000.0)
    if config["stocking.mode"].lower() == "dry_biomass":
        dry_biomass = dry_override
        animal_count = (
            dry_biomass*1000.0 /
            config["species.mean_dry_tissue_g"])
    else:
        dry_biomass = np.where(
            dry_override >= 0.0, dry_override, inferred_dry)
        animal_count = np.where(
            dry_override >= 0.0,
            dry_biomass*1000.0 /
            config["species.mean_dry_tissue_g"],
            animal_count)
    deposition = baseline_result["oxygen"][
        "weighted_deposition_kg_m2_day"]

    row = np.asarray([
        math.log(max(len(obstacles), 1)),
        math.log(max(area, 1e-12)),
        float(np.sum(major*minor))/area,
        float(np.mean(major)),
        _coefficient_of_variation(major),
        float(np.mean(minor)),
        _coefficient_of_variation(minor),
        float(np.mean(major/np.maximum(minor, 1e-12))),
        float(np.mean(height)),
        _coefficient_of_variation(height),
        float(np.mean(porosity)),
        float(np.std(porosity)),
        float(np.mean(np.clip(
            height/config["site.depth_m"], 0.0, 1.0))),
        nearest_spacing/max(float(np.mean(major)), 1e-12),
        float(np.dot(weights, speeds)),
        float(np.max(speeds)),
        _coefficient_of_variation(speeds),
        weighted_orientation,
        math.log1p(float(np.dot(weights, reynolds))),
        config["site.temperature_c"],
        config["site.salinity_psu"],
        math.log1p(config["site.chlorophyll_ug_l"]),
        math.log1p(config["site.tsm_mg_l"]),
        math.log1p(float(np.sum(animal_count))),
        math.log1p(float(np.sum(dry_biomass))),
        math.log1p(baseline_result["chlorophyll_capture_g_day"]),
        float(np.mean(baseline_result["probe_speed_ratio"])),
        float(np.mean(baseline_result["obstacle_food_fraction"])),
        baseline_result["oxygen"]["minimum_mg_l"],
        math.log1p(max(float(deposition), 0.0)),
    ], dtype=float)
    if len(row) != len(LAYOUT_FEATURE_NAMES):
        raise AssertionError("Internal layout feature schema mismatch.")
    return row, list(LAYOUT_FEATURE_NAMES)


def _ridge_fit(basis, targets, regularisation):
    identity = np.eye(basis.shape[1])
    matrix = basis.T@basis+regularisation*identity
    right_hand_side = basis.T@targets
    try:
        return np.linalg.solve(matrix, right_hand_side)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(matrix, right_hand_side, rcond=None)[0]


def fit_residual_ensemble(features, residual_targets, feature_names,
                          target_names, ensemble_size=12,
                          regularisation=1.0e-3, seed=20260728):
    """Fit deterministic bootstrap ridge models and return a model dictionary."""
    x = _as_2d(features, "features")
    y = _as_2d(residual_targets, "residual_targets")
    if len(x) != len(y):
        raise ValueError("features and residual_targets need equal row counts.")
    if len(x) < max(12, x.shape[1]+2):
        raise ValueError(
            "Too few training rows for even this screening surrogate.")
    feature_names = [str(value) for value in feature_names]
    target_names = [str(value) for value in target_names]
    if len(feature_names) != x.shape[1]:
        raise ValueError("feature_names length does not match features.")
    if len(target_names) != y.shape[1]:
        raise ValueError("target_names length does not match targets.")
    if len(set(feature_names)) != len(feature_names):
        raise ValueError("feature_names must be unique.")
    if len(set(target_names)) != len(target_names):
        raise ValueError("target_names must be unique.")

    x_mean = np.mean(x, axis=0)
    x_scale = np.std(x, axis=0)
    x_scale = np.where(x_scale > 1.0e-12, x_scale, 1.0)
    z = (x-x_mean)/x_scale
    basis = quadratic_basis(z)
    y_mean = np.mean(y, axis=0)
    centred_y = y-y_mean

    generator = np.random.default_rng(int(seed))
    coefficients = []
    row_count = len(x)
    for ensemble_index in range(max(2, int(ensemble_size))):
        indices = generator.integers(0, row_count, row_count)
        coefficients.append(_ridge_fit(
            basis[indices], centred_y[indices], float(regularisation)))
    coefficients = np.asarray(coefficients, dtype=float)

    train_members = np.einsum(
        "nb,ebt->ent", basis, coefficients)+y_mean[None, None, :]
    train_prediction = np.mean(train_members, axis=0)
    errors = train_prediction-y
    rmse = np.sqrt(np.mean(errors*errors, axis=0))
    mae = np.mean(np.abs(errors), axis=0)

    return {
        "schema_version": MODEL_SCHEMA_VERSION,
        "feature_names": feature_names,
        "target_names": target_names,
        "x_mean": x_mean,
        "x_scale": x_scale,
        "x_min": np.min(x, axis=0),
        "x_max": np.max(x, axis=0),
        "y_mean": y_mean,
        "coefficients": coefficients,
        "train_rmse": rmse,
        "train_mae": mae,
        "training_rows": row_count,
        "regularisation": float(regularisation),
        "seed": int(seed),
        "basis": "standardised_linear_plus_upper_quadratic",
    }


def predict_residual(features, model, maximum_standard_score=3.0,
                     range_padding_fraction=0.05):
    """Return mean residual, uncertainty, and per-row in-domain diagnostics."""
    x = _as_2d(features, "features")
    if int(model["schema_version"]) != MODEL_SCHEMA_VERSION:
        raise ValueError("Unsupported surrogate model schema.")
    if x.shape[1] != len(model["feature_names"]):
        raise ValueError("Feature count does not match surrogate package.")
    x_mean = np.asarray(model["x_mean"], dtype=float)
    x_scale = np.asarray(model["x_scale"], dtype=float)
    z = (x-x_mean)/x_scale
    basis = quadratic_basis(z)
    coefficients = np.asarray(model["coefficients"], dtype=float)
    y_mean = np.asarray(model["y_mean"], dtype=float)
    members = np.einsum(
        "nb,ebt->ent", basis, coefficients)+y_mean[None, None, :]
    mean = np.mean(members, axis=0)
    standard_deviation = np.std(members, axis=0, ddof=1)

    x_min = np.asarray(model["x_min"], dtype=float)
    x_max = np.asarray(model["x_max"], dtype=float)
    padding = np.maximum(
        (x_max-x_min)*float(range_padding_fraction), 1.0e-12)
    inside_range = np.all(
        (x >= x_min-padding) & (x <= x_max+padding), axis=1)
    maximum_z = np.max(np.abs(z), axis=1)
    in_domain = inside_range & (
        maximum_z <= float(maximum_standard_score))
    return {
        "mean": mean,
        "standard_deviation": standard_deviation,
        "in_domain": in_domain,
        "maximum_standard_score": maximum_z,
        "member_predictions": members,
    }


def save_model(path, model):
    """Write a portable NPZ package; intended for offline training only."""
    path = Path(path)
    metadata = {
        "schema_version": int(model["schema_version"]),
        "feature_names": list(model["feature_names"]),
        "target_names": list(model["target_names"]),
        "training_rows": int(model["training_rows"]),
        "regularisation": float(model["regularisation"]),
        "seed": int(model["seed"]),
        "basis": str(model["basis"]),
    }
    np.savez_compressed(
        path,
        metadata=np.asarray(json.dumps(metadata)),
        x_mean=np.asarray(model["x_mean"], dtype=float),
        x_scale=np.asarray(model["x_scale"], dtype=float),
        x_min=np.asarray(model["x_min"], dtype=float),
        x_max=np.asarray(model["x_max"], dtype=float),
        y_mean=np.asarray(model["y_mean"], dtype=float),
        coefficients=np.asarray(model["coefficients"], dtype=float),
        train_rmse=np.asarray(model["train_rmse"], dtype=float),
        train_mae=np.asarray(model["train_mae"], dtype=float),
    )


def load_model(path):
    """Load and validate a residual-surrogate NPZ package."""
    with np.load(Path(path), allow_pickle=False) as package:
        metadata = json.loads(str(package["metadata"].item()))
        required = {
            "x_mean", "x_scale", "x_min", "x_max", "y_mean",
            "coefficients", "train_rmse", "train_mae",
        }
        missing = required-set(package.files)
        if missing:
            raise ValueError(
                "Surrogate package is missing arrays: %s"
                % ", ".join(sorted(missing)))
        model = dict(metadata)
        for key in required:
            model[key] = np.asarray(package[key], dtype=float)
    if int(model.get("schema_version", -1)) != MODEL_SCHEMA_VERSION:
        raise ValueError("Unsupported surrogate model schema.")
    if not np.all(np.isfinite(model["coefficients"])):
        raise ValueError("Surrogate coefficients contain non-finite values.")
    return model
