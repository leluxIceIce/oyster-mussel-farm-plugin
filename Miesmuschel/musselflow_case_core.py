"""MusselFlow ecological-case parsing and solver-adapter utilities.

This module deliberately has no Rhino or Grasshopper dependency:

* one strict ecological grammar JSON;
* explicit timeline versus ensemble semantics;
* multiple flow-vector references;
* deterministic canonical serialization and hashing;
* compilation into the flat configuration consumed by the numerical core.

It does not validate the scientific model and cannot grant calibrated status.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path


SCHEMA_VERSION = "1.0.0"
COMPATIBLE_SCHEMA_VERSIONS = {"2.0.0-draft"}
MODEL_ROLE = "UNVALIDATED_LAYOUT_SCREENING"
RUNTIME_STATUS = "UNVALIDATED_SCREENING"


class CaseValidationError(ValueError):
    """Raised when one or more structural or numeric case errors exist."""

    def __init__(self, errors):
        self.errors = list(errors)
        super().__init__(" | ".join(self.errors))


class UnsupportedCaseError(ValueError):
    """Raised when valid case semantics are not implemented by an evaluator."""


TOP_KEYS = {
    "schema_version", "case_id", "model_role", "units", "forcing", "site",
    "species", "stocking", "structure", "hydrodynamics", "sediment",
    "harvest", "constraints", "objective", "validation", "machine_learning",
    "food", "evidence",
}
REQUIRED_TOP_KEYS = TOP_KEYS-{"machine_learning", "food"}
UNIT_KEYS = {
    "geometry", "flow_vector_length", "duration", "temperature", "salinity",
    "dissolved_oxygen", "chlorophyll_a", "suspended_particulate_matter",
}
BOUNDARY_KEYS = {
    "temperature_c", "salinity_psu", "dissolved_oxygen_mg_l",
    "chlorophyll_a_ug_l", "tsm_mg_l",
}
SITE_KEYS = {
    "depth_m", "initial_dissolved_oxygen_mg_l",
    "particulate_organic_fraction",
    "background_sediment_oxygen_demand_g_o2_m2_day",
    "pelagic_respiration_g_o2_m3_day",
    "primary_production_g_o2_m3_day", "reaeration_per_day",
    "vertical_exchange_per_day", "advective_exchange_efficiency",
}
SPECIES_REQUIRED_KEYS = {
    "taxon", "profile_status", "mean_dry_tissue_g", "size_cv", "clearance",
    "respiration", "retention_efficiency",
    "particulate_retention_efficiency", "assimilation_efficiency",
    "pseudofaeces_fraction", "activity_fraction", "salinity_activity_psu",
    "ammonia_mg_n_g_dw_h", "low_food_threshold_ug_l",
    "low_food_transition_ug_l", "valid_flow_max_m_s",
}
SPECIES_KEYS = SPECIES_REQUIRED_KEYS | {
    "retention_d50_um", "retention_slope_um",
    "ingestion_half_saturation_ug_l",
    "high_food_assimilation_threshold_ug_l",
    "high_food_assimilation_decay_ug_l",
    "assimilation_quality_half_saturation",
    "assimilation_reference_organic_fraction",
    "pseudofaeces_tsm_threshold_mg_l",
    "pseudofaeces_tsm_transition_mg_l", "pseudofaeces_max_fraction",
    "oxygen_zero_saturation_fraction", "oxygen_full_saturation_fraction",
    "current_clearance_start_m_s", "current_clearance_zero_m_s",
    "current_protection_group_size", "oxycalorific_kj_per_ml_o2",
    "tissue_energy_kj_g_dw",
}
ALLOMETRY_COMMON_KEYS = {
    "mass_exponent", "reference_temperature_c", "q10",
}
SALINITY_KEYS = {"zero_low", "full_low", "full_high", "zero_high"}
STOCKING_REQUIRED_KEYS = {
    "mode", "mussels_per_obstacle", "live_wet_g_per_individual",
    "annual_mortality_fraction",
}
STOCKING_KEYS = STOCKING_REQUIRED_KEYS | {"effective_aggregation_size"}
FOOD_KEYS = {
    "carbon_to_chlorophyll_mg_c_per_mg_chl", "organic_carbon_fraction",
    "small_particle_fraction", "small_particle_diameter_um",
    "large_particle_diameter_um", "detritus_fraction_of_organic",
    "detritus_preference", "detritus_assimilation_multiplier",
    "organic_energy_kj_g",
}
STRUCTURE_KEYS = {
    "porosity_per_obstacle", "twine_diameter_m", "drag_multiplier",
    "fallback_plan_size_m", "fallback_height_m",
}
SEDIMENT_KEYS = {
    "settling_velocity_m_s", "in_domain_deposition_fraction",
    "oxygen_demand_kg_o2_per_kg_organic", "decay_per_day",
    "resuspension_per_day", "mortality_deposition_fraction",
    "initial_organic_stock_kg",
}
HARVEST_KEYS = {
    "fraction_per_year", "turnovers_per_year", "n_kg_per_t_wet",
    "p_kg_per_t_wet",
}
HYDRO_KEYS = {
    "wake_spread", "food_plume_spread", "food_recovery_lengths",
    "minimum_speed_ratio", "kinematic_viscosity_m2_s",
    "water_density_kg_m3",
}
CONSTRAINT_KEYS = {
    "minimum_dissolved_oxygen_mg_l", "minimum_probe_speed_m_s",
    "maximum_biodeposition_kg_m2_day", "minimum_obstacle_clearance_m",
}
OBJECTIVE_KEYS = {
    "weights", "target_chlorophyll_capture_g_day",
    "target_biodeposition_kg_m2_day",
}
WEIGHT_KEYS = {
    "extraction", "flushing", "food_delivery", "uniformity", "oxygen",
    "low_deposition",
}
VALIDATION_KEYS = {
    "status", "evidence_manifest_id", "allow_extrapolation",
}
ML_KEYS = {
    "surrogate_enabled", "model_artifact_id", "training_dataset_id",
    "reinforcement_learning",
}
RL_KEYS = {"enabled", "intended_future_role", "static_layout_optimizer"}
EVIDENCE_KEYS = {
    "parameter", "citation", "equation", "validity_envelope", "doi"}
LEGACY_EVIDENCE_KEYS = {"parameter", "doi"}


def _load_source(source):
    if isinstance(source, dict):
        return copy.deepcopy(source)
    if isinstance(source, Path):
        text = source.read_text(encoding="utf-8")
    elif isinstance(source, str):
        text = source
    else:
        raise TypeError("case source must be a dict, JSON string, or pathlib.Path")
    try:
        value = json.loads(text, parse_constant=lambda value: (_ for _ in ()).throw(
            ValueError("non-finite JSON constant %s" % value)))
    except (json.JSONDecodeError, ValueError) as exception:
        raise CaseValidationError(["Invalid JSON: %s" % exception]) from exception
    if not isinstance(value, dict):
        raise CaseValidationError(["Case root must be a JSON object."])
    return value


def _shape(mapping, required, allowed, path, errors):
    if not isinstance(mapping, dict):
        errors.append("%s must be an object." % path)
        return False
    missing = sorted(set(required)-set(mapping))
    unknown = sorted(set(mapping)-set(allowed))
    if missing:
        errors.append("%s is missing: %s." % (path, ", ".join(missing)))
    if unknown:
        errors.append("%s has unknown keys: %s." % (path, ", ".join(unknown)))
    return not missing and not unknown


def _number(value, path, errors, lower=None, upper=None, strict_lower=False,
            nullable=False):
    if value is None and nullable:
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        errors.append("%s must be a number%s." % (
            path, " or null" if nullable else ""))
        return
    if not math.isfinite(float(value)):
        errors.append("%s must be finite." % path)
        return
    if lower is not None:
        invalid = value <= lower if strict_lower else value < lower
        if invalid:
            errors.append("%s must be %s %s." % (
                path, ">" if strict_lower else ">=", lower))
    if upper is not None and value > upper:
        errors.append("%s must be <= %s." % (path, upper))


def _integer(value, path, errors, lower=0):
    if isinstance(value, bool) or not isinstance(value, int):
        errors.append("%s must be an integer." % path)
    elif value < lower:
        errors.append("%s must be >= %d." % (path, lower))


def _text(value, path, errors):
    if not isinstance(value, str) or not value.strip():
        errors.append("%s must be non-empty text." % path)


def _scalar_or_array(value, path, errors, lower=0.0, upper=None):
    values = value if isinstance(value, list) else [value]
    if not values:
        errors.append("%s cannot be empty." % path)
        return
    for index, item in enumerate(values):
        _number(item, "%s[%d]" % (path, index), errors,
                lower=lower, upper=upper)


def _boundary(value, path, errors):
    if not _shape(value, BOUNDARY_KEYS, BOUNDARY_KEYS, path, errors):
        return
    _number(value["temperature_c"], path+".temperature_c", errors, -2, 40)
    _number(value["salinity_psu"], path+".salinity_psu", errors, 0, 42)
    for key in (
            "dissolved_oxygen_mg_l", "chlorophyll_a_ug_l", "tsm_mg_l"):
        _number(value[key], path+"."+key, errors, 0)


def _validate_forcing(forcing, flow_count, errors, warnings):
    if not isinstance(forcing, dict):
        errors.append("forcing must be an object.")
        return
    mode = forcing.get("mode")
    if mode == "timeline":
        allowed = {"mode", "vector_source", "steps", "repeat_count"}
        if not _shape(forcing, allowed, allowed, "forcing", errors):
            return
        if forcing["vector_source"] != "Grasshopper.flowVectors":
            errors.append("forcing.vector_source must be Grasshopper.flowVectors.")
        _integer(forcing["repeat_count"], "forcing.repeat_count", errors, 1)
        items = forcing["steps"]
        item_key = "steps"
        item_allowed = {
            "id", "order", "flow_vector_index", "duration_h", "boundary"}
    elif mode == "ensemble":
        allowed = {"mode", "vector_source", "states"}
        if not _shape(forcing, allowed, allowed, "forcing", errors):
            return
        if forcing["vector_source"] != "Grasshopper.flowVectors":
            errors.append("forcing.vector_source must be Grasshopper.flowVectors.")
        items = forcing["states"]
        item_key = "states"
        item_allowed = {
            "id", "flow_vector_index", "occurrence_probability",
            "duration_h", "boundary"}
    else:
        errors.append("forcing.mode must be timeline or ensemble.")
        return

    if not isinstance(items, list) or not items:
        errors.append("forcing.%s must be a non-empty array." % item_key)
        return

    identifiers = []
    orders = []
    probability_sum = 0.0
    for index, item in enumerate(items):
        path = "forcing.%s[%d]" % (item_key, index)
        if not _shape(item, item_allowed, item_allowed, path, errors):
            continue
        _text(item["id"], path+".id", errors)
        identifiers.append(item["id"])
        _integer(item["flow_vector_index"], path+".flow_vector_index", errors)
        if (flow_count is not None and
                isinstance(item["flow_vector_index"], int) and
                not isinstance(item["flow_vector_index"], bool) and
                item["flow_vector_index"] >= flow_count):
            errors.append(
                "%s.flow_vector_index=%d but flow_count=%d."
                % (path, item["flow_vector_index"], flow_count))
        _boundary(item["boundary"], path+".boundary", errors)
        if mode == "timeline":
            _integer(item["order"], path+".order", errors)
            orders.append(item["order"])
            _number(item["duration_h"], path+".duration_h", errors,
                    lower=0, strict_lower=True)
        else:
            _number(item["duration_h"], path+".duration_h", errors,
                    lower=0, strict_lower=True)
            _number(item["occurrence_probability"],
                    path+".occurrence_probability", errors, 0, 1)
            if (isinstance(item["occurrence_probability"], (int, float)) and
                    not isinstance(item["occurrence_probability"], bool) and
                    math.isfinite(float(item["occurrence_probability"]))):
                probability_sum += item["occurrence_probability"]

    duplicates = sorted({
        value for value in identifiers if identifiers.count(value) > 1})
    if duplicates:
        errors.append("Scenario IDs must be unique: %s." % ", ".join(duplicates))
    if mode == "timeline" and all(
            isinstance(value, int) and not isinstance(value, bool)
            for value in orders):
        expected = list(range(len(items)))
        if sorted(orders) != expected:
            errors.append(
                "Timeline order values must be exactly %s." % expected)
    if mode == "ensemble":
        if probability_sum <= 0:
            errors.append("Ensemble occurrence probabilities must have positive sum.")
        elif abs(probability_sum-1.0) > 1e-9:
            warnings.append(
                "Ensemble probabilities sum to %.12g and were normalized."
                % probability_sum)
            for item in items:
                if isinstance(item, dict) and "occurrence_probability" in item:
                    item["occurrence_probability"] /= probability_sum


def _validate_case(case, flow_count, obstacle_count):
    errors = []
    warnings = []
    if not _shape(case, REQUIRED_TOP_KEYS, TOP_KEYS, "case", errors):
        return errors, warnings
    schema_version = case["schema_version"]
    if schema_version != SCHEMA_VERSION:
        if schema_version in COMPATIBLE_SCHEMA_VERSIONS:
            warnings.append(
                "Legacy schema %s was accepted for migration; use %s."
                % (schema_version, SCHEMA_VERSION))
        else:
            errors.append("schema_version must be %s." % SCHEMA_VERSION)
    _text(case["case_id"], "case_id", errors)
    if case["model_role"] != MODEL_ROLE:
        errors.append("model_role must be %s." % MODEL_ROLE)

    units = case["units"]
    if _shape(units, UNIT_KEYS, UNIT_KEYS, "units", errors):
        expected_units = {
            "flow_vector_length": "m/s", "duration": "h",
            "temperature": "degC", "salinity": "PSU",
            "dissolved_oxygen": "mg/L", "chlorophyll_a": "ug/L",
            "suspended_particulate_matter": "mg/L",
        }
        _text(units["geometry"], "units.geometry", errors)
        for key, expected in expected_units.items():
            if units[key] != expected:
                errors.append("units.%s must be %s." % (key, expected))

    _validate_forcing(case["forcing"], flow_count, errors, warnings)

    site = case["site"]
    if _shape(site, SITE_KEYS, SITE_KEYS, "site", errors):
        _number(site["depth_m"], "site.depth_m", errors, 0, strict_lower=True)
        for key in SITE_KEYS-{"depth_m"}:
            upper = 1 if key in {
                "particulate_organic_fraction",
                "advective_exchange_efficiency"} else None
            _number(site[key], "site."+key, errors, 0, upper)

    food = case.get("food")
    if food is not None and _shape(food, FOOD_KEYS, FOOD_KEYS, "food", errors):
        for key in FOOD_KEYS:
            upper = 1.0 if key in {
                "organic_carbon_fraction", "small_particle_fraction",
                "detritus_fraction_of_organic", "detritus_preference",
                "detritus_assimilation_multiplier"} else None
            _number(
                food[key], "food."+key, errors, 0, upper,
                strict_lower=upper is None)

    species = case["species"]
    if _shape(species, SPECIES_REQUIRED_KEYS, SPECIES_KEYS, "species", errors):
        _text(species["taxon"], "species.taxon", errors)
        if species["profile_status"] != "screening_prior":
            errors.append("species.profile_status must be screening_prior.")
        _number(species["mean_dry_tissue_g"],
                "species.mean_dry_tissue_g", errors, 0, strict_lower=True)
        _number(species["size_cv"], "species.size_cv", errors, 0)
        for name, coefficient_key in (
                ("clearance", "a_l_h"), ("respiration", "a_ml_o2_h")):
            value = species[name]
            allowed = ALLOMETRY_COMMON_KEYS | {coefficient_key}
            if _shape(value, allowed, allowed, "species."+name, errors):
                _number(value[coefficient_key],
                        "species.%s.%s" % (name, coefficient_key),
                        errors, 0, strict_lower=True)
                _number(value["mass_exponent"],
                        "species.%s.mass_exponent" % name,
                        errors, 0, strict_lower=True)
                _number(value["reference_temperature_c"],
                        "species.%s.reference_temperature_c" % name, errors)
                _number(value["q10"], "species.%s.q10" % name,
                        errors, 0, strict_lower=True)
        for key in (
                "retention_efficiency", "particulate_retention_efficiency",
                "assimilation_efficiency", "pseudofaeces_fraction",
                "activity_fraction"):
            _number(species[key], "species."+key, errors, 0, 1)
        _number(species["ammonia_mg_n_g_dw_h"],
                "species.ammonia_mg_n_g_dw_h", errors, 0)
        _number(species["low_food_threshold_ug_l"],
                "species.low_food_threshold_ug_l", errors, 0)
        _number(species["low_food_transition_ug_l"],
                "species.low_food_transition_ug_l", errors, 0,
                strict_lower=True)
        _number(species["valid_flow_max_m_s"],
                "species.valid_flow_max_m_s", errors, 0, strict_lower=True)
        for key in (
                "retention_d50_um", "retention_slope_um",
                "ingestion_half_saturation_ug_l",
                "high_food_assimilation_decay_ug_l",
                "pseudofaeces_tsm_transition_mg_l",
                "current_clearance_zero_m_s", "current_protection_group_size",
                "oxycalorific_kj_per_ml_o2", "tissue_energy_kj_g_dw"):
            if key in species:
                _number(species[key], "species."+key, errors, 0,
                        strict_lower=True)
        for key in (
                "high_food_assimilation_threshold_ug_l",
                "pseudofaeces_tsm_threshold_mg_l",
                "current_clearance_start_m_s"):
            if key in species:
                _number(species[key], "species."+key, errors, 0)
        for key in (
                "assimilation_quality_half_saturation",
                "assimilation_reference_organic_fraction",
                "pseudofaeces_max_fraction",
                "oxygen_zero_saturation_fraction",
                "oxygen_full_saturation_fraction"):
            if key in species:
                _number(species[key], "species."+key, errors, 0, 1)
        if ("oxygen_zero_saturation_fraction" in species and
                "oxygen_full_saturation_fraction" in species and
                species["oxygen_zero_saturation_fraction"] >=
                species["oxygen_full_saturation_fraction"]):
            errors.append(
                "species.oxygen_zero_saturation_fraction must be smaller "
                "than species.oxygen_full_saturation_fraction.")
        if ("current_clearance_start_m_s" in species and
                "current_clearance_zero_m_s" in species and
                species["current_clearance_start_m_s"] >=
                species["current_clearance_zero_m_s"]):
            errors.append(
                "species.current_clearance_start_m_s must be smaller than "
                "species.current_clearance_zero_m_s.")
        if ("pseudofaeces_max_fraction" in species and
                species["pseudofaeces_max_fraction"] <
                species["pseudofaeces_fraction"]):
            errors.append(
                "species.pseudofaeces_max_fraction must be >= "
                "species.pseudofaeces_fraction.")
        salinity = species["salinity_activity_psu"]
        if _shape(salinity, SALINITY_KEYS, SALINITY_KEYS,
                  "species.salinity_activity_psu", errors):
            values = [salinity[key] for key in (
                "zero_low", "full_low", "full_high", "zero_high")]
            for key, value in zip(
                    ("zero_low", "full_low", "full_high", "zero_high"), values):
                _number(value, "species.salinity_activity_psu."+key, errors, 0)
            if all(isinstance(value, (int, float)) for value in values):
                if values != sorted(values):
                    errors.append(
                        "species.salinity_activity_psu values must be ordered.")

    stocking = case["stocking"]
    if _shape(stocking, STOCKING_REQUIRED_KEYS, STOCKING_KEYS,
              "stocking", errors):
        if stocking["mode"] != "animals":
            errors.append("stocking.mode currently supports only animals.")
        _scalar_or_array(
            stocking["mussels_per_obstacle"],
            "stocking.mussels_per_obstacle", errors)
        _number(stocking["live_wet_g_per_individual"],
                "stocking.live_wet_g_per_individual",
                errors, 0, strict_lower=True)
        _number(stocking["annual_mortality_fraction"],
                "stocking.annual_mortality_fraction", errors, 0, 1)
        if "effective_aggregation_size" in stocking:
            _number(stocking["effective_aggregation_size"],
                    "stocking.effective_aggregation_size", errors, 0,
                    strict_lower=True)

    structure = case["structure"]
    if _shape(structure, STRUCTURE_KEYS, STRUCTURE_KEYS, "structure", errors):
        _scalar_or_array(
            structure["porosity_per_obstacle"],
            "structure.porosity_per_obstacle", errors, 0, 0.99)
        for key in STRUCTURE_KEYS-{"porosity_per_obstacle"}:
            _number(structure[key], "structure."+key,
                    errors, 0, strict_lower=True)

    hydro = case["hydrodynamics"]
    if _shape(hydro, HYDRO_KEYS, HYDRO_KEYS, "hydrodynamics", errors):
        for key in HYDRO_KEYS:
            if key == "minimum_speed_ratio":
                _number(hydro[key], "hydrodynamics."+key, errors, 0, 1)
            else:
                _number(hydro[key], "hydrodynamics."+key,
                        errors, 0, strict_lower=True)

    sediment = case["sediment"]
    if _shape(sediment, SEDIMENT_KEYS, SEDIMENT_KEYS, "sediment", errors):
        _number(sediment["in_domain_deposition_fraction"],
                "sediment.in_domain_deposition_fraction", errors, 0, 1,
                nullable=True)
        _number(sediment["mortality_deposition_fraction"],
                "sediment.mortality_deposition_fraction", errors, 0, 1)
        for key in (
                "settling_velocity_m_s", "oxygen_demand_kg_o2_per_kg_organic",
                "decay_per_day", "resuspension_per_day",
                "initial_organic_stock_kg"):
            _number(sediment[key], "sediment."+key, errors, 0)

    harvest = case["harvest"]
    if _shape(harvest, HARVEST_KEYS, HARVEST_KEYS, "harvest", errors):
        _number(harvest["fraction_per_year"],
                "harvest.fraction_per_year", errors, 0, 1)
        for key in ("turnovers_per_year", "n_kg_per_t_wet", "p_kg_per_t_wet"):
            _number(harvest[key], "harvest."+key, errors, 0)

    constraints = case["constraints"]
    if _shape(constraints, CONSTRAINT_KEYS, CONSTRAINT_KEYS,
              "constraints", errors):
        for key in CONSTRAINT_KEYS:
            _number(constraints[key], "constraints."+key, errors, 0,
                    nullable=key != "minimum_obstacle_clearance_m")

    objective = case["objective"]
    if _shape(objective, OBJECTIVE_KEYS, OBJECTIVE_KEYS, "objective", errors):
        weights = objective["weights"]
        if _shape(weights, WEIGHT_KEYS, WEIGHT_KEYS,
                  "objective.weights", errors):
            for key in WEIGHT_KEYS:
                _number(weights[key], "objective.weights."+key, errors, 0)
            if all(isinstance(weights[key], (int, float))
                   for key in WEIGHT_KEYS) and sum(weights.values()) <= 0:
                errors.append("At least one objective weight must be positive.")
        _number(objective["target_chlorophyll_capture_g_day"],
                "objective.target_chlorophyll_capture_g_day",
                errors, 0, nullable=True)
        _number(objective["target_biodeposition_kg_m2_day"],
                "objective.target_biodeposition_kg_m2_day",
                errors, 0, strict_lower=True)

    validation = case["validation"]
    if _shape(validation, VALIDATION_KEYS, VALIDATION_KEYS,
              "validation", errors):
        if validation["status"] != RUNTIME_STATUS:
            errors.append(
                "validation.status is derived and must remain %s."
                % RUNTIME_STATUS)
        manifest = validation["evidence_manifest_id"]
        if manifest is not None:
            _text(manifest, "validation.evidence_manifest_id", errors)
            warnings.append(
                "evidence_manifest_id is recorded but this screening model "
                "does not verify it.")
        if not isinstance(validation["allow_extrapolation"], bool):
            errors.append("validation.allow_extrapolation must be boolean.")

    ml = case.get("machine_learning")
    if ml is not None and _shape(
            ml, ML_KEYS, ML_KEYS, "machine_learning", errors):
        warnings.append(
            "machine_learning is legacy inert metadata and does not affect "
            "the ecological calculation.")
        if ml["surrogate_enabled"] is not False:
            errors.append(
                "No validated surrogate is active; surrogate_enabled must "
                "remain false.")
        if ml["model_artifact_id"] is not None:
            errors.append("model_artifact_id must be null.")
        if ml["training_dataset_id"] is not None:
            errors.append("training_dataset_id must be null.")
        rl = ml["reinforcement_learning"]
        if _shape(rl, RL_KEYS, RL_KEYS,
                  "machine_learning.reinforcement_learning", errors):
            if rl["enabled"] is not False:
                errors.append("Reinforcement learning is not active.")
            _text(rl["intended_future_role"],
                  "machine_learning.reinforcement_learning.intended_future_role",
                  errors)
            _text(rl["static_layout_optimizer"],
                  "machine_learning.reinforcement_learning.static_layout_optimizer",
                  errors)

    evidence = case["evidence"]
    if not isinstance(evidence, list):
        errors.append("evidence must be an array.")
    else:
        for index, item in enumerate(evidence):
            path = "evidence[%d]" % index
            required_evidence = (
                EVIDENCE_KEYS if schema_version == SCHEMA_VERSION
                else LEGACY_EVIDENCE_KEYS)
            if _shape(item, required_evidence, EVIDENCE_KEYS, path, errors):
                _text(item["parameter"], path+".parameter", errors)
                for key in (
                        "citation", "equation", "validity_envelope"):
                    if key in item:
                        _text(item[key], path+"."+key, errors)
                _text(item["doi"], path+".doi", errors)
                if (isinstance(item["doi"], str) and
                        not item["doi"].startswith("10.")):
                    errors.append(path+".doi must begin with 10.")

    if obstacle_count is not None:
        _integer(obstacle_count, "obstacle_count", errors, 1)
        if isinstance(obstacle_count, int) and obstacle_count >= 1:
            for path, value in (
                    ("stocking.mussels_per_obstacle",
                     stocking.get("mussels_per_obstacle")
                     if isinstance(stocking, dict) else None),
                    ("structure.porosity_per_obstacle",
                     structure.get("porosity_per_obstacle")
                     if isinstance(structure, dict) else None)):
                if isinstance(value, list) and len(value) not in {
                        1, obstacle_count}:
                    errors.append(
                        "%s has %d values; expected 1 or %d."
                        % (path, len(value), obstacle_count))
    return errors, warnings


def parse_case(source, flow_count=None, obstacle_count=None):
    """Return ``(validated_case, warnings)`` or raise CaseValidationError."""
    case = _load_source(source)
    errors, warnings = _validate_case(case, flow_count, obstacle_count)
    if errors:
        raise CaseValidationError(errors)
    return case, warnings


def canonical_json(case):
    """Return deterministic compact JSON, rejecting NaN and infinity."""
    return json.dumps(
        case, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
        allow_nan=False)


def case_hash(case):
    """Return a SHA-256 hash of the canonical case JSON."""
    return hashlib.sha256(canonical_json(case).encode("utf-8")).hexdigest()


def _as_list(value, count, path):
    values = list(value) if isinstance(value, list) else [value]
    if len(values) == 1:
        return values*count
    if len(values) != count:
        raise CaseValidationError([
            "%s has %d values; expected 1 or %d."
            % (path, len(values), count)])
    return values


# Config keys whose values come from a forcing boundary block, and the boundary
# field that fills each one.  Timeline sets these once from its single shared
# boundary; the ensemble evaluator overrides them per independent state.
_BOUNDARY_TO_CONFIG = {
    "temperature_c": "site.temperature_c",
    "salinity_psu": "site.salinity_psu",
    "dissolved_oxygen_mg_l": "site.boundary_do_mg_l",
    "chlorophyll_a_ug_l": "site.chlorophyll_ug_l",
    "tsm_mg_l": "site.tsm_mg_l",
}


def _apply_boundary(config, boundary):
    """Write a forcing boundary block into the five boundary config keys."""
    for field, key in _BOUNDARY_TO_CONFIG.items():
        config[key] = boundary[field]
    return config


def _compile_shared(case, obstacle_count):
    """Build the flat solver config for every field except the per-scenario
    boundary chemistry and the scenario timing.

    The five boundary keys (temperature, salinity, boundary DO, chlorophyll,
    TSM) and the ``scenario.*`` keys keep their defaults so the caller can fill
    them: :func:`compile_timeline` from a single shared timeline boundary,
    or :func:`ensemble_state_configs` per independent ensemble state.
    """
    try:
        from .ecogrammar_core import DEFAULTS
    except (ImportError, ValueError):
        from musselflow_ecogrammar_core import DEFAULTS

    config = copy.deepcopy(DEFAULTS)
    site = case["site"]
    food = case.get("food", {})
    species = case["species"]
    clearance = species["clearance"]
    respiration = species["respiration"]
    salinity = species["salinity_activity_psu"]
    stocking = case["stocking"]
    structure = case["structure"]
    hydro = case["hydrodynamics"]
    sediment = case["sediment"]
    harvest = case["harvest"]
    constraints = case["constraints"]
    objective = case["objective"]

    config.update({
        "site.depth_m": site["depth_m"],
        "site.initial_do_mg_l": site["initial_dissolved_oxygen_mg_l"],
        "site.particulate_organic_fraction":
            site["particulate_organic_fraction"],
        "site.background_sod_g_o2_m2_day":
            site["background_sediment_oxygen_demand_g_o2_m2_day"],
        "site.pelagic_respiration_g_o2_m3_day":
            site["pelagic_respiration_g_o2_m3_day"],
        "site.primary_production_g_o2_m3_day":
            site["primary_production_g_o2_m3_day"],
        "site.reaeration_per_day": site["reaeration_per_day"],
        "site.vertical_exchange_per_day": site["vertical_exchange_per_day"],
        "site.advective_exchange_efficiency":
            site["advective_exchange_efficiency"],
        "food.carbon_to_chlorophyll_mg_c_per_mg_chl": food.get(
            "carbon_to_chlorophyll_mg_c_per_mg_chl",
            config["food.carbon_to_chlorophyll_mg_c_per_mg_chl"]),
        "food.organic_carbon_fraction": food.get(
            "organic_carbon_fraction",
            config["food.organic_carbon_fraction"]),
        "food.small_particle_fraction": food.get(
            "small_particle_fraction",
            config["food.small_particle_fraction"]),
        "food.small_particle_diameter_um": food.get(
            "small_particle_diameter_um",
            config["food.small_particle_diameter_um"]),
        "food.large_particle_diameter_um": food.get(
            "large_particle_diameter_um",
            config["food.large_particle_diameter_um"]),
        "food.detritus_fraction_of_organic": food.get(
            "detritus_fraction_of_organic",
            config["food.detritus_fraction_of_organic"]),
        "food.detritus_preference": food.get(
            "detritus_preference", config["food.detritus_preference"]),
        "food.detritus_assimilation_multiplier": food.get(
            "detritus_assimilation_multiplier",
            config["food.detritus_assimilation_multiplier"]),
        "food.organic_energy_kj_g": food.get(
            "organic_energy_kj_g", config["food.organic_energy_kj_g"]),
        "species.mean_dry_tissue_g": species["mean_dry_tissue_g"],
        "species.size_cv": species["size_cv"],
        "species.clearance_a_l_h": clearance["a_l_h"],
        "species.clearance_b": clearance["mass_exponent"],
        "species.clearance_ref_temp_c":
            clearance["reference_temperature_c"],
        "species.clearance_q10": clearance["q10"],
        "species.respiration_a_ml_o2_h": respiration["a_ml_o2_h"],
        "species.respiration_b": respiration["mass_exponent"],
        "species.respiration_ref_temp_c":
            respiration["reference_temperature_c"],
        "species.respiration_q10": respiration["q10"],
        "species.retention_efficiency": species["retention_efficiency"],
        "species.particulate_retention_efficiency":
            species["particulate_retention_efficiency"],
        "species.assimilation_efficiency":
            species["assimilation_efficiency"],
        "species.pseudofaeces_fraction": species["pseudofaeces_fraction"],
        "species.activity_fraction": species["activity_fraction"],
        "species.ammonia_mg_n_g_dw_h": species["ammonia_mg_n_g_dw_h"],
        "species.low_food_threshold_ug_l":
            species["low_food_threshold_ug_l"],
        "species.low_food_transition_ug_l":
            species["low_food_transition_ug_l"],
        "species.retention_d50_um": species.get(
            "retention_d50_um", config["species.retention_d50_um"]),
        "species.retention_slope_um": species.get(
            "retention_slope_um", config["species.retention_slope_um"]),
        "species.ingestion_half_saturation_ug_l": species.get(
            "ingestion_half_saturation_ug_l",
            config["species.ingestion_half_saturation_ug_l"]),
        "species.high_food_assimilation_threshold_ug_l": species.get(
            "high_food_assimilation_threshold_ug_l",
            config["species.high_food_assimilation_threshold_ug_l"]),
        "species.high_food_assimilation_decay_ug_l": species.get(
            "high_food_assimilation_decay_ug_l",
            config["species.high_food_assimilation_decay_ug_l"]),
        "species.assimilation_quality_half_saturation": species.get(
            "assimilation_quality_half_saturation",
            config["species.assimilation_quality_half_saturation"]),
        "species.assimilation_reference_organic_fraction": species.get(
            "assimilation_reference_organic_fraction",
            config["species.assimilation_reference_organic_fraction"]),
        "species.pseudofaeces_tsm_threshold_mg_l": species.get(
            "pseudofaeces_tsm_threshold_mg_l",
            config["species.pseudofaeces_tsm_threshold_mg_l"]),
        "species.pseudofaeces_tsm_transition_mg_l": species.get(
            "pseudofaeces_tsm_transition_mg_l",
            config["species.pseudofaeces_tsm_transition_mg_l"]),
        "species.pseudofaeces_max_fraction": species.get(
            "pseudofaeces_max_fraction",
            config["species.pseudofaeces_max_fraction"]),
        "species.oxygen_zero_saturation_fraction": species.get(
            "oxygen_zero_saturation_fraction",
            config["species.oxygen_zero_saturation_fraction"]),
        "species.oxygen_full_saturation_fraction": species.get(
            "oxygen_full_saturation_fraction",
            config["species.oxygen_full_saturation_fraction"]),
        "species.current_clearance_start_m_s": species.get(
            "current_clearance_start_m_s",
            config["species.current_clearance_start_m_s"]),
        "species.current_clearance_zero_m_s": species.get(
            "current_clearance_zero_m_s",
            config["species.current_clearance_zero_m_s"]),
        "species.current_protection_group_size": species.get(
            "current_protection_group_size",
            config["species.current_protection_group_size"]),
        "species.oxycalorific_kj_per_ml_o2": species.get(
            "oxycalorific_kj_per_ml_o2",
            config["species.oxycalorific_kj_per_ml_o2"]),
        "species.tissue_energy_kj_g_dw": species.get(
            "tissue_energy_kj_g_dw",
            config["species.tissue_energy_kj_g_dw"]),
        "species.valid_flow_max_m_s": species["valid_flow_max_m_s"],
        "species.salinity_zero_low_psu": salinity["zero_low"],
        "species.salinity_full_low_psu": salinity["full_low"],
        "species.salinity_full_high_psu": salinity["full_high"],
        "species.salinity_zero_high_psu": salinity["zero_high"],
        "stocking.mode": "animals",
        "stocking.mussels_per_obstacle": _as_list(
            stocking["mussels_per_obstacle"], obstacle_count,
            "stocking.mussels_per_obstacle"),
        "stocking.live_wet_g_per_individual":
            stocking["live_wet_g_per_individual"],
        "stocking.annual_mortality_fraction":
            stocking["annual_mortality_fraction"],
        "stocking.effective_aggregation_size": stocking.get(
            "effective_aggregation_size",
            config["stocking.effective_aggregation_size"]),
        "structure.porosity": _as_list(
            structure["porosity_per_obstacle"], obstacle_count,
            "structure.porosity_per_obstacle"),
        "structure.twine_diameter_m": structure["twine_diameter_m"],
        "structure.drag_multiplier": structure["drag_multiplier"],
        "structure.fallback_plan_size_m": structure["fallback_plan_size_m"],
        "structure.fallback_height_m": structure["fallback_height_m"],
        "hydrodynamics.wake_spread": hydro["wake_spread"],
        "hydrodynamics.food_plume_spread": hydro["food_plume_spread"],
        "hydrodynamics.food_recovery_lengths":
            hydro["food_recovery_lengths"],
        "hydrodynamics.min_speed_ratio": hydro["minimum_speed_ratio"],
        "hydrodynamics.kinematic_viscosity_m2_s":
            hydro["kinematic_viscosity_m2_s"],
        "hydrodynamics.water_density_kg_m3":
            hydro["water_density_kg_m3"],
        "sediment.settling_velocity_m_s":
            sediment["settling_velocity_m_s"],
        "sediment.in_domain_deposition_fraction":
            -1.0 if sediment["in_domain_deposition_fraction"] is None
            else sediment["in_domain_deposition_fraction"],
        "sediment.oxygen_demand_kg_o2_per_kg_organic":
            sediment["oxygen_demand_kg_o2_per_kg_organic"],
        "sediment.decay_per_day": sediment["decay_per_day"],
        "sediment.resuspension_per_day": sediment["resuspension_per_day"],
        "sediment.mortality_deposition_fraction":
            sediment["mortality_deposition_fraction"],
        "sediment.initial_organic_stock_kg":
            sediment["initial_organic_stock_kg"],
        "harvest.fraction_per_year": harvest["fraction_per_year"],
        "harvest.turnovers_per_year": harvest["turnovers_per_year"],
        "harvest.n_kg_per_t_wet": harvest["n_kg_per_t_wet"],
        "harvest.p_kg_per_t_wet": harvest["p_kg_per_t_wet"],
        "constraint.min_do_mg_l":
            -1.0 if constraints["minimum_dissolved_oxygen_mg_l"] is None
            else constraints["minimum_dissolved_oxygen_mg_l"],
        "constraint.min_probe_speed_m_s":
            -1.0 if constraints["minimum_probe_speed_m_s"] is None
            else constraints["minimum_probe_speed_m_s"],
        "constraint.max_biodeposition_kg_m2_day":
            -1.0 if constraints["maximum_biodeposition_kg_m2_day"] is None
            else constraints["maximum_biodeposition_kg_m2_day"],
        "constraint.min_obstacle_clearance_m":
            constraints["minimum_obstacle_clearance_m"],
        "objective.target_chlorophyll_capture_g_day":
            -1.0 if objective["target_chlorophyll_capture_g_day"] is None
            else objective["target_chlorophyll_capture_g_day"],
        "objective.target_biodeposition_kg_m2_day":
            objective["target_biodeposition_kg_m2_day"],
        "validation.geometry_calibrated": False,
        "validation.hydrodynamics_calibrated": False,
        "validation.biology_calibrated": False,
        "validation.oxygen_calibrated": False,
        "validation.sediment_calibrated": False,
        "validation.allow_extrapolation":
            case["validation"]["allow_extrapolation"],
    })
    for name, value in objective["weights"].items():
        config["weight."+name] = value
    return config


def compile_timeline(case, obstacle_count, flow_count):
    """Compile a timeline case into the numerical configuration.

    Returns ``(config, flow_indices, warnings)``.  Ensemble forcing is handled
    by :func:`compile_ensemble`; scenario-varying timeline boundary chemistry
    is still refused because the current oxygen/sediment sequence carries one
    shared boundary.
    """
    case, warnings = parse_case(
        case, flow_count=flow_count, obstacle_count=obstacle_count)
    forcing = case["forcing"]
    if forcing["mode"] != "timeline":
        raise UnsupportedCaseError(
            "compile_timeline handles timelines; use compile_ensemble "
            "for ensemble forcing.")
    steps = sorted(forcing["steps"], key=lambda item: item["order"])
    boundaries = [canonical_json(step["boundary"]) for step in steps]
    if len(set(boundaries)) != 1:
        raise UnsupportedCaseError(
            "The current oxygen/sediment sequence supports one shared boundary; "
            "scenario-varying timeline boundaries need the per-step evaluator.")

    try:
        from .ecogrammar_core import validate_config
    except (ImportError, ValueError):
        from musselflow_ecogrammar_core import validate_config

    config = _compile_shared(case, obstacle_count)
    _apply_boundary(config, steps[0]["boundary"])
    config["scenario.duration_h"] = [step["duration_h"] for step in steps]
    config["scenario.repeat_count"] = forcing["repeat_count"]
    durations = config["scenario.duration_h"]
    duration_sum = sum(durations)
    config["scenario.weights"] = [
        duration/duration_sum for duration in durations]

    flow_indices = [step["flow_vector_index"] for step in steps]
    adapter_warnings, adapter_errors = validate_config(
        config, obstacle_count=obstacle_count, scenario_count=len(steps))
    if adapter_errors:
        raise CaseValidationError(adapter_errors)
    warnings.extend(adapter_warnings)
    warnings.append(
        "Runtime status remains UNVALIDATED_SCREENING; this compiler does not "
        "verify calibration evidence.")
    return config, flow_indices, warnings


def compile_legacy_config(case, obstacle_count, flow_count):
    """Compatibility alias for saved scripts; use :func:`compile_timeline`."""
    return compile_timeline(case, obstacle_count, flow_count)


def compile_ensemble(case, obstacle_count, flow_count):
    """Compile an ensemble case into ``(base_config, states, warnings)``.

    Each ensemble state is an independent environmental condition evaluated on
    its own boundary chemistry.  ``base_config`` holds every field except the
    five boundary keys (left at defaults); :func:`ensemble_state_configs` fills
    them per state.  The evaluator combines objectives by occurrence
    probability and takes hard constraints from the worst state, so ensembles
    no longer need the fake sequential semantics the timeline bridge refused.
    """
    case, warnings = parse_case(
        case, flow_count=flow_count, obstacle_count=obstacle_count)
    forcing = case["forcing"]
    if forcing["mode"] != "ensemble":
        raise UnsupportedCaseError(
            "compile_ensemble handles ensemble forcing; use "
            "compile_timeline for timelines.")

    try:
        from .ecogrammar_core import validate_config
    except (ImportError, ValueError):
        from musselflow_ecogrammar_core import validate_config

    base_config = _compile_shared(case, obstacle_count)
    states = [
        {
            "id": state["id"],
            "flow_vector_index": state["flow_vector_index"],
            "occurrence_probability": state["occurrence_probability"],
            "duration_h": state["duration_h"],
            "boundary": state["boundary"],
        }
        for state in forcing["states"]]

    # Validate the shared (non-boundary) config once, using the first state's
    # boundary and a single scenario as representative; every state's boundary
    # ranges were already checked by parse_case.
    probe_config = _apply_boundary(
        copy.deepcopy(base_config), states[0]["boundary"])
    probe_config["scenario.duration_h"] = [states[0]["duration_h"]]
    probe_config["scenario.weights"] = [1.0]
    probe_config["scenario.repeat_count"] = 1
    adapter_warnings, adapter_errors = validate_config(
        probe_config, obstacle_count=obstacle_count, scenario_count=1)
    if adapter_errors:
        raise CaseValidationError(adapter_errors)
    warnings.extend(adapter_warnings)
    warnings.append(
        "Ensemble forcing: each state is evaluated independently under its own "
        "boundary; objectives are probability-weighted and hard constraints "
        "use the worst state.")
    warnings.append(
        "Runtime status remains UNVALIDATED_SCREENING; this compiler does not "
        "verify calibration evidence.")
    return base_config, states, warnings


def ensemble_state_configs(base_config, states):
    """Return per-state evaluation specs for :func:`evaluate_ensemble`.

    Each item is ``{"id", "flow_vector_index", "probability", "config"}`` where
    ``config`` is a deep copy of ``base_config`` with that state's boundary
    chemistry applied and a single-scenario window of the state's duration.
    """
    resolved = []
    for state in states:
        config = _apply_boundary(
            copy.deepcopy(base_config), state["boundary"])
        config["scenario.duration_h"] = [state["duration_h"]]
        config["scenario.weights"] = [1.0]
        config["scenario.repeat_count"] = 1
        resolved.append({
            "id": state["id"],
            "flow_vector_index": state["flow_vector_index"],
            "probability": state["occurrence_probability"],
            "config": config,
        })
    return resolved
