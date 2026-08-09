"""Strict ecological-grammar parser for the MusselFlow screening model.

Grasshopper Panels can supply either one multiline string or a list of lines:

    site.depth_m = 12
    site.temperature_c = 12
    structure.porosity = 0.70
    scenario.duration_h = 6, 6

The parser deliberately accepts a small, explicit language.  Natural-language
notes are allowed under ``note.*`` and ``source.*`` keys, but they never alter
the calculation.  Unknown computational keys are errors; this prevents a typo
from silently changing the meaning of an optimization run.

This module has no Rhino dependency and is suitable for unit testing.
"""

from __future__ import annotations

import copy
import math
import re


SCHEMA_VERSION = 1


# These are screening priors, not universal biological constants.  Every value
# that materially affects a result is exposed in the grammar and reported.
DEFAULTS = {
    "schema.version": SCHEMA_VERSION,
    "profile.species": "Mytilus_edulis_screening",
    "profile.net": "knotless_screening",

    # Site state and water-column box.
    "site.depth_m": 12.0,
    "site.temperature_c": 12.0,
    "site.salinity_psu": 20.0,
    "site.initial_do_mg_l": 9.0,
    "site.boundary_do_mg_l": 9.0,
    "site.chlorophyll_ug_l": 5.0,
    "site.tsm_mg_l": 3.0,
    "site.particulate_organic_fraction": 0.35,
    "site.background_sod_g_o2_m2_day": 0.8,
    "site.pelagic_respiration_g_o2_m3_day": 0.05,
    "site.primary_production_g_o2_m3_day": 0.0,
    "site.reaeration_per_day": 0.20,
    "site.vertical_exchange_per_day": 0.50,
    "site.advective_exchange_efficiency": 0.25,

    # Food composition and particle classes. These remain site-calibration
    # priors; chlorophyll alone is not an energy or nutrient measurement.
    "food.carbon_to_chlorophyll_mg_c_per_mg_chl": 50.0,
    "food.organic_carbon_fraction": 0.45,
    "food.small_particle_fraction": 0.20,
    "food.small_particle_diameter_um": 2.0,
    "food.large_particle_diameter_um": 15.0,
    "food.detritus_fraction_of_organic": 0.50,
    "food.detritus_preference": 0.80,
    "food.detritus_assimilation_multiplier": 0.80,
    "food.organic_energy_kj_g": 20.0,

    # Adult Mytilus edulis priors at approximately 10-13 C and 30 psu.
    "species.mean_dry_tissue_g": 0.20,
    "species.size_cv": 0.30,
    "species.clearance_a_l_h": 7.45,
    "species.clearance_b": 0.66,
    "species.clearance_ref_temp_c": 12.0,
    "species.clearance_q10": 1.50,
    "species.respiration_a_ml_o2_h": 0.475,
    "species.respiration_b": 0.663,
    "species.respiration_ref_temp_c": 12.0,
    "species.respiration_q10": 2.00,
    "species.retention_efficiency": 0.80,
    "species.particulate_retention_efficiency": 0.80,
    "species.assimilation_efficiency": 0.70,
    "species.pseudofaeces_fraction": 0.10,
    "species.ammonia_mg_n_g_dw_h": 0.015,
    "species.activity_fraction": 1.00,
    "species.low_food_threshold_ug_l": 0.70,
    "species.low_food_transition_ug_l": 0.15,
    "species.retention_d50_um": 5.0,
    "species.retention_slope_um": 1.5,
    "species.ingestion_half_saturation_ug_l": 17.0,
    "species.high_food_assimilation_threshold_ug_l": 17.0,
    "species.high_food_assimilation_decay_ug_l": 8.0,
    "species.assimilation_quality_half_saturation": 0.15,
    "species.assimilation_reference_organic_fraction": 0.35,
    "species.pseudofaeces_tsm_threshold_mg_l": 4.0,
    "species.pseudofaeces_tsm_transition_mg_l": 2.0,
    "species.pseudofaeces_max_fraction": 0.60,
    "species.oxygen_zero_saturation_fraction": 0.20,
    "species.oxygen_full_saturation_fraction": 0.70,
    "species.current_clearance_start_m_s": 0.20,
    "species.current_clearance_zero_m_s": 0.60,
    "species.current_protection_group_size": 20.0,
    "species.oxycalorific_kj_per_ml_o2": 0.0201,
    "species.tissue_energy_kj_g_dw": 23.5,
    # A broad, editable screening envelope.  Local populations must replace it.
    "species.salinity_zero_low_psu": 4.0,
    "species.salinity_full_low_psu": 15.0,
    "species.salinity_full_high_psu": 32.0,
    "species.salinity_zero_high_psu": 42.0,
    "species.valid_flow_max_m_s": 1.40,

    # Stocking.  A scalar is broadcast; a comma list maps by obstacle index.
    "stocking.mode": "animals",
    "stocking.mussels_per_obstacle": [1000.0],
    "stocking.dry_tissue_kg_per_obstacle": [-1.0],
    "stocking.live_wet_g_per_individual": 20.0,
    "stocking.annual_mortality_fraction": 0.10,
    "stocking.effective_aggregation_size": 20.0,

    # Porous structure and reduced wake.
    "structure.porosity": [0.70],
    "structure.twine_diameter_m": 0.003,
    "structure.drag_multiplier": 1.0,
    "structure.fallback_plan_size_m": 0.10,
    "structure.fallback_height_m": 1.00,
    "hydrodynamics.wake_spread": 0.12,
    "hydrodynamics.food_plume_spread": 0.16,
    "hydrodynamics.food_recovery_lengths": 8.0,
    "hydrodynamics.min_speed_ratio": 0.02,
    "hydrodynamics.kinematic_viscosity_m2_s": 1.30e-6,
    "hydrodynamics.water_density_kg_m3": 1020.0,

    # Flow vectors are supplied geometrically.  These lists define how the
    # vectors are combined and for how long oxygen/deposit states are advanced.
    "scenario.weights": [1.0],
    "scenario.duration_h": [24.0],
    "scenario.repeat_count": 1,

    # Organic deposition and oxygen demand.  These are especially
    # site-dependent and deliberately carry large validation warnings.
    "sediment.settling_velocity_m_s": 0.006,
    "sediment.in_domain_deposition_fraction": -1.0,
    "sediment.oxygen_demand_kg_o2_per_kg_organic": 1.0,
    "sediment.decay_per_day": 0.05,
    "sediment.resuspension_per_day": 0.0,
    "sediment.mortality_deposition_fraction": 1.0,
    "sediment.initial_organic_stock_kg": 0.0,

    # Harvest accounting is a management scenario, not predicted growth.
    "harvest.fraction_per_year": 0.80,
    "harvest.turnovers_per_year": 1.0,
    "harvest.n_kg_per_t_wet": 13.7,
    "harvest.p_kg_per_t_wet": 0.9,

    # Hard constraints.  Negative values disable a constraint.  No legal
    # threshold is silently asserted by the program.
    "constraint.min_do_mg_l": -1.0,
    "constraint.min_probe_speed_m_s": -1.0,
    "constraint.max_biodeposition_kg_m2_day": -1.0,
    "constraint.min_obstacle_clearance_m": 0.0,

    # Scalar fitness remains a user value judgement.  Physical outputs and
    # individual objective scores are also returned for Pareto workflows.
    "objective.target_chlorophyll_capture_g_day": -1.0,
    "objective.target_biodeposition_kg_m2_day": 0.05,
    "weight.extraction": 0.30,
    "weight.flushing": 0.20,
    "weight.food_delivery": 0.15,
    "weight.uniformity": 0.10,
    "weight.oxygen": 0.15,
    "weight.low_deposition": 0.10,

    # Governance / run envelope.
    "validation.geometry_calibrated": False,
    "validation.hydrodynamics_calibrated": False,
    "validation.biology_calibrated": False,
    "validation.oxygen_calibrated": False,
    "validation.sediment_calibrated": False,
    "validation.allow_extrapolation": False,
}


LIST_KEYS = {
    "stocking.mussels_per_obstacle",
    "stocking.dry_tissue_kg_per_obstacle",
    "structure.porosity",
    "scenario.weights",
    "scenario.duration_h",
}

STRING_KEYS = {
    "profile.species",
    "profile.net",
    "stocking.mode",
}

BOOL_KEYS = {
    "validation.geometry_calibrated",
    "validation.hydrodynamics_calibrated",
    "validation.biology_calibrated",
    "validation.oxygen_calibrated",
    "validation.sediment_calibrated",
    "validation.allow_extrapolation",
}

INT_KEYS = {
    "schema.version",
    "scenario.repeat_count",
}


PROVENANCE = {
    "species.clearance": {
        "citation": "Moehlenberg and Riisgaard (1979), adult Mytilus edulis",
        "equation": "clearance_L_h = 7.45 * dry_tissue_g ** 0.66",
        "envelope": "adult dry tissue 0.011-1.361 g; approximately 10-13 C, 30 psu",
        "url": "https://doi.org/10.1242/bio.062024",
    },
    "species.respiration": {
        "citation": "Hamburger et al. (1983), adult Mytilus edulis",
        "equation": "respiration_mL_O2_h = 0.475 * dry_tissue_g ** 0.663",
        "envelope": "adult dry tissue; approximately 10-13 C, 30 psu",
        "url": "https://doi.org/10.1242/bio.062024",
    },
    "species.current": {
        "citation": "Nielsen and Vismann (2014)",
        "equation": "No universal multiplier; model reports values above study envelope",
        "envelope": "0.05-1.4 m/s; aggregation strongly changed response",
        "url": (
            "https://orbit.dtu.dk/en/publications/"
            "clearance-rate-of-mytilus-edulis-l-as-a-function-of-current-veloc/"
        ),
    },
    "species.retention": {
        "citation": "Strohmeier et al. (2012)",
        "equation": "Direct editable retention efficiency",
        "envelope": "Natural seston; strong particle-size and temporal variability",
        "url": "https://doi.org/10.1016/j.jembe.2011.11.006",
    },
    "structure.drag": {
        "citation": "Lader et al. (2009), Aquaculture Net Drag Force and Added Mass",
        "equation": "Cd_solid = 1 + 1.37*S + 0.78*S^2 for knotless net",
        "envelope": "Net-panel empirical relation; S is solidity",
        "url": "https://doi.org/10.1016/j.aquaeng.2009.04.003",
    },
    "oxygen.solubility": {
        "citation": "Garcia and Gordon (1992)",
        "equation": "O2 saturation from temperature and salinity",
        "envelope": "-2 to 40 C and 0 to 42 salinity",
        "url": "https://doi.org/10.4319/lo.1992.37.6.1307",
    },
    "nutrient.turnover": {
        "citation": "Jansen et al. (2012)",
        "equation": "Direct editable excretion and tissue-accounting priors",
        "envelope": "Annual study; rates covaried with food and temperature",
        "url": "https://doi.org/10.1016/j.jembe.2011.11.009",
    },
    "harvest.content": {
        "citation": "Danish mitigation-mussel demonstration values",
        "equation": "13.7 kg N and 0.9 kg P per harvested wet tonne",
        "envelope": "Management accounting prior; replace with measured harvest composition",
        "url": "https://doi.org/10.3389/fmars.2019.00698",
    },
}


_TRUE = {"true", "yes", "on", "1"}
_FALSE = {"false", "no", "off", "0"}
_KEY_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


def _normalise_lines(value):
    """Return Panel-like input as a flat list of text lines."""
    if value is None:
        return []
    if isinstance(value, str):
        items = [value]
    else:
        try:
            items = list(value)
        except TypeError:
            items = [value]
    lines = []
    for item in items:
        if item is None:
            continue
        lines.extend(str(item).replace("\r\n", "\n").replace("\r", "\n").split("\n"))
    return lines


def _strip_inline_comment(text):
    """Strip # comments while preserving quoted text."""
    quote = None
    output = []
    for character in text:
        if character in {"'", '"'}:
            if quote is None:
                quote = character
            elif quote == character:
                quote = None
        if character == "#" and quote is None:
            break
        output.append(character)
    return "".join(output).strip()


def _parse_bool(value):
    normalised = value.strip().lower()
    if normalised in _TRUE:
        return True
    if normalised in _FALSE:
        return False
    raise ValueError("expected true/false")


def _parse_number(value):
    number = float(value.strip())
    if not math.isfinite(number):
        raise ValueError("number must be finite")
    return number


def _parse_value(key, raw):
    raw = raw.strip()
    if key in STRING_KEYS:
        if ((raw.startswith('"') and raw.endswith('"')) or
                (raw.startswith("'") and raw.endswith("'"))):
            raw = raw[1:-1]
        if not raw:
            raise ValueError("text value cannot be empty")
        return raw
    if key in BOOL_KEYS:
        return _parse_bool(raw)
    if key in INT_KEYS:
        number = _parse_number(raw)
        if int(number) != number:
            raise ValueError("expected an integer")
        return int(number)
    if key in LIST_KEYS:
        parts = [part.strip() for part in raw.split(",")]
        if not parts or any(not part for part in parts):
            raise ValueError("expected a comma-separated numeric list")
        return [_parse_number(part) for part in parts]
    return _parse_number(raw)


def _check_range(config, key, lower=None, upper=None, strict_lower=False):
    value = config[key]
    values = value if isinstance(value, list) else [value]
    errors = []
    for item in values:
        if lower is not None:
            invalid = item <= lower if strict_lower else item < lower
            if invalid:
                operator = ">" if strict_lower else ">="
                errors.append("%s must be %s %s." % (key, operator, lower))
                break
        if upper is not None and item > upper:
            errors.append("%s must be <= %s." % (key, upper))
            break
    return errors


def broadcast(values, count, key):
    """Broadcast a one-item numeric list or verify an exact per-item list."""
    values = list(values)
    if count <= 0:
        return []
    if len(values) == 1:
        return values * count
    if len(values) != count:
        raise ValueError(
            "%s has %d values; expected 1 or %d." % (key, len(values), count))
    return values


def match_scenarios(values, count, key):
    """Broadcast or repeat a scenario list deterministically."""
    values = list(values)
    if count <= 0:
        return []
    if len(values) == 1:
        return values * count
    if len(values) != count:
        raise ValueError(
            "%s has %d values; expected 1 or %d flow scenarios."
            % (key, len(values), count))
    return values


def validate_config(config, obstacle_count=None, scenario_count=None):
    """Return ``(warnings, errors)`` for a parsed flat configuration."""
    warnings = []
    errors = []

    if config["schema.version"] != SCHEMA_VERSION:
        errors.append(
            "schema.version=%s is unsupported; this node requires %s."
            % (config["schema.version"], SCHEMA_VERSION))

    for key in (
            "site.depth_m",
            "species.mean_dry_tissue_g",
            "species.clearance_a_l_h",
            "species.clearance_b",
            "species.clearance_q10",
            "species.respiration_a_ml_o2_h",
            "species.respiration_b",
            "species.respiration_q10",
            "stocking.live_wet_g_per_individual",
            "structure.twine_diameter_m",
            "structure.drag_multiplier",
            "structure.fallback_plan_size_m",
            "structure.fallback_height_m",
            "hydrodynamics.wake_spread",
            "hydrodynamics.food_plume_spread",
            "hydrodynamics.food_recovery_lengths",
            "hydrodynamics.kinematic_viscosity_m2_s",
            "hydrodynamics.water_density_kg_m3",
            "species.low_food_transition_ug_l",
            "food.carbon_to_chlorophyll_mg_c_per_mg_chl",
            "food.small_particle_diameter_um",
            "food.large_particle_diameter_um",
            "food.organic_energy_kj_g",
            "species.retention_d50_um",
            "species.retention_slope_um",
            "species.ingestion_half_saturation_ug_l",
            "species.high_food_assimilation_decay_ug_l",
            "species.pseudofaeces_tsm_transition_mg_l",
            "species.current_clearance_zero_m_s",
            "species.current_protection_group_size",
            "species.oxycalorific_kj_per_ml_o2",
            "species.tissue_energy_kj_g_dw",
            "stocking.effective_aggregation_size",
            "species.valid_flow_max_m_s",
            "objective.target_biodeposition_kg_m2_day"):
        errors.extend(_check_range(config, key, 0.0, strict_lower=True))

    for key in (
            "site.salinity_psu",
            "site.initial_do_mg_l",
            "site.boundary_do_mg_l",
            "site.chlorophyll_ug_l",
            "site.tsm_mg_l",
            "site.background_sod_g_o2_m2_day",
            "site.pelagic_respiration_g_o2_m3_day",
            "site.primary_production_g_o2_m3_day",
            "site.reaeration_per_day",
            "site.vertical_exchange_per_day",
            "species.size_cv",
            "species.low_food_threshold_ug_l",
            "species.high_food_assimilation_threshold_ug_l",
            "species.pseudofaeces_tsm_threshold_mg_l",
            "species.current_clearance_start_m_s",
            "species.ammonia_mg_n_g_dw_h",
            "species.salinity_zero_low_psu",
            "species.salinity_full_low_psu",
            "species.salinity_full_high_psu",
            "species.salinity_zero_high_psu",
            "sediment.settling_velocity_m_s",
            "sediment.oxygen_demand_kg_o2_per_kg_organic",
            "sediment.decay_per_day",
            "sediment.resuspension_per_day",
            "sediment.initial_organic_stock_kg",
            "harvest.n_kg_per_t_wet",
            "harvest.p_kg_per_t_wet",
            "harvest.turnovers_per_year",
            "constraint.min_obstacle_clearance_m"):
        errors.extend(_check_range(config, key, 0.0))

    for key in (
            "site.particulate_organic_fraction",
            "site.advective_exchange_efficiency",
            "food.organic_carbon_fraction",
            "food.small_particle_fraction",
            "food.detritus_fraction_of_organic",
            "food.detritus_preference",
            "food.detritus_assimilation_multiplier",
            "species.retention_efficiency",
            "species.particulate_retention_efficiency",
            "species.assimilation_efficiency",
            "species.pseudofaeces_fraction",
            "species.activity_fraction",
            "species.assimilation_quality_half_saturation",
            "species.assimilation_reference_organic_fraction",
            "species.pseudofaeces_max_fraction",
            "species.oxygen_zero_saturation_fraction",
            "species.oxygen_full_saturation_fraction",
            "stocking.annual_mortality_fraction",
            "sediment.mortality_deposition_fraction",
            "harvest.fraction_per_year"):
        errors.extend(_check_range(config, key, 0.0, 1.0))
    errors.extend(_check_range(
        config, "hydrodynamics.min_speed_ratio", 0.0, 1.0))

    errors.extend(_check_range(config, "structure.porosity", 0.0, 0.99))
    errors.extend(_check_range(
        config, "stocking.mussels_per_obstacle", 0.0))
    errors.extend(_check_range(
        config, "stocking.dry_tissue_kg_per_obstacle", -1.0))
    errors.extend(_check_range(config, "scenario.weights", 0.0))
    errors.extend(_check_range(
        config, "scenario.duration_h", 0.0, strict_lower=True))
    deposition_fraction = config["sediment.in_domain_deposition_fraction"]
    if (deposition_fraction < -1.0 or deposition_fraction > 1.0):
        errors.append(
            "sediment.in_domain_deposition_fraction must be -1 (automatic) "
            "or between 0 and 1.")
    invalid_dry_values = [
        value for value in config["stocking.dry_tissue_kg_per_obstacle"]
        if -1.0 < value < 0.0]
    if invalid_dry_values:
        errors.append(
            "stocking.dry_tissue_kg_per_obstacle values must be -1 "
            "(unused) or non-negative.")

    if config["scenario.repeat_count"] < 1:
        errors.append("scenario.repeat_count must be >= 1.")

    if config["stocking.mode"].lower() not in {"animals", "dry_biomass"}:
        errors.append("stocking.mode must be animals or dry_biomass.")

    salinity_values = [
        config["species.salinity_zero_low_psu"],
        config["species.salinity_full_low_psu"],
        config["species.salinity_full_high_psu"],
        config["species.salinity_zero_high_psu"],
    ]
    if salinity_values != sorted(salinity_values):
        errors.append(
            "Species salinity breakpoints must be ordered zero_low <= "
            "full_low <= full_high <= zero_high.")

    if (config["species.oxygen_zero_saturation_fraction"] >=
            config["species.oxygen_full_saturation_fraction"]):
        errors.append(
            "species.oxygen_zero_saturation_fraction must be smaller than "
            "species.oxygen_full_saturation_fraction.")
    if (config["species.current_clearance_start_m_s"] >=
            config["species.current_clearance_zero_m_s"]):
        errors.append(
            "species.current_clearance_start_m_s must be smaller than "
            "species.current_clearance_zero_m_s.")
    if (config["species.pseudofaeces_max_fraction"] <
            config["species.pseudofaeces_fraction"]):
        errors.append(
            "species.pseudofaeces_max_fraction must be >= the baseline "
            "species.pseudofaeces_fraction.")

    weights = [
        config[key] for key in (
            "weight.extraction",
            "weight.flushing",
            "weight.food_delivery",
            "weight.uniformity",
            "weight.oxygen",
            "weight.low_deposition")
    ]
    if any(weight < 0.0 for weight in weights):
        errors.append("Objective weights must be non-negative.")
    if sum(weights) <= 0.0:
        errors.append("At least one objective weight must be positive.")

    if obstacle_count is not None:
        for key in (
                "stocking.mussels_per_obstacle",
                "stocking.dry_tissue_kg_per_obstacle",
                "structure.porosity"):
            try:
                broadcast(config[key], obstacle_count, key)
            except ValueError as exception:
                errors.append(str(exception))

    if scenario_count is not None:
        for key in ("scenario.weights", "scenario.duration_h"):
            try:
                match_scenarios(config[key], scenario_count, key)
            except ValueError as exception:
                errors.append(str(exception))

    if sum(config["scenario.weights"]) <= 0.0:
        errors.append("scenario.weights must contain at least one positive value.")

    if config["site.temperature_c"] < -2.0 or config["site.temperature_c"] > 40.0:
        warnings.append(
            "Temperature is outside the Garcia-Gordon oxygen-solubility envelope.")
    if config["site.salinity_psu"] > 42.0:
        warnings.append(
            "Salinity is outside the Garcia-Gordon oxygen-solubility envelope.")
    if config["site.salinity_psu"] < config["species.salinity_full_low_psu"]:
        warnings.append(
            "Salinity reduces activity under the generic profile; replace the "
            "profile with measurements for the local Mytilus population.")
    net_profile = config["profile.net"].lower()
    if "knotless" not in net_profile and "knotted" not in net_profile:
        warnings.append(
            "Unknown net profile name; the numerical core will fall back to "
            "the knotless drag relation.")

    for key in (
            "validation.geometry_calibrated",
            "validation.hydrodynamics_calibrated",
            "validation.biology_calibrated",
            "validation.oxygen_calibrated",
            "validation.sediment_calibrated"):
        if not config[key]:
            warnings.append("%s is false." % key)

    if config["constraint.min_do_mg_l"] < 0.0:
        warnings.append(
            "No minimum dissolved-oxygen constraint is active.")
    if config["constraint.max_biodeposition_kg_m2_day"] < 0.0:
        warnings.append(
            "No maximum biodeposition constraint is active.")
    if config["objective.target_chlorophyll_capture_g_day"] <= 0.0:
        warnings.append(
            "No absolute chlorophyll-capture target is set; extraction score "
            "will use the current design's theoretical clearance capacity.")
    if config["sediment.in_domain_deposition_fraction"] < 0.0:
        warnings.append(
            "In-domain deposition fraction will use a simple settling-distance "
            "screen, not a seabed transport model.")

    return warnings, errors


def parse_grammar(lines, obstacle_count=None, scenario_count=None):
    """Parse Panel text and return ``(config, warnings, errors, notes)``.

    ``note.*`` and ``source.*`` entries are preserved in ``notes`` but do not
    enter the computational configuration.
    """
    config = copy.deepcopy(DEFAULTS)
    warnings = []
    errors = []
    notes = {}
    seen = set()

    for line_number, source_line in enumerate(_normalise_lines(lines), 1):
        line = _strip_inline_comment(source_line).strip()
        if not line:
            continue
        if "=" not in line:
            errors.append(
                "Grammar line %d must use key = value: %s"
                % (line_number, source_line))
            continue
        key, raw = [part.strip() for part in line.split("=", 1)]
        if not _KEY_PATTERN.match(key):
            errors.append("Grammar line %d has an invalid key: %s"
                          % (line_number, key))
            continue
        if key.startswith("note.") or key.startswith("source."):
            notes[key] = raw.strip().strip("'\"")
            continue
        if key not in DEFAULTS:
            errors.append(
                "Grammar line %d uses unknown key '%s'." % (line_number, key))
            continue
        if key in seen:
            warnings.append(
                "Grammar key '%s' appears more than once; the last value wins."
                % key)
        seen.add(key)
        try:
            config[key] = _parse_value(key, raw)
        except (TypeError, ValueError) as exception:
            errors.append(
                "Grammar line %d (%s): %s" % (line_number, key, exception))

    validation_warnings, validation_errors = validate_config(
        config, obstacle_count=obstacle_count, scenario_count=scenario_count)
    warnings.extend(validation_warnings)
    errors.extend(validation_errors)
    return config, warnings, errors, notes


def resolved_lists(config, obstacle_count, scenario_count):
    """Return a copy with all obstacle/scenario lists explicitly expanded."""
    resolved = copy.deepcopy(config)
    for key in (
            "stocking.mussels_per_obstacle",
            "stocking.dry_tissue_kg_per_obstacle",
            "structure.porosity"):
        resolved[key] = broadcast(config[key], obstacle_count, key)
    for key in ("scenario.weights", "scenario.duration_h"):
        resolved[key] = match_scenarios(config[key], scenario_count, key)
    total = sum(resolved["scenario.weights"])
    if total <= 0.0:
        raise ValueError("scenario.weights must have positive sum.")
    resolved["scenario.weights"] = [
        value / total for value in resolved["scenario.weights"]]
    return resolved


def grammar_template():
    """Return a concise, editable starter grammar for a Grasshopper Panel."""
    return [
        "# Legacy compatibility panel grammar (use musselflow_grammar.json)",
        "schema.version = 1",
        "profile.species = Mytilus_edulis_screening",
        "profile.net = knotless_screening",
        "",
        "# Site and forcing (replace with measurements)",
        "site.depth_m = 12",
        "site.temperature_c = 12",
        "site.salinity_psu = 20",
        "site.initial_do_mg_l = 9",
        "site.boundary_do_mg_l = 9",
        "site.chlorophyll_ug_l = 5",
        "site.tsm_mg_l = 3",
        "site.particulate_organic_fraction = 0.35",
        "# Scalars broadcast to every connected flow vector",
        "scenario.weights = 1",
        "scenario.duration_h = 6",
        "",
        "# Stock and net; scalar values broadcast to every obstacle",
        "stocking.mode = animals",
        "stocking.mussels_per_obstacle = 1000",
        "species.mean_dry_tissue_g = 0.20",
        "species.size_cv = 0.30",
        "structure.porosity = 0.70",
        "structure.twine_diameter_m = 0.003",
        "structure.fallback_plan_size_m = 0.10",
        "structure.fallback_height_m = 1.00",
        "",
        "# Biological rates are priors until locally calibrated",
        "species.retention_efficiency = 0.80",
        "species.particulate_retention_efficiency = 0.80",
        "species.assimilation_efficiency = 0.70",
        "species.pseudofaeces_fraction = 0.10",
        "",
        "# Activate site/permit-specific hard limits explicitly",
        "constraint.min_do_mg_l = -1",
        "constraint.min_probe_speed_m_s = -1",
        "constraint.max_biodeposition_kg_m2_day = -1",
        "constraint.min_obstacle_clearance_m = 0",
        "",
        "# Optional absolute design target",
        "objective.target_chlorophyll_capture_g_day = -1",
        "",
        "# Calibration declarations remain false until evidence exists",
        "validation.geometry_calibrated = false",
        "validation.hydrodynamics_calibrated = false",
        "validation.biology_calibrated = false",
        "validation.oxygen_calibrated = false",
        "validation.sediment_calibrated = false",
    ]
