"""Fast, physics-guided ecological screening core for MusselFlow.

This is a reduced-order model for comparing candidate farm layouts.  It is not
Navier-Stokes CFD, a sediment diagenesis model, a population-growth model, or a
permit assessment.  It does, however, keep dimensions and mass-balance terms
explicit so that its outputs can be calibrated and falsified.

Numeric inputs
--------------
obstacles : (N, 6) array
    ``x_m, y_m, major_m, minor_m, yaw_rad, height_m``.
domain_polygon : (K, 2) array
    Closed or open list of planar boundary vertices in metres.
probes : (P, 2) array
    Probe coordinates in metres.  If empty, obstacle centres explicitly become
    the fallback probes.
flow_vectors : (S, 2) array
    Uniform-current scenarios.  Vector length is speed in m/s and direction is
    plan-view direction.  This is not a spatial vector field.
config : dict
    A validated flat configuration from ``musselflow_ecogrammar_core``.

The main function returns transparent physical rates, scenario arrays,
objective scores, constraints, warnings, and one constraint-dominating fitness
for Galapagos.
"""

from __future__ import annotations

import math

import numpy as np

from musselflow_ecogrammar_core import resolved_lists


O2_MG_PER_ML = 1.42903
SECONDS_PER_DAY = 86400.0
LITRES_PER_CUBIC_METRE = 1000.0


def oxygen_saturation_mg_l(temperature_c, salinity_psu):
    """Garcia & Gordon (1992) oxygen solubility at one atmosphere.

    Returns mg/L.  The fitted equation's published envelope is approximately
    -2 to 40 C and salinity 0 to 42.
    """
    temperature_c = float(temperature_c)
    salinity_psu = float(salinity_psu)
    scaled_temperature = math.log(
        (298.15-temperature_c)/(273.15+temperature_c))
    a = (2.00907, 3.22014, 4.05010, 4.94457, -0.256847, 3.88767)
    b = (-0.00624523, -0.00737614, -0.0103410, -0.00817083)
    c = -4.88682e-7
    log_ml_l = sum(
        coefficient*scaled_temperature**power
        for power, coefficient in enumerate(a))
    log_ml_l += salinity_psu*sum(
        coefficient*scaled_temperature**power
        for power, coefficient in enumerate(b))
    log_ml_l += c*salinity_psu*salinity_psu
    return math.exp(log_ml_l)*O2_MG_PER_ML


def lognormal_moment(mean, coefficient_of_variation, exponent):
    """Return E[X**exponent] for a lognormal size distribution."""
    mean = float(mean)
    cv = float(coefficient_of_variation)
    exponent = float(exponent)
    if mean <= 0.0:
        return 0.0
    return mean**exponent * (
        1.0+cv*cv)**(0.5*exponent*(exponent-1.0))


def salinity_activity(config):
    """Editable trapezoidal activity prior; never a universal species law."""
    salinity = config["site.salinity_psu"]
    zero_low = config["species.salinity_zero_low_psu"]
    full_low = config["species.salinity_full_low_psu"]
    full_high = config["species.salinity_full_high_psu"]
    zero_high = config["species.salinity_zero_high_psu"]
    if salinity <= zero_low or salinity >= zero_high:
        return 0.0
    if salinity < full_low:
        return (salinity-zero_low)/max(full_low-zero_low, 1e-12)
    if salinity <= full_high:
        return 1.0
    return (zero_high-salinity)/max(zero_high-full_high, 1e-12)


def food_activity(chlorophyll_ug_l, config):
    """Smooth valve-activity transition around an editable low-food value."""
    threshold = config["species.low_food_threshold_ug_l"]
    width = max(config["species.low_food_transition_ug_l"], 1e-9)
    argument = np.clip(
        (np.asarray(chlorophyll_ug_l, dtype=float)-threshold)/width,
        -60.0, 60.0)
    return 1.0/(1.0+np.exp(-argument))


def particle_retention(diameter_um, config):
    """Editable sigmoid retention curve for one representative particle size."""
    diameter = np.asarray(diameter_um, dtype=float)
    d50 = config["species.retention_d50_um"]
    width = max(config["species.retention_slope_um"], 1e-12)
    argument = np.clip((diameter-d50)/width, -60.0, 60.0)
    return 1.0/(1.0+np.exp(-argument))


def oxygen_activity(config):
    """Reduced filtration response to boundary oxygen saturation."""
    saturation = oxygen_saturation_mg_l(
        config["site.temperature_c"], config["site.salinity_psu"])
    fraction = config["site.boundary_do_mg_l"]/max(saturation, 1e-12)
    zero = config["species.oxygen_zero_saturation_fraction"]
    full = config["species.oxygen_full_saturation_fraction"]
    return float(np.clip((fraction-zero)/max(full-zero, 1e-12), 0.0, 1.0))


def current_activity(speed_m_s, config):
    """Aggregation-aware current response translated from Nielsen & Vismann.

    Small groups decline between the two editable speed breakpoints. The
    configured effective aggregation size interpolates toward the protected
    group response. This remains a screening curve, not a universal law.
    """
    start = config["species.current_clearance_start_m_s"]
    zero = config["species.current_clearance_zero_m_s"]
    small_group = float(np.clip(
        (zero-float(speed_m_s))/max(zero-start, 1e-12), 0.0, 1.0))
    aggregation = config["stocking.effective_aggregation_size"]
    protected_at = config["species.current_protection_group_size"]
    protection = float(np.clip(
        (aggregation-3.0)/max(protected_at-3.0, 1e-12), 0.0, 1.0))
    return small_group+(1.0-small_group)*protection


def feeding_state(chlorophyll_ug_l, tsm_mg_l, speed_m_s, config):
    """Return transparent DEB-lite feeding multipliers for one local state."""
    small_fraction = config["food.small_particle_fraction"]
    small_retention = float(particle_retention(
        config["food.small_particle_diameter_um"], config))
    large_retention = float(particle_retention(
        config["food.large_particle_diameter_um"], config))
    phyto_retention = config["species.retention_efficiency"] * (
        small_fraction*small_retention +
        (1.0-small_fraction)*large_retention)

    detritus = config["food.detritus_fraction_of_organic"]
    particulate_retention = (
        config["species.particulate_retention_efficiency"] *
        ((1.0-detritus)+detritus*config["food.detritus_preference"]))

    low_food = float(food_activity(chlorophyll_ug_l, config))
    half_saturation = config["species.ingestion_half_saturation_ug_l"]
    saturation = half_saturation/max(
        half_saturation+float(chlorophyll_ug_l), 1e-12)
    oxygen = oxygen_activity(config)
    current = current_activity(speed_m_s, config)
    clearance_activity = low_food*saturation*oxygen*current

    threshold = config["species.pseudofaeces_tsm_threshold_mg_l"]
    transition = config["species.pseudofaeces_tsm_transition_mg_l"]
    rejection_gate = 1.0/(1.0+math.exp(-float(np.clip(
        (float(tsm_mg_l)-threshold)/max(transition, 1e-12),
        -60.0, 60.0))))
    pseudo_min = config["species.pseudofaeces_fraction"]
    pseudo_max = config["species.pseudofaeces_max_fraction"]
    pseudofaeces = pseudo_min+(pseudo_max-pseudo_min)*rejection_gate

    organic_fraction = config["site.particulate_organic_fraction"]
    quality_half = config["species.assimilation_quality_half_saturation"]
    reference = config["species.assimilation_reference_organic_fraction"]
    quality = (
        organic_fraction/max(organic_fraction+quality_half, 1e-12) /
        max(reference/max(reference+quality_half, 1e-12), 1e-12))
    quality = float(np.clip(quality, 0.0, 1.0))
    detritus_multiplier = (
        (1.0-detritus) +
        detritus*config["food.detritus_assimilation_multiplier"])
    high_threshold = config[
        "species.high_food_assimilation_threshold_ug_l"]
    high_decay = config["species.high_food_assimilation_decay_ug_l"]
    excess = max(float(chlorophyll_ug_l)-high_threshold, 0.0)
    high_food = math.exp(-excess/max(high_decay, 1e-12))
    assimilation = float(np.clip(
        config["species.assimilation_efficiency"]*quality *
        detritus_multiplier*high_food, 0.0, 1.0))
    return {
        "clearance_activity": clearance_activity,
        "low_food_activity": low_food,
        "ingestion_saturation": saturation,
        "oxygen_activity": oxygen,
        "current_activity": current,
        "phyto_retention": float(np.clip(phyto_retention, 0.0, 1.0)),
        "particulate_retention": float(np.clip(
            particulate_retention, 0.0, 1.0)),
        "pseudofaeces_fraction": float(np.clip(pseudofaeces, 0.0, 1.0)),
        "assimilation_efficiency": assimilation,
    }


def polygon_area(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or len(points) < 3 or points.shape[1] != 2:
        return 0.0
    x = points[:, 0]
    y = points[:, 1]
    return 0.5*abs(float(np.dot(x, np.roll(y, -1)) -
                         np.dot(y, np.roll(x, -1))))


def _normalise(vector, fallback=(1.0, 0.0)):
    vector = np.asarray(vector, dtype=float)
    length = float(np.linalg.norm(vector))
    if length <= 1e-15:
        return np.asarray(fallback, dtype=float), 0.0
    return vector/length, length


def _projected_geometry(obstacles, flow):
    """Return flow-aligned centres and projected rectangle dimensions."""
    centres = obstacles[:, :2]
    major = obstacles[:, 2]
    minor = obstacles[:, 3]
    yaw = obstacles[:, 4]
    major_axis = np.column_stack((np.cos(yaw), np.sin(yaw)))
    minor_axis = np.column_stack((-np.sin(yaw), np.cos(yaw)))
    cross = np.array([-flow[1], flow[0]])
    stream_length = (
        major*np.abs(major_axis@flow) + minor*np.abs(minor_axis@flow))
    cross_width = (
        major*np.abs(major_axis@cross) + minor*np.abs(minor_axis@cross))
    source_u = centres@flow
    source_v = centres@cross
    return source_u, source_v, stream_length, cross_width, cross


def _net_drag_outline(solidity, net_profile, multiplier):
    """Return effective outline-area drag coefficient.

    Lader et al.'s quoted knotless/knotted relations describe drag relative to
    solid twine area.  Multiplication by solidity converts to outline area for
    this screening implementation.
    """
    solidity = np.asarray(solidity, dtype=float)
    if "knotted" in str(net_profile).lower() and "knotless" not in str(
            net_profile).lower():
        coefficient_solid = 1.0+1.89*solidity+2.34*solidity*solidity
    else:
        coefficient_solid = 1.0+1.37*solidity+0.78*solidity*solidity
    return coefficient_solid*solidity*float(multiplier)


def _wake_field(source_u, source_v, source_width, thrust, targets,
                flow, cross, wake_spread, min_speed_ratio):
    """Evaluate a Gaussian Jensen-style wake at arbitrary plan targets."""
    targets = np.asarray(targets, dtype=float).reshape((-1, 2))
    if len(targets) == 0:
        return np.empty(0), np.empty(0)
    target_u = targets@flow
    target_v = targets@cross
    downstream = target_u[:, None]-source_u[None, :]
    lateral_distance = target_v[:, None]-source_v[None, :]
    active = downstream > 1e-9

    positive_x = np.maximum(downstream, 0.0)
    diameter = np.maximum(source_width[None, :], 1e-9)
    sigma = np.maximum(
        0.5*diameter+wake_spread*positive_x, 1e-9)
    lateral = np.exp(-0.5*(lateral_distance/sigma)**2)
    centre_deficit = (
        1.0-np.sqrt(np.maximum(1.0-thrust[None, :], 0.0)))
    decay = (1.0+2.0*wake_spread*positive_x/diameter)**-2
    pair_deficit = np.where(
        active, np.clip(centre_deficit*decay*lateral, 0.0, 0.95), 0.0)
    speed_ratio = np.clip(
        np.prod(1.0-pair_deficit, axis=1),
        min_speed_ratio, 1.0)

    # Direction change remains a visualization/screening proxy because this
    # reduced wake does not solve transverse momentum or continuity.
    side = np.sign(lateral_distance)
    lateral_ratio = np.sum(
        np.where(
            active,
            side*pair_deficit*np.minimum(0.30, 0.15*diameter/sigma),
            0.0),
        axis=1)
    deflection = np.degrees(np.arctan2(
        np.abs(lateral_ratio), np.maximum(speed_ratio, min_speed_ratio)))
    return speed_ratio, np.clip(deflection, 0.0, 90.0)


def _food_at_targets(source_u, source_v, source_width, removal_fraction,
                     targets, flow, cross, spread, recovery_lengths):
    """Return remaining food fraction from sequential source removals."""
    targets = np.asarray(targets, dtype=float).reshape((-1, 2))
    if len(targets) == 0:
        return np.empty(0)
    target_u = targets@flow
    target_v = targets@cross
    downstream = target_u[:, None]-source_u[None, :]
    lateral_distance = target_v[:, None]-source_v[None, :]
    active = downstream > 1e-9
    positive_x = np.maximum(downstream, 0.0)
    width = np.maximum(source_width[None, :], 1e-9)
    sigma = np.maximum(0.5*width+spread*positive_x, 1e-9)
    lateral = np.exp(-0.5*(lateral_distance/sigma)**2)
    recovery = np.exp(
        -positive_x/np.maximum(recovery_lengths*width, 1e-9))
    pair_removal = np.where(
        active,
        np.clip(removal_fraction[None, :]*lateral*recovery, 0.0, 0.98),
        0.0)
    return np.clip(np.prod(1.0-pair_removal, axis=1), 0.0, 1.0)


def _stocking(config, obstacle_count):
    """Resolve animal count, dry biomass and standing live wet biomass."""
    mean_dry_g = config["species.mean_dry_tissue_g"]
    mode = config["stocking.mode"].lower()
    animals = np.asarray(
        config["stocking.mussels_per_obstacle"], dtype=float)
    dry_override = np.asarray(
        config["stocking.dry_tissue_kg_per_obstacle"], dtype=float)
    if mode == "dry_biomass":
        if np.any(dry_override < 0.0):
            raise ValueError(
                "stocking.mode=dry_biomass requires non-negative "
                "stocking.dry_tissue_kg_per_obstacle values.")
        dry_kg = dry_override
        animals = dry_kg*1000.0/max(mean_dry_g, 1e-12)
    else:
        dry_kg = animals*mean_dry_g/1000.0
        use_override = dry_override >= 0.0
        if np.any(use_override):
            dry_kg = np.where(use_override, dry_override, dry_kg)
            animals = np.where(
                use_override,
                dry_kg*1000.0/max(mean_dry_g, 1e-12),
                animals)
    wet_kg = (
        animals*config["stocking.live_wet_g_per_individual"]/1000.0)
    if len(animals) != obstacle_count:
        raise ValueError("Resolved stocking list does not match obstacles.")
    return animals, dry_kg, wet_kg


def _population_rates(config, animal_count):
    """Return per-obstacle maximum clearance and respiration rates."""
    mean_dry = config["species.mean_dry_tissue_g"]
    size_cv = config["species.size_cv"]
    salinity_factor = salinity_activity(config)
    common_activity = config["species.activity_fraction"]*salinity_factor

    clearance_moment = lognormal_moment(
        mean_dry, size_cv, config["species.clearance_b"])
    clearance_temperature = config["species.clearance_q10"]**(
        (config["site.temperature_c"] -
         config["species.clearance_ref_temp_c"])/10.0)
    clearance_l_h = (
        animal_count*config["species.clearance_a_l_h"]*clearance_moment *
        clearance_temperature*common_activity)

    respiration_moment = lognormal_moment(
        mean_dry, size_cv, config["species.respiration_b"])
    respiration_temperature = config["species.respiration_q10"]**(
        (config["site.temperature_c"] -
         config["species.respiration_ref_temp_c"])/10.0)
    respiration_ml_h = (
        animal_count*config["species.respiration_a_ml_o2_h"] *
        respiration_moment*respiration_temperature*common_activity)
    return clearance_l_h, respiration_ml_h, salinity_factor


def _scenario(config, obstacles, probes, domain_polygon, flow_vector,
              animal_count, maximum_clearance_l_h,
              hydraulic_profile=None):
    """Evaluate one uniform-current scenario."""
    fallback = (1.0, 0.0)
    flow, speed = _normalise(flow_vector, fallback)
    source_u, source_v, stream_length, cross_width, cross = (
        _projected_geometry(obstacles, flow))
    domain_u = domain_polygon@flow
    domain_v = domain_polygon@cross
    domain_stream_length = max(float(np.ptp(domain_u)), 1e-9)
    domain_cross_width = max(float(np.ptp(domain_v)), 1e-9)

    porosity = np.asarray(config["structure.porosity"], dtype=float)
    solidity = 1.0-porosity
    outline_drag = _net_drag_outline(
        solidity, config["profile.net"],
        config["structure.drag_multiplier"])
    vertical_blockage = np.clip(
        obstacles[:, 5]/config["site.depth_m"], 0.0, 1.0)
    envelope_area = np.maximum(
        cross_width*obstacles[:, 5], 1e-12)
    frontal_area = envelope_area.copy()
    frontal_fill = np.ones(len(obstacles), dtype=float)
    if hydraulic_profile is not None:
        frontal_area = np.asarray(
            hydraulic_profile["frontal_area_m2"], dtype=float)
        frontal_fill = np.asarray(
            hydraulic_profile["frontal_fill"], dtype=float)
        if frontal_area.shape != (len(obstacles),):
            raise ValueError(
                "hydraulic frontal_area_m2 must match obstacle count.")
        if frontal_fill.shape != (len(obstacles),):
            raise ValueError(
                "hydraulic frontal_fill must match obstacle count.")
        if (not np.all(np.isfinite(frontal_area)) or
                not np.all(np.isfinite(frontal_fill))):
            raise ValueError("hydraulic profile contains non-finite values.")
        frontal_area = np.clip(frontal_area, 1e-12, envelope_area)
        frontal_fill = np.clip(frontal_fill, 0.01, 1.0)
    thrust = np.clip(
        1.0-np.exp(
            -outline_drag*vertical_blockage*frontal_fill),
        0.0, 0.95)

    centre_ratio, centre_deflection = _wake_field(
        source_u, source_v, cross_width, thrust, obstacles[:, :2],
        flow, cross, config["hydrodynamics.wake_spread"],
        config["hydrodynamics.min_speed_ratio"])
    probe_ratio, probe_deflection = _wake_field(
        source_u, source_v, cross_width, thrust, probes,
        flow, cross, config["hydrodynamics.wake_spread"],
        config["hydrodynamics.min_speed_ratio"])

    local_food = np.ones(len(obstacles), dtype=float)
    local_particulate = np.ones(len(obstacles), dtype=float)
    removal_fraction = np.zeros(len(obstacles), dtype=float)
    particulate_removal_fraction = np.zeros(len(obstacles), dtype=float)
    effective_cleared_m3_s = np.zeros(len(obstacles), dtype=float)
    captured_chlorophyll_mg_s = np.zeros(len(obstacles), dtype=float)
    captured_particulate_kg_day = np.zeros(len(obstacles), dtype=float)
    active_clearance_l_h_values = np.zeros(len(obstacles), dtype=float)
    low_food_values = np.zeros(len(obstacles), dtype=float)
    ingestion_saturation_values = np.zeros(len(obstacles), dtype=float)
    oxygen_activity_values = np.zeros(len(obstacles), dtype=float)
    current_activity_values = np.zeros(len(obstacles), dtype=float)
    phyto_retention_values = np.zeros(len(obstacles), dtype=float)
    particulate_retention_values = np.zeros(len(obstacles), dtype=float)
    pseudofaeces_values = np.zeros(len(obstacles), dtype=float)
    assimilation_values = np.zeros(len(obstacles), dtype=float)
    order = np.argsort(source_u, kind="stable")
    processed = []

    boundary_chlorophyll_mg_m3 = config["site.chlorophyll_ug_l"]
    tsm_g_m3 = config["site.tsm_mg_l"]
    for obstacle_index in order:
        if processed:
            previous = np.asarray(processed, dtype=int)
            target = obstacles[obstacle_index:obstacle_index+1, :2]
            incoming = _food_at_targets(
                source_u[previous], source_v[previous],
                cross_width[previous], removal_fraction[previous],
                target, flow, cross,
                config["hydrodynamics.food_plume_spread"],
                config["hydrodynamics.food_recovery_lengths"])[0]
            particulate_incoming = _food_at_targets(
                source_u[previous], source_v[previous],
                cross_width[previous],
                particulate_removal_fraction[previous],
                target, flow, cross,
                config["hydrodynamics.food_plume_spread"],
                config["hydrodynamics.food_recovery_lengths"])[0]
        else:
            incoming = 1.0
            particulate_incoming = 1.0
        local_food[obstacle_index] = incoming
        local_particulate[obstacle_index] = particulate_incoming

        local_chlorophyll = boundary_chlorophyll_mg_m3*incoming
        local_tsm = tsm_g_m3*particulate_incoming
        feeding = feeding_state(
            local_chlorophyll, local_tsm,
            speed*centre_ratio[obstacle_index], config)
        active_clearance_l_h = (
            maximum_clearance_l_h[obstacle_index] *
            feeding["clearance_activity"])
        active_clearance_l_h_values[obstacle_index] = active_clearance_l_h
        low_food_values[obstacle_index] = feeding["low_food_activity"]
        ingestion_saturation_values[obstacle_index] = feeding[
            "ingestion_saturation"]
        oxygen_activity_values[obstacle_index] = feeding["oxygen_activity"]
        current_activity_values[obstacle_index] = feeding["current_activity"]
        phyto_retention_values[obstacle_index] = feeding["phyto_retention"]
        particulate_retention_values[obstacle_index] = feeding[
            "particulate_retention"]
        pseudofaeces_values[obstacle_index] = feeding[
            "pseudofaeces_fraction"]
        assimilation_values[obstacle_index] = feeding[
            "assimilation_efficiency"]
        clearance_m3_s = (
            active_clearance_l_h /
            (LITRES_PER_CUBIC_METRE*3600.0))
        outline_area = frontal_area[obstacle_index]
        advective_flux_m3_s = (
            speed*centre_ratio[obstacle_index]*outline_area *
            porosity[obstacle_index])
        if advective_flux_m3_s > 1e-15 and clearance_m3_s > 0.0:
            removal = 1.0-math.exp(
                -clearance_m3_s*feeding["phyto_retention"] /
                advective_flux_m3_s)
            removal = min(max(removal, 0.0), 0.98)
            particulate_removal = 1.0-math.exp(
                -clearance_m3_s*feeding["particulate_retention"] /
                advective_flux_m3_s)
            particulate_removal = min(
                max(particulate_removal, 0.0), 0.98)
            effective = advective_flux_m3_s*removal
        else:
            removal = 0.0
            particulate_removal = 0.0
            effective = 0.0
        removal_fraction[obstacle_index] = removal
        particulate_removal_fraction[obstacle_index] = particulate_removal
        effective_cleared_m3_s[obstacle_index] = effective
        captured_chlorophyll_mg_s[obstacle_index] = (
            effective*local_chlorophyll)
        captured_particulate_kg_day[obstacle_index] = (
            advective_flux_m3_s*particulate_removal *
            (tsm_g_m3*particulate_incoming)*SECONDS_PER_DAY/1000.0)
        processed.append(obstacle_index)

    probe_food = _food_at_targets(
        source_u, source_v, cross_width, removal_fraction, probes,
        flow, cross, config["hydrodynamics.food_plume_spread"],
        config["hydrodynamics.food_recovery_lengths"])
    probe_particulate = _food_at_targets(
        source_u, source_v, cross_width, particulate_removal_fraction,
        probes, flow, cross, config["hydrodynamics.food_plume_spread"],
        config["hydrodynamics.food_recovery_lengths"])

    # Domain-scale mass cap: a reduced plume model must not capture more than
    # the chlorophyll or suspended mass advected through the domain section.
    domain_flux_m3_s = speed*domain_cross_width*config["site.depth_m"]
    chlorophyll_inflow_mg_s = (
        domain_flux_m3_s*boundary_chlorophyll_mg_m3)
    particulate_inflow_kg_day = (
        domain_flux_m3_s*tsm_g_m3*SECONDS_PER_DAY/1000.0)
    captured_chlorophyll_total = float(
        np.sum(captured_chlorophyll_mg_s))
    captured_particulate_total = float(
        np.sum(captured_particulate_kg_day))
    chlorophyll_scale = (
        1.0 if captured_chlorophyll_total <= 1e-30 else
        min(1.0, chlorophyll_inflow_mg_s/captured_chlorophyll_total))
    particulate_scale = (
        1.0 if captured_particulate_total <= 1e-30 else
        min(1.0, particulate_inflow_kg_day/captured_particulate_total))
    mass_balance_scale = min(chlorophyll_scale, particulate_scale)
    if chlorophyll_scale < 1.0:
        captured_chlorophyll_mg_s *= chlorophyll_scale
        effective_cleared_m3_s *= chlorophyll_scale
    if particulate_scale < 1.0:
        captured_particulate_kg_day *= particulate_scale

    organic_capture_kg_day = (
        captured_particulate_kg_day *
        config["site.particulate_organic_fraction"])
    pseudofaeces_organic_kg_day = (
        organic_capture_kg_day*pseudofaeces_values)
    ingested_organic_kg_day = (
        organic_capture_kg_day-pseudofaeces_organic_kg_day)
    assimilated_organic_kg_day = (
        ingested_organic_kg_day*assimilation_values)
    faeces_organic_kg_day = (
        ingested_organic_kg_day-assimilated_organic_kg_day)
    biodeposit_organic_kg_day = (
        pseudofaeces_organic_kg_day+faeces_organic_kg_day)
    organic_carbon_fraction = config["food.organic_carbon_fraction"]
    filtered_carbon_kg_day = (
        organic_capture_kg_day*organic_carbon_fraction)
    assimilated_carbon_kg_day = (
        assimilated_organic_kg_day*organic_carbon_fraction)
    phytoplankton_carbon_capture_kg_day = (
        captured_chlorophyll_mg_s*SECONDS_PER_DAY/1000.0 *
        config["food.carbon_to_chlorophyll_mg_c_per_mg_chl"] / 1000.0)

    # Deposition transport speed uses the deterministic obstacle-centre control
    # set, never user probes: diagnostic probe placement must not change farm
    # physics (deposition, oxygen, capture) or farm-wide objectives.
    mean_transport_speed = speed*float(np.mean(centre_ratio))
    settling_distance = (
        mean_transport_speed*config["site.depth_m"] /
        max(config["sediment.settling_velocity_m_s"], 1e-12))
    specified_deposition = config["sediment.in_domain_deposition_fraction"]
    if specified_deposition >= 0.0:
        deposition_fraction = min(specified_deposition, 1.0)
    elif mean_transport_speed <= 1e-15:
        deposition_fraction = 1.0
    else:
        deposition_fraction = 1.0-math.exp(
            -domain_stream_length/max(settling_distance, 1e-12))

    twine_reynolds = (
        speed*config["structure.twine_diameter_m"] /
        config["hydrodynamics.kinematic_viscosity_m2_s"])

    return {
        "speed_m_s": speed,
        "direction": flow,
        "cross": cross,
        "domain_stream_length_m": domain_stream_length,
        "domain_cross_width_m": domain_cross_width,
        "centre_speed_ratio": centre_ratio,
        "centre_speed_m_s": centre_ratio*speed,
        "centre_deflection_deg": centre_deflection,
        "local_food_fraction": local_food,
        "local_particulate_fraction": local_particulate,
        "removal_fraction": removal_fraction,
        "particulate_removal_fraction": particulate_removal_fraction,
        "probe_speed_ratio": probe_ratio,
        "probe_speed_m_s": probe_ratio*speed,
        "probe_food_fraction": probe_food,
        "probe_particulate_fraction": probe_particulate,
        "probe_deflection_deg": probe_deflection,
        "effective_cleared_m3_day": (
            effective_cleared_m3_s*SECONDS_PER_DAY),
        "chlorophyll_capture_g_day": (
            captured_chlorophyll_mg_s*SECONDS_PER_DAY/1000.0),
        "particulate_capture_kg_day": captured_particulate_kg_day,
        "organic_capture_kg_day": organic_capture_kg_day,
        "ingested_organic_kg_day": ingested_organic_kg_day,
        "assimilated_organic_kg_day": assimilated_organic_kg_day,
        "pseudofaeces_organic_kg_day": pseudofaeces_organic_kg_day,
        "faeces_organic_kg_day": faeces_organic_kg_day,
        "biodeposit_organic_kg_day": biodeposit_organic_kg_day,
        "filtered_carbon_kg_day": filtered_carbon_kg_day,
        "assimilated_carbon_kg_day": assimilated_carbon_kg_day,
        "phytoplankton_carbon_capture_kg_day":
            phytoplankton_carbon_capture_kg_day,
        "active_clearance_l_h": active_clearance_l_h_values,
        "low_food_activity": low_food_values,
        "ingestion_saturation": ingestion_saturation_values,
        "oxygen_activity": oxygen_activity_values,
        "current_activity": current_activity_values,
        "phyto_retention": phyto_retention_values,
        "particulate_retention": particulate_retention_values,
        "pseudofaeces_fraction": pseudofaeces_values,
        "assimilation_efficiency": assimilation_values,
        "deposition_fraction": deposition_fraction,
        "twine_reynolds": twine_reynolds,
        "mass_balance_scale": mass_balance_scale,
        "solidity": solidity,
        "outline_drag_coefficient": outline_drag,
        "thrust_proxy": thrust,
        "frontal_area_m2": frontal_area,
        "frontal_fill": frontal_fill,
    }


def _weighted_stack(scenarios, key, weights):
    values = np.stack([scenario[key] for scenario in scenarios], axis=0)
    return np.tensordot(weights, values, axes=(0, 0)), values


def _oxygen_sequence(config, scenarios, weights, domain_area,
                     respiration_kg_o2_day,
                     mortality_organic_kg_day):
    """Advance a transparent one-box oxygen/deposit screening balance."""
    depth = config["site.depth_m"]
    volume = domain_area*depth
    if volume <= 0.0:
        raise ValueError("Domain water volume must be positive.")
    saturation = oxygen_saturation_mg_l(
        config["site.temperature_c"], config["site.salinity_psu"])
    oxygen = config["site.initial_do_mg_l"]
    boundary_oxygen = config["site.boundary_do_mg_l"]
    stock = config["sediment.initial_organic_stock_kg"]
    decay_rate = config["sediment.decay_per_day"]
    resuspension_rate = config["sediment.resuspension_per_day"]
    total_stock_loss_rate = decay_rate+resuspension_rate
    oxygen_demand = config[
        "sediment.oxygen_demand_kg_o2_per_kg_organic"]
    durations = config["scenario.duration_h"]
    repeat_count = config["scenario.repeat_count"]
    series = [oxygen]
    deposit_series = [stock]
    minimum = oxygen
    scenario_end = []
    scenario_deposition_rate = []

    static_sink_kg_day = (
        respiration_kg_o2_day +
        config["site.pelagic_respiration_g_o2_m3_day"]*volume/1000.0 +
        config["site.background_sod_g_o2_m2_day"]*domain_area/1000.0)

    for repeat_index in range(repeat_count):
        for scenario_index, scenario in enumerate(scenarios):
            duration_days = durations[scenario_index]/24.0
            deposition_rate = (
                float(np.sum(scenario["biodeposit_organic_kg_day"])) *
                scenario["deposition_fraction"] +
                mortality_organic_kg_day *
                config["sediment.mortality_deposition_fraction"])
            scenario_deposition_rate.append(deposition_rate)

            if total_stock_loss_rate > 1e-15:
                decay_factor = math.exp(
                    -total_stock_loss_rate*duration_days)
                stock_end = (
                    stock*decay_factor +
                    deposition_rate/total_stock_loss_rate *
                    (1.0-decay_factor))
            else:
                stock_end = stock+deposition_rate*duration_days
            removed_mass = max(
                stock+deposition_rate*duration_days-stock_end, 0.0)
            decayed_mass = (
                removed_mass*decay_rate/total_stock_loss_rate
                if total_stock_loss_rate > 1e-15 else 0.0)
            deposit_oxygen_kg_day = (
                decayed_mass*oxygen_demand/max(duration_days, 1e-12))
            total_sink_kg_day = static_sink_kg_day+deposit_oxygen_kg_day
            sink_mg_l_day = total_sink_kg_day/volume*1000.0

            # Oxygen advection uses the obstacle-centre control set, never user
            # probes (see _scenario): probes are diagnostic only.
            mean_speed_ratio = float(
                np.mean(scenario["centre_speed_ratio"]))
            throughflow_m3_s = (
                scenario["speed_m_s"]*mean_speed_ratio *
                scenario["domain_cross_width_m"]*depth)
            advective_rate = (
                throughflow_m3_s/volume*SECONDS_PER_DAY *
                config["site.advective_exchange_efficiency"])
            vertical_rate = config["site.vertical_exchange_per_day"]
            reaeration_rate = config["site.reaeration_per_day"]
            total_exchange = advective_rate+vertical_rate+reaeration_rate
            production = config["site.primary_production_g_o2_m3_day"]

            forcing = (
                (advective_rate+vertical_rate)*boundary_oxygen +
                reaeration_rate*saturation+production-sink_mg_l_day)
            if total_exchange > 1e-15:
                equilibrium = forcing/total_exchange
                oxygen_end = (
                    equilibrium+(oxygen-equilibrium) *
                    math.exp(-total_exchange*duration_days))
            else:
                oxygen_end = oxygen+forcing*duration_days
            oxygen = max(0.0, oxygen_end)
            stock = max(0.0, stock_end)
            minimum = min(minimum, oxygen)
            series.append(oxygen)
            deposit_series.append(stock)
            if repeat_index == 0:
                scenario_end.append(oxygen)

    # The weighted daily deposition rate is a more useful ranking quantity than
    # the repeated temporal list.
    weighted_deposition = float(np.dot(
        weights,
        [float(np.sum(s["biodeposit_organic_kg_day"])) *
         s["deposition_fraction"] +
         mortality_organic_kg_day *
         config["sediment.mortality_deposition_fraction"]
         for s in scenarios]))
    return {
        "saturation_mg_l": saturation,
        "minimum_mg_l": minimum,
        "final_mg_l": oxygen,
        "series_mg_l": np.asarray(series, dtype=float),
        "scenario_end_mg_l": np.asarray(scenario_end, dtype=float),
        "organic_stock_series_kg": np.asarray(deposit_series, dtype=float),
        "final_organic_stock_kg": stock,
        "weighted_deposition_kg_day": weighted_deposition,
        "weighted_deposition_kg_m2_day": weighted_deposition/domain_area,
        "static_sink_kg_o2_day": static_sink_kg_day,
    }


def _coefficient_of_variation_score(values):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return 0.0
    mean = float(np.mean(values))
    if mean <= 1e-15:
        return 0.0
    return float(np.clip(1.0-np.std(values)/mean, 0.0, 1.0))


def evaluate_layout(obstacles, domain_polygon, probes, flow_vectors, config,
                    boundary_score=1.0, collision_score=1.0,
                    hydraulic_profiles=None):
    """Evaluate a layout and return a deterministic screening result."""
    obstacles = np.asarray(obstacles, dtype=float)
    domain_polygon = np.asarray(domain_polygon, dtype=float)
    probes = np.asarray(probes, dtype=float)
    flow_vectors = np.asarray(flow_vectors, dtype=float)

    if obstacles.ndim != 2 or obstacles.shape[1] != 6 or len(obstacles) == 0:
        raise ValueError("obstacles must be a non-empty (N,6) array.")
    if not np.all(np.isfinite(obstacles)):
        raise ValueError("obstacles contain non-finite values.")
    if np.any(obstacles[:, 2:4] <= 0.0) or np.any(obstacles[:, 5] <= 0.0):
        raise ValueError("Obstacle major, minor and height dimensions must be > 0.")
    if domain_polygon.ndim != 2 or domain_polygon.shape[1] != 2:
        raise ValueError("domain_polygon must be a (K,2) array.")
    domain_area = polygon_area(domain_polygon)
    if domain_area <= 0.0:
        raise ValueError("domain_polygon must enclose positive area.")
    if probes.size == 0:
        probes = obstacles[:, :2].copy()
        probe_source = "obstacle_centres"
    else:
        probes = probes.reshape((-1, 2))
        probe_source = "user_probes"
    if not np.all(np.isfinite(probes)):
        raise ValueError("probes contain non-finite values.")
    if flow_vectors.size == 0:
        raise ValueError("At least one flow vector in m/s is required.")
    flow_vectors = flow_vectors.reshape((-1, 2))
    if not np.all(np.isfinite(flow_vectors)):
        raise ValueError("flow_vectors contain non-finite values.")
    if hydraulic_profiles is None:
        hydraulic_profiles = [None]*len(flow_vectors)
    elif len(hydraulic_profiles) != len(flow_vectors):
        raise ValueError(
            "hydraulic_profiles must contain one profile per flow vector.")

    config = resolved_lists(config, len(obstacles), len(flow_vectors))
    weights = np.asarray(config["scenario.weights"], dtype=float)
    animal_count, dry_biomass_kg, wet_biomass_kg = _stocking(
        config, len(obstacles))
    maximum_clearance_l_h, respiration_ml_h, salinity_factor = (
        _population_rates(config, animal_count))

    scenarios = [
        _scenario(
            config, obstacles, probes, domain_polygon, vector,
            animal_count, maximum_clearance_l_h,
            hydraulic_profile=hydraulic_profiles[scenario_index])
        for scenario_index, vector in enumerate(flow_vectors)
    ]

    probe_speed_m_s, scenario_probe_speed_m_s = _weighted_stack(
        scenarios, "probe_speed_m_s", weights)
    probe_speed_ratio, scenario_probe_speed_ratio = _weighted_stack(
        scenarios, "probe_speed_ratio", weights)
    probe_food, scenario_probe_food = _weighted_stack(
        scenarios, "probe_food_fraction", weights)
    probe_deflection, scenario_probe_deflection = _weighted_stack(
        scenarios, "probe_deflection_deg", weights)
    obstacle_speed_m_s, scenario_obstacle_speed_m_s = _weighted_stack(
        scenarios, "centre_speed_m_s", weights)
    obstacle_speed_ratio, scenario_obstacle_speed_ratio = _weighted_stack(
        scenarios, "centre_speed_ratio", weights)
    obstacle_food, scenario_obstacle_food = _weighted_stack(
        scenarios, "local_food_fraction", weights)
    obstacle_removal, scenario_obstacle_removal = _weighted_stack(
        scenarios, "removal_fraction", weights)

    rate_keys = (
        "effective_cleared_m3_day",
        "chlorophyll_capture_g_day",
        "particulate_capture_kg_day",
        "organic_capture_kg_day",
        "ingested_organic_kg_day",
        "assimilated_organic_kg_day",
        "pseudofaeces_organic_kg_day",
        "faeces_organic_kg_day",
        "biodeposit_organic_kg_day",
        "filtered_carbon_kg_day",
        "assimilated_carbon_kg_day",
        "phytoplankton_carbon_capture_kg_day",
        "active_clearance_l_h",
        "low_food_activity",
        "ingestion_saturation",
        "oxygen_activity",
        "current_activity",
        "phyto_retention",
        "particulate_retention",
        "pseudofaeces_fraction",
        "assimilation_efficiency",
    )
    weighted_rates = {}
    scenario_rates = {}
    for key in rate_keys:
        weighted, stacked = _weighted_stack(scenarios, key, weights)
        weighted_rates[key] = weighted
        scenario_rates[key] = stacked

    respiration_kg_o2_day_by_obstacle = (
        respiration_ml_h*O2_MG_PER_ML*24.0/1.0e6)
    respiration_kg_o2_day = float(
        np.sum(respiration_kg_o2_day_by_obstacle))
    ammonia_kg_n_day_by_obstacle = (
        dry_biomass_kg*1000.0 *
        config["species.ammonia_mg_n_g_dw_h"]*24.0/1.0e6)
    ammonia_kg_n_day = float(np.sum(ammonia_kg_n_day_by_obstacle))
    mortality_organic_kg_day = (
        float(np.sum(dry_biomass_kg)) *
        config["stocking.annual_mortality_fraction"]/365.0)

    oxygen = _oxygen_sequence(
        config, scenarios, weights, domain_area, respiration_kg_o2_day,
        mortality_organic_kg_day)

    annual_survival = 1.0-config["stocking.annual_mortality_fraction"]
    harvested_wet_t_year = (
        float(np.sum(wet_biomass_kg))/1000.0 *
        config["harvest.fraction_per_year"]*annual_survival *
        config["harvest.turnovers_per_year"])
    harvest_n_kg_year = (
        harvested_wet_t_year*config["harvest.n_kg_per_t_wet"])
    harvest_p_kg_year = (
        harvested_wet_t_year*config["harvest.p_kg_per_t_wet"])

    capture_g_day = float(np.sum(
        weighted_rates["chlorophyll_capture_g_day"]))
    particulate_kg_day = float(np.sum(
        weighted_rates["particulate_capture_kg_day"]))
    cleared_m3_day = float(np.sum(
        weighted_rates["effective_cleared_m3_day"]))
    biodeposit_kg_day = float(np.sum(
        weighted_rates["biodeposit_organic_kg_day"]))
    assimilated_kg_day = float(np.sum(
        weighted_rates["assimilated_organic_kg_day"]))
    ingested_kg_day = float(np.sum(
        weighted_rates["ingested_organic_kg_day"]))
    pseudofaeces_kg_day = float(np.sum(
        weighted_rates["pseudofaeces_organic_kg_day"]))
    faeces_kg_day = float(np.sum(
        weighted_rates["faeces_organic_kg_day"]))
    filtered_carbon_kg_day = float(np.sum(
        weighted_rates["filtered_carbon_kg_day"]))
    assimilated_carbon_kg_day = float(np.sum(
        weighted_rates["assimilated_carbon_kg_day"]))
    phytoplankton_carbon_kg_day = float(np.sum(
        weighted_rates["phytoplankton_carbon_capture_kg_day"]))
    assimilated_energy_kj_day_by_obstacle = (
        weighted_rates["assimilated_organic_kg_day"]*1000.0 *
        config["food.organic_energy_kj_g"])
    respiration_energy_kj_day_by_obstacle = (
        respiration_ml_h*24.0 *
        config["species.oxycalorific_kj_per_ml_o2"])
    scope_for_growth_kj_day_by_obstacle = (
        assimilated_energy_kj_day_by_obstacle-
        respiration_energy_kj_day_by_obstacle)
    potential_growth_g_dw_day_by_obstacle = (
        np.maximum(scope_for_growth_kj_day_by_obstacle, 0.0) /
        config["species.tissue_energy_kj_g_dw"])

    maximum_clearance_m3_day = (
        float(np.sum(maximum_clearance_l_h)) /
        LITRES_PER_CUBIC_METRE*24.0)
    theoretical_capture_g_day = (
        maximum_clearance_m3_day *
        config["species.retention_efficiency"] *
        config["site.chlorophyll_ug_l"]/1000.0)
    target_capture = config[
        "objective.target_chlorophyll_capture_g_day"]
    if target_capture <= 0.0:
        target_capture = max(theoretical_capture_g_day, 1e-12)
    extraction_score = float(np.clip(
        capture_g_day/target_capture, 0.0, 1.0))
    # Flushing and uniformity are farm-wide objectives and must be independent
    # of diagnostic probe placement, so they read the obstacle-centre control
    # set, not user probes. (With no user probes the two are identical.)
    flushing_score = float(np.clip(
        np.mean(obstacle_speed_ratio), 0.0, 1.0))
    food_delivery_score = float(np.clip(
        np.mean(obstacle_food), 0.0, 1.0))
    uniformity_score = 0.5*(
        _coefficient_of_variation_score(obstacle_food) +
        _coefficient_of_variation_score(obstacle_speed_ratio))
    oxygen_reference = max(
        config["site.initial_do_mg_l"],
        config["site.boundary_do_mg_l"], 1e-12)
    oxygen_score = float(np.clip(
        oxygen["minimum_mg_l"]/oxygen_reference, 0.0, 1.0))
    deposition_target = max(
        config["objective.target_biodeposition_kg_m2_day"], 1e-12)
    low_deposition_score = float(np.clip(
        1.0-oxygen["weighted_deposition_kg_m2_day"]/deposition_target,
        0.0, 1.0))

    objectives = {
        "extraction": extraction_score,
        "flushing": flushing_score,
        "food_delivery": food_delivery_score,
        "uniformity": uniformity_score,
        "oxygen": oxygen_score,
        "low_deposition": low_deposition_score,
    }
    objective_weights = {
        name: config["weight."+name] for name in objectives}
    total_weight = sum(objective_weights.values())
    raw_objective = sum(
        objective_weights[name]*value
        for name, value in objectives.items())/total_weight

    speed_values = np.linalg.norm(flow_vectors, axis=1)
    mean_dry = config["species.mean_dry_tissue_g"]
    solidities = 1.0-np.asarray(config["structure.porosity"])
    envelope_violations = []
    envelope_soft_score = 1.0
    maximum_valid_speed = config["species.valid_flow_max_m_s"]
    if np.any(speed_values > maximum_valid_speed):
        envelope_violations.append("current_speed")
        envelope_soft_score = min(
            envelope_soft_score,
            maximum_valid_speed/max(float(np.max(speed_values)), 1e-12))
    if mean_dry < 0.011:
        envelope_violations.append("adult_dry_tissue_mass")
        envelope_soft_score = min(
            envelope_soft_score, mean_dry/0.011)
    elif mean_dry > 1.361:
        envelope_violations.append("adult_dry_tissue_mass")
        envelope_soft_score = min(
            envelope_soft_score, 1.361/mean_dry)
    if np.any((solidities < 0.2) | (solidities > 0.8)):
        envelope_violations.append("net_solidity")
        lower_score = np.where(
            solidities < 0.2, solidities/0.2, 1.0)
        upper_score = np.where(
            solidities > 0.8, 0.8/np.maximum(solidities, 1e-12), 1.0)
        envelope_soft_score = min(
            envelope_soft_score,
            float(np.min(np.minimum(lower_score, upper_score))))
    temperature = config["site.temperature_c"]
    salinity = config["site.salinity_psu"]
    if temperature < -2.0 or temperature > 40.0 or salinity > 42.0:
        envelope_violations.append("oxygen_solubility")
        envelope_soft_score = min(envelope_soft_score, 0.5)

    constraint_margins = {
        "boundary": float(boundary_score)-0.999,
        "collision": float(collision_score)-0.999,
    }
    soft_scores = [
        float(np.clip(boundary_score, 0.0, 1.0)),
        float(np.clip(collision_score, 0.0, 1.0)),
    ]
    if (envelope_violations and
            not config["validation.allow_extrapolation"]):
        constraint_margins["validation_envelope"] = -float(
            len(envelope_violations))
        soft_scores.append(float(np.clip(
            envelope_soft_score, 0.0, 1.0)))
    minimum_do = config["constraint.min_do_mg_l"]
    if minimum_do >= 0.0:
        margin = oxygen["minimum_mg_l"]-minimum_do
        constraint_margins["minimum_do_mg_l"] = margin
        soft_scores.append(float(np.clip(
            oxygen["minimum_mg_l"]/max(minimum_do, 1e-12), 0.0, 1.0)))
    minimum_speed = config["constraint.min_probe_speed_m_s"]
    if minimum_speed >= 0.0:
        observed_minimum_speed = float(np.min(scenario_probe_speed_m_s))
        constraint_margins["minimum_probe_speed_m_s"] = (
            observed_minimum_speed-minimum_speed)
        soft_scores.append(float(np.clip(
            observed_minimum_speed/max(minimum_speed, 1e-12), 0.0, 1.0)))
    maximum_deposition = config[
        "constraint.max_biodeposition_kg_m2_day"]
    if maximum_deposition >= 0.0:
        observed_deposition = oxygen["weighted_deposition_kg_m2_day"]
        constraint_margins["maximum_biodeposition_kg_m2_day"] = (
            maximum_deposition-observed_deposition)
        soft_scores.append(float(np.clip(
            maximum_deposition/max(observed_deposition, 1e-12), 0.0, 1.0)))

    feasible = all(value >= 0.0 for value in constraint_margins.values())
    soft_feasibility = float(np.prod(soft_scores))
    # Constraint-domination encoding: every feasible design scores in [0.25,1]
    # and every infeasible design in [0,0.25].  This lets Galapagos approach a
    # feasible region without allowing a high ecological score to outweigh a
    # configured hard limit.
    fitness = (
        0.25+0.75*raw_objective if feasible
        else 0.25*soft_feasibility)

    warnings = []
    if np.any(speed_values > config["species.valid_flow_max_m_s"]):
        warnings.append(
            "One or more currents exceed species.valid_flow_max_m_s; "
            "clearance behavior is outside the cited current study envelope.")
    if np.any(speed_values == 0.0):
        warnings.append(
            "A zero-speed scenario has no intrinsic direction; +X was used "
            "only to define coordinates and its advective capture is zero.")
    if mean_dry < 0.011 or mean_dry > 1.361:
        warnings.append(
            "Mean dry tissue mass is outside the cited adult clearance "
            "allometry envelope (0.011-1.361 g).")
    if np.any((solidities < 0.2) | (solidities > 0.8)):
        warnings.append(
            "Net solidity is outside the 0.2-0.8 range associated with the "
            "selected empirical net-drag prior.")
    if envelope_violations:
        warnings.append(
            "Outside declared use envelope: %s.%s"
            % (", ".join(envelope_violations),
               " Extrapolation was allowed by grammar."
               if config["validation.allow_extrapolation"]
               else " Candidate is constraint-infeasible."))
    if any(scenario["mass_balance_scale"] < 0.999999 for scenario in scenarios):
        warnings.append(
            "A domain inflow mass cap was active; the reduced plume model "
            "otherwise attempted to capture more material than entered.")
    if probe_source == "obstacle_centres":
        warnings.append(
            "No user probes were supplied: obstacle centres are the explicit "
            "fallback probe set.")
    if not feasible:
        warnings.append(
            "At least one configured hard constraint is violated.")
    if config["profile.species"].lower() != "mytilus_edulis_screening":
        warnings.append(
            "The species profile name is custom; verify every biological "
            "coefficient and evidence envelope.")
    warnings.append(
        "DEB-lite feeding and potential growth are transparent screening "
        "proxies, not a calibrated DEB state model or harvest forecast.")

    calibration_keys = (
        "validation.geometry_calibrated",
        "validation.hydrodynamics_calibrated",
        "validation.biology_calibrated",
        "validation.oxygen_calibrated",
        "validation.sediment_calibrated",
    )
    calibration_fraction = sum(
        1.0 for key in calibration_keys if config[key])/len(calibration_keys)
    model_status = (
        "CALIBRATED_WITHIN_DECLARED_ENVELOPE"
        if calibration_fraction == 1.0
        else "UNVALIDATED_SCREENING")

    return {
        "fitness": float(np.clip(fitness, 0.0, 1.0)),
        "raw_objective": float(np.clip(raw_objective, 0.0, 1.0)),
        "feasible": bool(feasible),
        "soft_feasibility": soft_feasibility,
        "constraint_margins": constraint_margins,
        "objectives": objectives,
        "objective_weights": objective_weights,
        "probe_source": probe_source,
        "probe_points": probes,
        "probe_speed_m_s": probe_speed_m_s,
        "probe_speed_ratio": probe_speed_ratio,
        "probe_food_fraction": probe_food,
        "probe_deflection_deg": probe_deflection,
        "scenario_probe_speed_m_s": scenario_probe_speed_m_s,
        "scenario_probe_speed_ratio": scenario_probe_speed_ratio,
        "scenario_probe_food_fraction": scenario_probe_food,
        "scenario_probe_deflection_deg": scenario_probe_deflection,
        "obstacle_speed_m_s": obstacle_speed_m_s,
        "obstacle_speed_ratio": obstacle_speed_ratio,
        "obstacle_food_fraction": obstacle_food,
        "obstacle_removal_fraction": obstacle_removal,
        "scenario_obstacle_speed_m_s": scenario_obstacle_speed_m_s,
        "scenario_obstacle_speed_ratio": scenario_obstacle_speed_ratio,
        "scenario_obstacle_food_fraction": scenario_obstacle_food,
        "scenario_obstacle_removal_fraction": scenario_obstacle_removal,
        "effective_cleared_m3_day_by_obstacle": weighted_rates[
            "effective_cleared_m3_day"],
        "chlorophyll_capture_g_day_by_obstacle": weighted_rates[
            "chlorophyll_capture_g_day"],
        "particulate_capture_kg_day_by_obstacle": weighted_rates[
            "particulate_capture_kg_day"],
        "ingested_organic_kg_day_by_obstacle": weighted_rates[
            "ingested_organic_kg_day"],
        "assimilated_organic_kg_day_by_obstacle": weighted_rates[
            "assimilated_organic_kg_day"],
        "pseudofaeces_organic_kg_day_by_obstacle": weighted_rates[
            "pseudofaeces_organic_kg_day"],
        "faeces_organic_kg_day_by_obstacle": weighted_rates[
            "faeces_organic_kg_day"],
        "biodeposit_organic_kg_day_by_obstacle": weighted_rates[
            "biodeposit_organic_kg_day"],
        "active_clearance_l_h_by_obstacle": weighted_rates[
            "active_clearance_l_h"],
        "phyto_retention_by_obstacle": weighted_rates[
            "phyto_retention"],
        "particulate_retention_by_obstacle": weighted_rates[
            "particulate_retention"],
        "pseudofaeces_fraction_by_obstacle": weighted_rates[
            "pseudofaeces_fraction"],
        "assimilation_efficiency_by_obstacle": weighted_rates[
            "assimilation_efficiency"],
        "oxygen_activity_by_obstacle": weighted_rates[
            "oxygen_activity"],
        "current_activity_by_obstacle": weighted_rates[
            "current_activity"],
        "scenario_rates": scenario_rates,
        "effective_cleared_m3_day": cleared_m3_day,
        "chlorophyll_capture_g_day": capture_g_day,
        "particulate_capture_kg_day": particulate_kg_day,
        "ingested_organic_kg_day": ingested_kg_day,
        "assimilated_organic_kg_day": assimilated_kg_day,
        "pseudofaeces_organic_kg_day": pseudofaeces_kg_day,
        "faeces_organic_kg_day": faeces_kg_day,
        "biodeposit_organic_kg_day": biodeposit_kg_day,
        "filtered_carbon_kg_day": filtered_carbon_kg_day,
        "assimilated_carbon_kg_day": assimilated_carbon_kg_day,
        "phytoplankton_carbon_capture_kg_day":
            phytoplankton_carbon_kg_day,
        "assimilated_energy_kj_day": float(np.sum(
            assimilated_energy_kj_day_by_obstacle)),
        "respiration_energy_kj_day": float(np.sum(
            respiration_energy_kj_day_by_obstacle)),
        "scope_for_growth_kj_day": float(np.sum(
            scope_for_growth_kj_day_by_obstacle)),
        "potential_growth_g_dw_day": float(np.sum(
            potential_growth_g_dw_day_by_obstacle)),
        "scope_for_growth_kj_day_by_obstacle":
            scope_for_growth_kj_day_by_obstacle,
        "potential_growth_g_dw_day_by_obstacle":
            potential_growth_g_dw_day_by_obstacle,
        "mussel_respiration_kg_o2_day": respiration_kg_o2_day,
        "mussel_respiration_kg_o2_day_by_obstacle":
            respiration_kg_o2_day_by_obstacle,
        "ammonia_excretion_kg_n_day": ammonia_kg_n_day,
        "ammonia_excretion_kg_n_day_by_obstacle":
            ammonia_kg_n_day_by_obstacle,
        "mortality_organic_kg_day": mortality_organic_kg_day,
        "oxygen": oxygen,
        "standing_live_wet_biomass_kg": float(np.sum(wet_biomass_kg)),
        "standing_dry_tissue_biomass_kg": float(np.sum(dry_biomass_kg)),
        "harvested_wet_t_year": harvested_wet_t_year,
        "harvest_n_kg_year": harvest_n_kg_year,
        "harvest_p_kg_year": harvest_p_kg_year,
        "theoretical_chlorophyll_capture_g_day":
            theoretical_capture_g_day,
        "scenario_weights": weights,
        "scenario_speeds_m_s": speed_values,
        "scenario_durations_h": np.asarray(
            config["scenario.duration_h"], dtype=float),
        "scenario_deposition_fraction": np.asarray(
            [scenario["deposition_fraction"] for scenario in scenarios]),
        "scenario_twine_reynolds": np.asarray(
            [scenario["twine_reynolds"] for scenario in scenarios]),
        "scenario_mass_balance_scale": np.asarray(
            [scenario["mass_balance_scale"] for scenario in scenarios]),
        "scenario_frontal_area_m2": np.asarray(
            [scenario["frontal_area_m2"] for scenario in scenarios]),
        "scenario_frontal_fill": np.asarray(
            [scenario["frontal_fill"] for scenario in scenarios]),
        "salinity_activity": salinity_factor,
        "calibration_fraction": calibration_fraction,
        "model_status": model_status,
        "warnings": warnings,
    }


def evaluate_ensemble(obstacles, domain_polygon, probes, flow_vectors,
                      state_configs, boundary_score=1.0, collision_score=1.0):
    """Evaluate an unordered ensemble of independent environmental states.

    ``state_configs`` is a list of dicts, one per state, each with a validated
    flat ``config`` (its own boundary chemistry already applied), a
    ``flow_vector_index`` into ``flow_vectors``, and an occurrence
    ``probability``.  Every state is evaluated independently with
    :func:`evaluate_layout`; farm objectives and screening rates are combined by
    probability, while hard constraints and the minimum dissolved oxygen take
    the worst state.  The aggregate is deterministic and independent of state
    order.  Each state's full result is returned as well -- the natural seam for
    a future residual surrogate to correct one condition at a time.
    """
    flow_vectors = np.asarray(flow_vectors, dtype=float).reshape((-1, 2))
    if not state_configs:
        raise ValueError("evaluate_ensemble requires at least one state.")

    per_state = []
    probabilities = []
    for state in state_configs:
        index = int(state["flow_vector_index"])
        if index < 0 or index >= len(flow_vectors):
            raise ValueError(
                "Ensemble state flow_vector_index %d is out of range." % index)
        result = evaluate_layout(
            obstacles, domain_polygon, probes,
            flow_vectors[index:index+1], state["config"],
            boundary_score=boundary_score, collision_score=collision_score)
        per_state.append((state, result))
        probabilities.append(float(state["probability"]))

    probability_sum = sum(probabilities)
    if probability_sum <= 0.0:
        raise ValueError("Ensemble probabilities must have a positive sum.")
    probabilities = [value/probability_sum for value in probabilities]

    objective_names = list(per_state[0][1]["objectives"])
    objectives = {
        name: sum(
            weight*float(result["objectives"][name])
            for weight, (_, result) in zip(probabilities, per_state))
        for name in objective_names}
    objective_weights = dict(per_state[0][1]["objective_weights"])
    total_weight = sum(objective_weights.values())
    raw_objective = sum(
        objective_weights[name]*objectives[name]
        for name in objective_names)/total_weight

    # Hard constraints: the worst (minimum) margin across the states that define
    # each constraint.  An envelope constraint that only some states trigger
    # therefore still fails the whole ensemble.
    constraint_margins = {}
    for _, result in per_state:
        for name, margin in result["constraint_margins"].items():
            margin = float(margin)
            if name not in constraint_margins or margin < constraint_margins[
                    name]:
                constraint_margins[name] = margin
    feasible = all(margin >= 0.0 for margin in constraint_margins.values())

    soft_feasibility = min(
        float(result["soft_feasibility"]) for _, result in per_state)
    fitness = (
        0.25+0.75*raw_objective if feasible
        else 0.25*soft_feasibility)

    def probability_weighted(getter):
        return sum(
            weight*float(getter(result))
            for weight, (_, result) in zip(probabilities, per_state))

    minimum_do = min(
        float(result["oxygen"]["minimum_mg_l"]) for _, result in per_state)
    status = (
        "OUT_OF_ENVELOPE"
        if constraint_margins.get("validation_envelope", 0.0) < 0.0
        else "INFEASIBLE" if not feasible
        else "UNVALIDATED_SCREENING")

    return {
        "mode": "ensemble",
        "status": status,
        "model_status": per_state[0][1]["model_status"],
        "fitness": float(np.clip(fitness, 0.0, 1.0)),
        "raw_objective": float(np.clip(raw_objective, 0.0, 1.0)),
        "feasible": bool(feasible),
        "soft_feasibility": soft_feasibility,
        "objectives": objectives,
        "objective_weights": objective_weights,
        "constraint_margins": constraint_margins,
        "probabilities": probabilities,
        "effective_cleared_m3_day": probability_weighted(
            lambda result: result["effective_cleared_m3_day"]),
        "chlorophyll_capture_g_day": probability_weighted(
            lambda result: result["chlorophyll_capture_g_day"]),
        "particulate_capture_kg_day": probability_weighted(
            lambda result: result["particulate_capture_kg_day"]),
        "biodeposit_organic_kg_day": probability_weighted(
            lambda result: result["biodeposit_organic_kg_day"]),
        "minimum_do_mg_l": minimum_do,
        "weighted_deposition_kg_m2_day": probability_weighted(
            lambda result: result["oxygen"]["weighted_deposition_kg_m2_day"]),
        "states": [
            {
                "id": state["id"],
                "flow_vector_index": int(state["flow_vector_index"]),
                "probability": probability,
                "result": result,
            }
            for probability, (state, result) in zip(
                probabilities, per_state)],
        "warnings": sorted({
            warning
            for _, result in per_state
            for warning in result["warnings"]}),
    }
