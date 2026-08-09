"""MusselFlow ecological grammar book — the single authoritative catalog.

This module is the "grammar book": one place that describes *every* factor the
screening model understands.  For each grammar key it records

* the default screening value and its unit;
* a plain-language summary;
* the governing equation and scientific validity envelope (where one exists);
* the literature citation and DOI (where one exists);
* whether the strict ecological case exposes the field, derives it, or keeps a
  compatibility default.

Design rules that keep this honest:

* It is **pure data with no Rhino, NumPy, or solver dependency**, so it can be
  imported anywhere and rendered into a Grasshopper Panel or a document.
* It **changes no runtime behavior**.  Nothing computes from it yet; the
  numerical core still reads ``musselflow_ecogrammar_core.DEFAULTS``.
* ``test_musselflow_grammar_book.py`` pins every entry to that ``DEFAULTS``
  dictionary in both directions, so the book cannot drift out of sync with the
  code that actually runs, and a new grammar key cannot be added to the core
  without also being documented here.
* Values are labelled **screening priors**, not universal constants.  Where a
  value is a model tuning parameter with no single literature source, its
  citation is ``None`` and it is described as such.  No DOI is ever invented.

The catalog and ``musselflow_grammar.json`` together define the ecological
interface. A future evaluator or optimizer can read units and envelopes from
one place instead of hard-coding them.
"""

from __future__ import annotations


# --- Ecological-case coverage status ---------------------------------------
# How musselflow_grammar.json + musselflow_case_core.py treat each field.
EXPOSED = "exposed"       # user-settable directly in the case JSON
BOUNDARY = "boundary"     # user-settable via a forcing-step boundary block
DERIVED = "derived"       # computed by the runtime; never user-set
DEFAULTED = "defaulted"   # compatibility-only metadata/default
DROPPED = "dropped"       # capability not represented by the current case

CASE_STATUSES = {EXPOSED, BOUNDARY, DERIVED, DEFAULTED, DROPPED}

# Fields the ecological case cannot currently drive.
CASE_MISSING = {DEFAULTED, DROPPED}


# --- Citation registry ------------------------------------------------------
# Shared provenance so keys that come from the same study reference one record.
# Text and URLs are taken from the project's own PROVENANCE table and
# MUSSELFLOW_PHYSICS_INTELLIGENCE_ARCHITECTURE.md; none are fabricated.
CITATIONS = {
    "clearance": {
        "citation": "Moehlenberg & Riisgaard (1979), adult Mytilus edulis; "
                    "synthesis Riisgaard et al. (2025)",
        "equation": "clearance_L_h = 7.45 * dry_tissue_g ** 0.66",
        "envelope": "adult dry tissue 0.011-1.361 g; approximately 10-13 C, "
                    "30 PSU",
        "url": "https://doi.org/10.1242/bio.062024",
    },
    "respiration": {
        "citation": "Hamburger et al. (1983), adult Mytilus edulis; "
                    "synthesis Riisgaard et al. (2025)",
        "equation": "respiration_mL_O2_h = 0.475 * dry_tissue_g ** 0.663",
        "envelope": "adult dry tissue; approximately 10-13 C, 30 PSU",
        "url": "https://doi.org/10.1242/bio.062024",
    },
    "current": {
        "citation": "Nielsen & Vismann (2014)",
        "equation": "No universal multiplier; speeds above the envelope are "
                    "flagged, not silently penalised",
        "envelope": "0.05-1.4 m/s; aggregation strongly changed the response "
                    "(groups of 20 held clearance where groups of 3 did not)",
        "url": "https://doi.org/10.2983/035.033.0214",
    },
    "retention": {
        "citation": "Strohmeier et al. (2012)",
        "equation": "Direct editable retention efficiency (no size resolution)",
        "envelope": "Natural seston; strong particle-size and temporal "
                    "variability (peak retention near 30-35 um)",
        "url": "https://doi.org/10.1016/j.jembe.2011.11.006",
    },
    "net_drag": {
        "citation": "Lader et al. (2009), aquaculture net drag",
        "equation": "Cd_solid = 1 + 1.37*S + 0.78*S**2 for a knotless net, "
                    "with S = 1 - porosity",
        "envelope": "Net-panel empirical relation; depends on Reynolds number, "
                    "mesh geometry, and angle",
        "url": "https://doi.org/10.1016/j.aquaeng.2009.04.003",
    },
    "oxygen_solubility": {
        "citation": "Garcia & Gordon (1992)",
        "equation": "O2 saturation from temperature and salinity",
        "envelope": "-2 to 40 C and 0 to 42 PSU",
        "url": "https://doi.org/10.4319/lo.1992.37.6.1307",
    },
    "ammonia": {
        "citation": "Jansen et al. (2012)",
        "equation": "Direct editable excretion prior",
        "envelope": "reported annual range approximately 0.3-2.2 umol N per g "
                    "dry tissue per hour; covaried with food and temperature",
        "url": "https://doi.org/10.1016/j.jembe.2011.11.009",
    },
    "biodeposition": {
        "citation": "Callier et al. (2006)",
        "equation": "Settling-distance screen from depth, transport speed, and "
                    "settling velocity; first-order sediment box",
        "envelope": "Biodeposit production and dispersion vary by mussel age, "
                    "particle size, current, and site; not a seabed model",
        "url": "https://doi.org/10.3354/meps322129",
    },
    "harvest": {
        "citation": "Danish mitigation-mussel values (Taylor et al. 2019); "
                    "UK rope-culture survey (Mascorda Cabre et al. 2021)",
        "equation": "13.7 kg N and 0.9 kg P per harvested wet tonne (Danish); "
                    "8.50 +/- 0.59 kg N and 0.95 +/- 0.07 kg P per live tonne "
                    "(UK) -- composition must stay editable",
        "envelope": "Management accounting prior; NOT predicted growth",
        "url": "https://doi.org/10.3389/fmars.2019.00698",
    },
    "deb_food": {
        "citation": "Maar et al. (2023), coupled hydrodynamic-"
                    "biogeochemical-DEB mussel-farm model",
        "equation": "Size-selective retention, saturating ingestion, "
                    "food-quality assimilation, and separate faeces and "
                    "pseudofaeces pathways",
        "envelope": "Screening translation of published process structure; "
                    "the constants below remain editable priors and require "
                    "site/species calibration",
        "url": "https://doi.org/10.1016/j.scitotenv.2023.164168",
    },
    "oxygen_food": {
        "citation": "Kamermans & Saurel (2022)",
        "equation": "Editable oxygen-activity ramp multiplied into clearance",
        "envelope": "Mytilus edulis responses depend jointly on oxygen, "
                    "temperature and food; this reduced ramp is not a fitted "
                    "dose-response curve",
        "url": "https://doi.org/10.1051/alr/2022001",
    },
}


def _entry(group, default, unit, kind, summary,
           cite=None, case=EXPOSED, bounds=None, bounds_note=None):
    """Build one catalog entry.  ``cite`` is a key into :data:`CITATIONS`."""
    return {
        "group": group,
        "default": default,
        "unit": unit,
        "kind": kind,           # number | int | bool | string | list
        "summary": summary,
        "cite": cite,
        "case": case,
        "bounds": bounds,       # informational (lower, upper); None if open
        "bounds_note": bounds_note,
    }


# --- The grammar book -------------------------------------------------------
# Ordered by section.  Defaults MUST equal musselflow_ecogrammar_core.DEFAULTS;
# the test guards this in both directions.
GRAMMAR = {
    # Schema and profile labels ---------------------------------------------
    "schema.version": _entry(
        "profile", 1, "-", "int",
        "Grammar schema version accepted by the numerical core.",
        case=DERIVED),
    "profile.species": _entry(
        "profile", "Mytilus_edulis_screening", "-", "string",
        "Compatibility label for the built-in species prior set. The case "
        "uses species.taxon instead.",
        case=DEFAULTED),
    "profile.net": _entry(
        "profile", "knotless_screening", "-", "string",
        "Compatibility label for the net-drag prior. The current solver uses "
        "the documented knotless-net relation.",
        cite="net_drag", case=DEFAULTED),

    # Site state and water-column box ---------------------------------------
    "site.depth_m": _entry(
        "site", 12.0, "m", "number",
        "Mean water-column depth of the farm box.",
        case=EXPOSED, bounds=(0.0, None), bounds_note="strictly > 0"),
    "site.temperature_c": _entry(
        "site", 12.0, "degC", "number",
        "Water temperature; drives Q10 activity and oxygen solubility.",
        cite="oxygen_solubility", case=BOUNDARY, bounds=(-2.0, 40.0)),
    "site.salinity_psu": _entry(
        "site", 20.0, "PSU", "number",
        "Salinity; drives the trapezoidal activity profile and O2 solubility.",
        cite="oxygen_solubility", case=BOUNDARY, bounds=(0.0, 42.0)),
    "site.initial_do_mg_l": _entry(
        "site", 9.0, "mg/L", "number",
        "Initial dissolved oxygen in the water-column box.",
        case=EXPOSED, bounds=(0.0, None)),
    "site.boundary_do_mg_l": _entry(
        "site", 9.0, "mg/L", "number",
        "Dissolved oxygen advected in from the domain boundary.",
        case=BOUNDARY, bounds=(0.0, None)),
    "site.chlorophyll_ug_l": _entry(
        "site", 5.0, "ug/L", "number",
        "Boundary chlorophyll-a concentration (phytoplankton food proxy).",
        case=BOUNDARY, bounds=(0.0, None)),
    "site.tsm_mg_l": _entry(
        "site", 3.0, "mg/L", "number",
        "Total suspended particulate matter at the boundary.",
        case=BOUNDARY, bounds=(0.0, None)),
    "site.particulate_organic_fraction": _entry(
        "site", 0.35, "-", "number",
        "Organic fraction of captured particulate matter.",
        case=EXPOSED, bounds=(0.0, 1.0)),
    "site.background_sod_g_o2_m2_day": _entry(
        "site", 0.8, "g O2/m2/day", "number",
        "Background sediment oxygen demand not attributable to the farm.",
        case=EXPOSED, bounds=(0.0, None),
        bounds_note="site prior; replace with benthic measurements"),
    "site.pelagic_respiration_g_o2_m3_day": _entry(
        "site", 0.05, "g O2/m3/day", "number",
        "Water-column (non-mussel) respiration oxygen sink.",
        case=EXPOSED, bounds=(0.0, None)),
    "site.primary_production_g_o2_m3_day": _entry(
        "site", 0.0, "g O2/m3/day", "number",
        "Photosynthetic oxygen source in the box (0 = ignored).",
        case=EXPOSED, bounds=(0.0, None)),
    "site.reaeration_per_day": _entry(
        "site", 0.20, "1/day", "number",
        "First-order surface reaeration rate toward O2 saturation.",
        case=EXPOSED, bounds=(0.0, None)),
    "site.vertical_exchange_per_day": _entry(
        "site", 0.50, "1/day", "number",
        "Vertical mixing coefficient toward boundary DO (stratification "
        "compressed to one term).",
        case=EXPOSED, bounds=(0.0, None)),
    "site.advective_exchange_efficiency": _entry(
        "site", 0.25, "-", "number",
        "Fraction of the box flushed advectively per transit.",
        case=EXPOSED, bounds=(0.0, 1.0)),

    # Food composition and particle classes --------------------------------
    "food.carbon_to_chlorophyll_mg_c_per_mg_chl": _entry(
        "food", 50.0, "mg C/mg Chl-a", "number",
        "Carbon-to-chlorophyll ratio used only for transparent carbon "
        "accounting; it is strongly variable and must be calibrated.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),
    "food.organic_carbon_fraction": _entry(
        "food", 0.45, "kg C/kg organic dry matter", "number",
        "Carbon fraction of captured organic dry matter.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, 1.0)),
    "food.small_particle_fraction": _entry(
        "food", 0.20, "-", "number",
        "Fraction of chlorophyll-bearing particles represented by the small "
        "particle class.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, 1.0)),
    "food.small_particle_diameter_um": _entry(
        "food", 2.0, "um", "number",
        "Representative diameter of the poorly retained small-food class.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),
    "food.large_particle_diameter_um": _entry(
        "food", 15.0, "um", "number",
        "Representative diameter of the efficiently retained food class.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),
    "food.detritus_fraction_of_organic": _entry(
        "food", 0.50, "-", "number",
        "Fraction of particulate organic matter treated as detritus.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, 1.0)),
    "food.detritus_preference": _entry(
        "food", 0.80, "-", "number",
        "Relative selection/retention multiplier for detrital particles.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, 1.0)),
    "food.detritus_assimilation_multiplier": _entry(
        "food", 0.80, "-", "number",
        "Relative assimilation multiplier for detritus versus live food.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, 1.0)),
    "food.organic_energy_kj_g": _entry(
        "food", 20.0, "kJ/g organic dry matter", "number",
        "Energy density used for a transparent scope-for-growth proxy.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),

    # Species priors --------------------------------------------------------
    "species.mean_dry_tissue_g": _entry(
        "species", 0.20, "g", "number",
        "Mean individual dry tissue mass; sets allometric clearance and "
        "respiration.",
        cite="clearance", case=EXPOSED, bounds=(0.0, None),
        bounds_note="strictly > 0; envelope 0.011-1.361 g"),
    "species.size_cv": _entry(
        "species", 0.30, "-", "number",
        "Coefficient of variation of the lognormal size distribution; the "
        "allometric moment is evaluated analytically.",
        cite="clearance", case=EXPOSED, bounds=(0.0, None)),
    "species.clearance_a_l_h": _entry(
        "species", 7.45, "L/h", "number",
        "Clearance-rate allometric coefficient a in a * W**b.",
        cite="clearance", case=EXPOSED, bounds=(0.0, None)),
    "species.clearance_b": _entry(
        "species", 0.66, "-", "number",
        "Clearance-rate allometric exponent b.",
        cite="clearance", case=EXPOSED, bounds=(0.0, None)),
    "species.clearance_ref_temp_c": _entry(
        "species", 12.0, "degC", "number",
        "Reference temperature for the clearance Q10 correction.",
        cite="clearance", case=EXPOSED),
    "species.clearance_q10": _entry(
        "species", 1.50, "-", "number",
        "Editable screening Q10 temperature multiplier for clearance; the "
        "cited allometry does not establish this value.",
        case=EXPOSED, bounds=(0.0, None)),
    "species.respiration_a_ml_o2_h": _entry(
        "species", 0.475, "mL O2/h", "number",
        "Respiration allometric coefficient a in a * W**b.",
        cite="respiration", case=EXPOSED, bounds=(0.0, None)),
    "species.respiration_b": _entry(
        "species", 0.663, "-", "number",
        "Respiration allometric exponent b.",
        cite="respiration", case=EXPOSED, bounds=(0.0, None)),
    "species.respiration_ref_temp_c": _entry(
        "species", 12.0, "degC", "number",
        "Reference temperature for the respiration Q10 correction.",
        cite="respiration", case=EXPOSED),
    "species.respiration_q10": _entry(
        "species", 2.00, "-", "number",
        "Editable screening Q10 temperature multiplier for respiration; the "
        "cited allometry does not establish this value.",
        case=EXPOSED, bounds=(0.0, None)),
    "species.retention_efficiency": _entry(
        "species", 0.80, "-", "number",
        "Retention efficiency for phytoplankton/chlorophyll particles.",
        cite="retention", case=EXPOSED, bounds=(0.0, 1.0)),
    "species.particulate_retention_efficiency": _entry(
        "species", 0.80, "-", "number",
        "Retention efficiency for general suspended particulate matter.",
        cite="retention", case=EXPOSED, bounds=(0.0, 1.0)),
    "species.assimilation_efficiency": _entry(
        "species", 0.70, "-", "number",
        "Fraction of ingested organic matter assimilated into tissue.",
        case=EXPOSED, bounds=(0.0, 1.0)),
    "species.pseudofaeces_fraction": _entry(
        "species", 0.10, "-", "number",
        "Fraction of captured material rejected as pseudofaeces.",
        case=EXPOSED, bounds=(0.0, 1.0)),
    "species.ammonia_mg_n_g_dw_h": _entry(
        "species", 0.015, "mg N/g dw/h", "number",
        "Ammonia (dissolved inorganic N) excretion rate per dry tissue.",
        cite="ammonia", case=EXPOSED, bounds=(0.0, None)),
    "species.activity_fraction": _entry(
        "species", 1.00, "-", "number",
        "Baseline fraction of time the animals actively filter.",
        case=EXPOSED, bounds=(0.0, 1.0)),
    "species.low_food_threshold_ug_l": _entry(
        "species", 0.70, "ug/L", "number",
        "Chlorophyll below which filtration begins to down-regulate.",
        case=EXPOSED, bounds=(0.0, None)),
    "species.low_food_transition_ug_l": _entry(
        "species", 0.15, "ug/L", "number",
        "Width of the low-food down-regulation transition.",
        case=EXPOSED, bounds=(0.0, None), bounds_note="strictly > 0"),
    "species.retention_d50_um": _entry(
        "species", 5.0, "um", "number",
        "Particle diameter at 50 percent of the editable retention curve.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),
    "species.retention_slope_um": _entry(
        "species", 1.5, "um", "number",
        "Width of the particle-size retention transition.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),
    "species.ingestion_half_saturation_ug_l": _entry(
        "species", 17.0, "ug Chl-a/L", "number",
        "Half-saturation scale used to prevent clearance-based ingestion "
        "from increasing without bound at high food.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),
    "species.high_food_assimilation_threshold_ug_l": _entry(
        "species", 17.0, "ug Chl-a/L", "number",
        "Food concentration above which assimilation begins to decline.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),
    "species.high_food_assimilation_decay_ug_l": _entry(
        "species", 8.0, "ug Chl-a/L", "number",
        "Decay scale for reduced assimilation above the high-food threshold.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),
    "species.assimilation_quality_half_saturation": _entry(
        "species", 0.15, "organic fraction", "number",
        "Half-saturation of the food-quality assimilation multiplier.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, 1.0)),
    "species.assimilation_reference_organic_fraction": _entry(
        "species", 0.35, "organic fraction", "number",
        "Reference organic fraction at which baseline assimilation is "
        "recovered.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, 1.0)),
    "species.pseudofaeces_tsm_threshold_mg_l": _entry(
        "species", 4.0, "mg/L", "number",
        "TSM concentration above which rejection as pseudofaeces increases.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),
    "species.pseudofaeces_tsm_transition_mg_l": _entry(
        "species", 2.0, "mg/L", "number",
        "Width of the TSM-dependent pseudofaeces transition.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, None)),
    "species.pseudofaeces_max_fraction": _entry(
        "species", 0.60, "-", "number",
        "Upper screening bound on the captured fraction rejected as "
        "pseudofaeces.",
        cite="deb_food", case=EXPOSED, bounds=(0.0, 1.0)),
    "species.oxygen_zero_saturation_fraction": _entry(
        "species", 0.20, "fraction O2 saturation", "number",
        "Oxygen saturation fraction at/below which filtration is zero.",
        cite="oxygen_food", case=EXPOSED, bounds=(0.0, 1.0)),
    "species.oxygen_full_saturation_fraction": _entry(
        "species", 0.70, "fraction O2 saturation", "number",
        "Oxygen saturation fraction at/above which filtration is unaffected.",
        cite="oxygen_food", case=EXPOSED, bounds=(0.0, 1.0)),
    "species.current_clearance_start_m_s": _entry(
        "species", 0.20, "m/s", "number",
        "Current speed above which unprotected small aggregations begin to "
        "lose clearance.",
        cite="current", case=EXPOSED, bounds=(0.0, None)),
    "species.current_clearance_zero_m_s": _entry(
        "species", 0.60, "m/s", "number",
        "Current speed at which the small-aggregation response reaches zero.",
        cite="current", case=EXPOSED, bounds=(0.0, None)),
    "species.current_protection_group_size": _entry(
        "species", 20.0, "animals", "number",
        "Aggregation size assigned full current protection in the screening "
        "interpolation.",
        cite="current", case=EXPOSED, bounds=(0.0, None)),
    "species.oxycalorific_kj_per_ml_o2": _entry(
        "species", 0.0201, "kJ/mL O2", "number",
        "Oxycalorific conversion used in the energy-balance proxy.",
        case=EXPOSED, bounds=(0.0, None)),
    "species.tissue_energy_kj_g_dw": _entry(
        "species", 23.5, "kJ/g dry tissue", "number",
        "Tissue-energy conversion used only for potential-growth screening.",
        case=EXPOSED, bounds=(0.0, None)),
    "species.salinity_zero_low_psu": _entry(
        "species", 4.0, "PSU", "number",
        "Salinity at/below which activity is zero (low end of trapezoid).",
        case=EXPOSED, bounds=(0.0, None),
        bounds_note="editable screening prior; calibrate locally"),
    "species.salinity_full_low_psu": _entry(
        "species", 15.0, "PSU", "number",
        "Salinity above which activity is full (low breakpoint).",
        case=EXPOSED, bounds=(0.0, None)),
    "species.salinity_full_high_psu": _entry(
        "species", 32.0, "PSU", "number",
        "Salinity below which activity is full (high breakpoint).",
        case=EXPOSED, bounds=(0.0, None)),
    "species.salinity_zero_high_psu": _entry(
        "species", 42.0, "PSU", "number",
        "Salinity at/above which activity is zero (high end of trapezoid).",
        case=EXPOSED, bounds=(0.0, None)),
    "species.valid_flow_max_m_s": _entry(
        "species", 1.40, "m/s", "number",
        "Upper current speed of the cited clearance study envelope. Currents "
        "above it are constraint-infeasible unless allow_extrapolation is set.",
        cite="current", case=EXPOSED, bounds=(0.0, None),
        bounds_note="strictly > 0"),

    # Stocking --------------------------------------------------------------
    "stocking.mode": _entry(
        "stocking", "animals", "-", "string",
        "Whether stock is given as animal counts or dry biomass. The current "
        "strict case accepts only 'animals'.",
        case=EXPOSED, bounds_note="animals | dry_biomass"),
    "stocking.mussels_per_obstacle": _entry(
        "stocking", [1000.0], "count", "list",
        "Animal count per obstacle (scalar broadcasts to every obstacle).",
        case=EXPOSED, bounds=(0.0, None)),
    "stocking.dry_tissue_kg_per_obstacle": _entry(
        "stocking", [-1.0], "kg", "list",
        "Dry tissue mass per obstacle for the dry_biomass mode. The current "
        "strict case does not expose that mode.",
        case=DROPPED, bounds_note="-1 = unused, else >= 0"),
    "stocking.live_wet_g_per_individual": _entry(
        "stocking", 20.0, "g", "number",
        "Live wet mass per individual, for harvest accounting.",
        case=EXPOSED, bounds=(0.0, None), bounds_note="strictly > 0"),
    "stocking.annual_mortality_fraction": _entry(
        "stocking", 0.10, "1/year", "number",
        "Annual mortality fraction; feeds organic loading and survival.",
        case=EXPOSED, bounds=(0.0, 1.0)),
    "stocking.effective_aggregation_size": _entry(
        "stocking", 20.0, "animals", "number",
        "Effective group size controlling current-speed protection of "
        "clearance; this is not the total obstacle population.",
        cite="current", case=EXPOSED, bounds=(0.0, None)),

    # Porous structure and reduced wake -------------------------------------
    "structure.porosity": _entry(
        "structure", [0.70], "-", "list",
        "Open-area fraction per obstacle (solidity S = 1 - porosity).",
        cite="net_drag", case=EXPOSED, bounds=(0.0, 0.99)),
    "structure.twine_diameter_m": _entry(
        "structure", 0.003, "m", "number",
        "Net twine diameter; sets the reported twine Reynolds number.",
        cite="net_drag", case=EXPOSED, bounds=(0.0, None), bounds_note=">0"),
    "structure.drag_multiplier": _entry(
        "structure", 1.0, "-", "number",
        "Scalar multiplier on the net-drag prior for sensitivity studies.",
        cite="net_drag", case=EXPOSED, bounds=(0.0, None), bounds_note=">0"),
    "structure.fallback_plan_size_m": _entry(
        "structure", 0.10, "m", "number",
        "Plan dimension assumed when an obstacle's descriptor is degenerate.",
        case=EXPOSED, bounds=(0.0, None), bounds_note=">0"),
    "structure.fallback_height_m": _entry(
        "structure", 1.00, "m", "number",
        "Height assumed when an obstacle's descriptor is degenerate.",
        case=EXPOSED, bounds=(0.0, None), bounds_note=">0"),

    # Hydrodynamic screening parameters -------------------------------------
    "hydrodynamics.wake_spread": _entry(
        "hydrodynamics", 0.12, "-", "number",
        "Lateral spread rate of the algebraic Jensen-style wake. Model tuning "
        "parameter, not a measured constant.",
        case=EXPOSED, bounds=(0.0, None), bounds_note=">0"),
    "hydrodynamics.food_plume_spread": _entry(
        "hydrodynamics", 0.16, "-", "number",
        "Lateral spread rate of the food-depletion plume. Model parameter.",
        case=EXPOSED, bounds=(0.0, None), bounds_note=">0"),
    "hydrodynamics.food_recovery_lengths": _entry(
        "hydrodynamics", 8.0, "-", "number",
        "Downstream length scales over which food recovers. Model parameter.",
        case=EXPOSED, bounds=(0.0, None), bounds_note=">0"),
    "hydrodynamics.min_speed_ratio": _entry(
        "hydrodynamics", 0.02, "-", "number",
        "Floor on local/free-stream speed ratio to keep the algebra finite.",
        case=EXPOSED, bounds=(0.0, 1.0)),
    "hydrodynamics.kinematic_viscosity_m2_s": _entry(
        "hydrodynamics", 1.30e-6, "m2/s", "number",
        "Seawater kinematic viscosity for Reynolds numbers.",
        case=EXPOSED, bounds=(0.0, None), bounds_note=">0"),
    "hydrodynamics.water_density_kg_m3": _entry(
        "hydrodynamics", 1020.0, "kg/m3", "number",
        "Seawater density reserved for physical drag-force reporting; the "
        "current wake ratio does not consume it.",
        case=EXPOSED, bounds=(0.0, None), bounds_note=">0"),

    # Scenario combination --------------------------------------------------
    "scenario.weights": _entry(
        "scenario", [1.0], "-", "list",
        "Relative weight per current scenario. In the strict case these are derived from "
        "step durations, not set directly.",
        case=DERIVED, bounds=(0.0, None)),
    "scenario.duration_h": _entry(
        "scenario", [24.0], "h", "list",
        "Duration each current scenario is held while advancing oxygen and "
        "sediment stocks.",
        case=EXPOSED, bounds=(0.0, None), bounds_note="strictly > 0"),
    "scenario.repeat_count": _entry(
        "scenario", 1, "-", "int",
        "Number of times the ordered timeline repeats.",
        case=EXPOSED, bounds=(1, None)),

    # Organic deposition and sediment oxygen demand -------------------------
    "sediment.settling_velocity_m_s": _entry(
        "sediment", 0.006, "m/s", "number",
        "Biodeposit settling velocity for the settling-distance screen.",
        cite="biodeposition", case=EXPOSED, bounds=(0.0, None)),
    "sediment.in_domain_deposition_fraction": _entry(
        "sediment", -1.0, "-", "number",
        "Fraction of biodeposits landing inside the domain. -1 = compute from "
        "a settling-distance screen instead of a seabed transport model.",
        cite="biodeposition", case=EXPOSED,
        bounds_note="-1 (automatic) or 0..1"),
    "sediment.oxygen_demand_kg_o2_per_kg_organic": _entry(
        "sediment", 1.0, "kg O2/kg", "number",
        "Oxygen demand per kilogram of deposited organic matter.",
        cite="biodeposition", case=EXPOSED, bounds=(0.0, None)),
    "sediment.decay_per_day": _entry(
        "sediment", 0.05, "1/day", "number",
        "First-order decay of the deposited organic stock.",
        cite="biodeposition", case=EXPOSED, bounds=(0.0, None)),
    "sediment.resuspension_per_day": _entry(
        "sediment", 0.0, "1/day", "number",
        "First-order resuspension of the deposited organic stock.",
        cite="biodeposition", case=EXPOSED, bounds=(0.0, None)),
    "sediment.mortality_deposition_fraction": _entry(
        "sediment", 1.0, "-", "number",
        "Fraction of mortality organic loss deposited to the seabed.",
        cite="biodeposition", case=EXPOSED, bounds=(0.0, 1.0)),
    "sediment.initial_organic_stock_kg": _entry(
        "sediment", 0.0, "kg", "number",
        "Initial deposited organic stock at the start of a run.",
        cite="biodeposition", case=EXPOSED, bounds=(0.0, None)),

    # Harvest accounting (management scenario, not predicted growth) --------
    "harvest.fraction_per_year": _entry(
        "harvest", 0.80, "1/year", "number",
        "Fraction of standing stock harvested per year.",
        cite="harvest", case=EXPOSED, bounds=(0.0, 1.0)),
    "harvest.turnovers_per_year": _entry(
        "harvest", 1.0, "1/year", "number",
        "Harvest turnovers per year for the accounting scenario.",
        cite="harvest", case=EXPOSED, bounds=(0.0, None)),
    "harvest.n_kg_per_t_wet": _entry(
        "harvest", 13.7, "kg N/t", "number",
        "Nitrogen removed per harvested wet tonne (editable composition).",
        cite="harvest", case=EXPOSED, bounds=(0.0, None)),
    "harvest.p_kg_per_t_wet": _entry(
        "harvest", 0.9, "kg P/t", "number",
        "Phosphorus removed per harvested wet tonne (editable composition).",
        cite="harvest", case=EXPOSED, bounds=(0.0, None)),

    # Hard constraints (negative disables) ----------------------------------
    "constraint.min_do_mg_l": _entry(
        "constraint", -1.0, "mg/L", "number",
        "Minimum modeled dissolved oxygen. Negative disables the constraint; "
        "no legal threshold is asserted by the program.",
        case=EXPOSED, bounds_note="-1 disables, else >= 0"),
    "constraint.min_probe_speed_m_s": _entry(
        "constraint", -1.0, "m/s", "number",
        "Minimum required probe/flushing speed. Negative disables.",
        case=EXPOSED, bounds_note="-1 disables, else >= 0"),
    "constraint.max_biodeposition_kg_m2_day": _entry(
        "constraint", -1.0, "kg/m2/day", "number",
        "Maximum allowed biodeposition. Negative disables.",
        case=EXPOSED, bounds_note="-1 disables, else >= 0"),
    "constraint.min_obstacle_clearance_m": _entry(
        "constraint", 0.0, "m", "number",
        "Minimum spacing between obstacles (collision constraint).",
        case=EXPOSED, bounds=(0.0, None)),

    # Objective target and weights ------------------------------------------
    "objective.target_chlorophyll_capture_g_day": _entry(
        "objective", -1.0, "g/day", "number",
        "Absolute chlorophyll-capture target. Negative uses the design's own "
        "theoretical clearance capacity for the extraction score.",
        case=EXPOSED, bounds_note="-1 = auto, else >= 0"),
    "objective.target_biodeposition_kg_m2_day": _entry(
        "objective", 0.05, "kg/m2/day", "number",
        "Reference biodeposition used to normalise the low-deposition score.",
        case=EXPOSED, bounds=(0.0, None), bounds_note="strictly > 0"),
    "weight.extraction": _entry(
        "objective", 0.30, "-", "number",
        "Objective weight: particulate/chlorophyll extraction.",
        case=EXPOSED, bounds=(0.0, None)),
    "weight.flushing": _entry(
        "objective", 0.20, "-", "number",
        "Objective weight: flushing/water exchange.",
        case=EXPOSED, bounds=(0.0, None)),
    "weight.food_delivery": _entry(
        "objective", 0.15, "-", "number",
        "Objective weight: food delivery to obstacles.",
        case=EXPOSED, bounds=(0.0, None)),
    "weight.uniformity": _entry(
        "objective", 0.10, "-", "number",
        "Objective weight: uniformity of food delivery.",
        case=EXPOSED, bounds=(0.0, None)),
    "weight.oxygen": _entry(
        "objective", 0.15, "-", "number",
        "Objective weight: dissolved-oxygen safety proxy.",
        case=EXPOSED, bounds=(0.0, None)),
    "weight.low_deposition": _entry(
        "objective", 0.10, "-", "number",
        "Objective weight: low biodeposition proxy.",
        case=EXPOSED, bounds=(0.0, None)),

    # Governance / run envelope ---------------------------------------------
    # Calibration flags never change equations; they are governance
    # declarations and are derived (forced false), never user-set.
    "validation.geometry_calibrated": _entry(
        "validation", False, "-", "bool",
        "Governance flag: geometry descriptors validated against known shapes.",
        case=DERIVED),
    "validation.hydrodynamics_calibrated": _entry(
        "validation", False, "-", "bool",
        "Governance flag: wake/drag validated against current measurements.",
        case=DERIVED),
    "validation.biology_calibrated": _entry(
        "validation", False, "-", "bool",
        "Governance flag: clearance/capture validated against measurements.",
        case=DERIVED),
    "validation.oxygen_calibrated": _entry(
        "validation", False, "-", "bool",
        "Governance flag: oxygen box validated over tidal/low-flow sequences.",
        case=DERIVED),
    "validation.sediment_calibrated": _entry(
        "validation", False, "-", "bool",
        "Governance flag: deposition validated against traps and benthic flux.",
        case=DERIVED),
    "validation.allow_extrapolation": _entry(
        "validation", False, "-", "bool",
        "If true, use-envelope violations become soft warnings instead of a "
        "hard infeasible constraint. A deliberate, declared choice.",
        case=EXPOSED),
}


# --- Accessors --------------------------------------------------------------
def keys():
    """Return all grammar keys in catalog order."""
    return list(GRAMMAR)


def entry(key):
    """Return the full catalog entry for one key."""
    return GRAMMAR[key]


def defaults():
    """Return ``{key: default}`` derived from the book (a fresh dict)."""
    return {key: item["default"] for key, item in GRAMMAR.items()}


def groups():
    """Return grammar keys grouped by section, preserving catalog order."""
    grouped = {}
    for key, item in GRAMMAR.items():
        grouped.setdefault(item["group"], []).append(key)
    return grouped


def provenance(key):
    """Return the resolved citation record for a key, or ``None``."""
    cite = GRAMMAR[key]["cite"]
    return dict(CITATIONS[cite]) if cite is not None else None


def case_gaps():
    """Return keys the current strict case cannot drive."""
    return [key for key, item in GRAMMAR.items() if item["case"] in CASE_MISSING]


def render_markdown():
    """Render the grammar book as a human-readable Markdown document."""
    lines = ["# MusselFlow ecological grammar book", ""]
    lines.append(
        "Screening priors, not universal constants. Every active value is "
        "editable and reported. `Case` records how the strict ecological "
        "grammar supplies each solver field.")
    lines.append("")
    status_note = {
        EXPOSED: "settable in musselflow_grammar.json",
        BOUNDARY: "settable via a forcing-step boundary",
        DERIVED: "computed by the runtime",
        DEFAULTED: "compatibility metadata/default",
        DROPPED: "not represented by the current case",
    }
    for group, group_keys in groups().items():
        lines.append("## %s" % group)
        lines.append("")
        for key in group_keys:
            item = GRAMMAR[key]
            lines.append("### `%s`" % key)
            lines.append("")
            lines.append("- **default:** `%r` %s" % (
                item["default"], item["unit"]))
            lines.append("- **summary:** %s" % item["summary"])
            if item["bounds_note"]:
                lines.append("- **range:** %s" % item["bounds_note"])
            elif item["bounds"] is not None:
                lower, upper = item["bounds"]
                lines.append("- **range:** %s .. %s" % (
                    "-inf" if lower is None else lower,
                    "+inf" if upper is None else upper))
            lines.append("- **Case:** %s (%s)" % (
                item["case"], status_note[item["case"]]))
            record = provenance(key)
            if record is not None:
                lines.append("- **equation:** %s" % record["equation"])
                lines.append("- **envelope:** %s" % record["envelope"])
                lines.append("- **source:** %s — %s" % (
                    record["citation"], record["url"]))
            lines.append("")
    return "\n".join(lines)
