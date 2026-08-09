# MusselFlow ecological grammar book

Screening priors, not universal constants. Every active value is editable and reported. `Case` records how the strict ecological grammar supplies each solver field.

## profile

### `schema.version`

- **default:** `1` -
- **summary:** Grammar schema version accepted by the numerical core.
- **Case:** derived (computed by the runtime)

### `profile.species`

- **default:** `'Mytilus_edulis_screening'` -
- **summary:** Compatibility label for the built-in species prior set. The case uses species.taxon instead.
- **Case:** defaulted (compatibility metadata/default)

### `profile.net`

- **default:** `'knotless_screening'` -
- **summary:** Compatibility label for the net-drag prior. The current solver uses the documented knotless-net relation.
- **Case:** defaulted (compatibility metadata/default)
- **equation:** Cd_solid = 1 + 1.37*S + 0.78*S**2 for a knotless net, with S = 1 - porosity
- **envelope:** Net-panel empirical relation; depends on Reynolds number, mesh geometry, and angle
- **source:** Lader et al. (2009), aquaculture net drag — https://doi.org/10.1016/j.aquaeng.2009.04.003

## site

### `site.depth_m`

- **default:** `12.0` m
- **summary:** Mean water-column depth of the farm box.
- **range:** strictly > 0
- **Case:** exposed (settable in musselflow_grammar.json)

### `site.temperature_c`

- **default:** `12.0` degC
- **summary:** Water temperature; drives Q10 activity and oxygen solubility.
- **range:** -2.0 .. 40.0
- **Case:** boundary (settable via a forcing-step boundary)
- **equation:** O2 saturation from temperature and salinity
- **envelope:** -2 to 40 C and 0 to 42 PSU
- **source:** Garcia & Gordon (1992) — https://doi.org/10.4319/lo.1992.37.6.1307

### `site.salinity_psu`

- **default:** `20.0` PSU
- **summary:** Salinity; drives the trapezoidal activity profile and O2 solubility.
- **range:** 0.0 .. 42.0
- **Case:** boundary (settable via a forcing-step boundary)
- **equation:** O2 saturation from temperature and salinity
- **envelope:** -2 to 40 C and 0 to 42 PSU
- **source:** Garcia & Gordon (1992) — https://doi.org/10.4319/lo.1992.37.6.1307

### `site.initial_do_mg_l`

- **default:** `9.0` mg/L
- **summary:** Initial dissolved oxygen in the water-column box.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `site.boundary_do_mg_l`

- **default:** `9.0` mg/L
- **summary:** Dissolved oxygen advected in from the domain boundary.
- **range:** 0.0 .. +inf
- **Case:** boundary (settable via a forcing-step boundary)

### `site.chlorophyll_ug_l`

- **default:** `5.0` ug/L
- **summary:** Boundary chlorophyll-a concentration (phytoplankton food proxy).
- **range:** 0.0 .. +inf
- **Case:** boundary (settable via a forcing-step boundary)

### `site.tsm_mg_l`

- **default:** `3.0` mg/L
- **summary:** Total suspended particulate matter at the boundary.
- **range:** 0.0 .. +inf
- **Case:** boundary (settable via a forcing-step boundary)

### `site.particulate_organic_fraction`

- **default:** `0.35` -
- **summary:** Organic fraction of captured particulate matter.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)

### `site.background_sod_g_o2_m2_day`

- **default:** `0.8` g O2/m2/day
- **summary:** Background sediment oxygen demand not attributable to the farm.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** O2 saturation from temperature and salinity
- **envelope:** -2 to 40 C and 0 to 42 PSU
- **source:** Garcia & Gordon (1992) — https://doi.org/10.4319/lo.1992.37.6.1307

### `site.pelagic_respiration_g_o2_m3_day`

- **default:** `0.05` g O2/m3/day
- **summary:** Water-column (non-mussel) respiration oxygen sink.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `site.primary_production_g_o2_m3_day`

- **default:** `0.0` g O2/m3/day
- **summary:** Photosynthetic oxygen source in the box (0 = ignored).
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `site.reaeration_per_day`

- **default:** `0.2` 1/day
- **summary:** First-order surface reaeration rate toward O2 saturation.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `site.vertical_exchange_per_day`

- **default:** `0.5` 1/day
- **summary:** Vertical mixing coefficient toward boundary DO (stratification compressed to one term).
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `site.advective_exchange_efficiency`

- **default:** `0.25` -
- **summary:** Fraction of the box flushed advectively per transit.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)

## species

### `species.mean_dry_tissue_g`

- **default:** `0.2` g
- **summary:** Mean individual dry tissue mass; sets allometric clearance and respiration.
- **range:** strictly > 0; envelope 0.011-1.361 g
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** clearance_L_h = 7.45 * dry_tissue_g ** 0.66
- **envelope:** adult dry tissue 0.011-1.361 g; approximately 10-13 C, 30 PSU
- **source:** Moehlenberg & Riisgaard (1979), adult Mytilus edulis; synthesis Riisgaard et al. (2025) — https://doi.org/10.1242/bio.062024

### `species.size_cv`

- **default:** `0.3` -
- **summary:** Coefficient of variation of the lognormal size distribution; the allometric moment is evaluated analytically.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** clearance_L_h = 7.45 * dry_tissue_g ** 0.66
- **envelope:** adult dry tissue 0.011-1.361 g; approximately 10-13 C, 30 PSU
- **source:** Moehlenberg & Riisgaard (1979), adult Mytilus edulis; synthesis Riisgaard et al. (2025) — https://doi.org/10.1242/bio.062024

### `species.clearance_a_l_h`

- **default:** `7.45` L/h
- **summary:** Clearance-rate allometric coefficient a in a * W**b.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** clearance_L_h = 7.45 * dry_tissue_g ** 0.66
- **envelope:** adult dry tissue 0.011-1.361 g; approximately 10-13 C, 30 PSU
- **source:** Moehlenberg & Riisgaard (1979), adult Mytilus edulis; synthesis Riisgaard et al. (2025) — https://doi.org/10.1242/bio.062024

### `species.clearance_b`

- **default:** `0.66` -
- **summary:** Clearance-rate allometric exponent b.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** clearance_L_h = 7.45 * dry_tissue_g ** 0.66
- **envelope:** adult dry tissue 0.011-1.361 g; approximately 10-13 C, 30 PSU
- **source:** Moehlenberg & Riisgaard (1979), adult Mytilus edulis; synthesis Riisgaard et al. (2025) — https://doi.org/10.1242/bio.062024

### `species.clearance_ref_temp_c`

- **default:** `12.0` degC
- **summary:** Reference temperature for the clearance Q10 correction.
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** clearance_L_h = 7.45 * dry_tissue_g ** 0.66
- **envelope:** adult dry tissue 0.011-1.361 g; approximately 10-13 C, 30 PSU
- **source:** Moehlenberg & Riisgaard (1979), adult Mytilus edulis; synthesis Riisgaard et al. (2025) — https://doi.org/10.1242/bio.062024

### `species.clearance_q10`

- **default:** `1.5` -
- **summary:** Q10 temperature multiplier for clearance.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** clearance_L_h = 7.45 * dry_tissue_g ** 0.66
- **envelope:** adult dry tissue 0.011-1.361 g; approximately 10-13 C, 30 PSU
- **source:** Moehlenberg & Riisgaard (1979), adult Mytilus edulis; synthesis Riisgaard et al. (2025) — https://doi.org/10.1242/bio.062024

### `species.respiration_a_ml_o2_h`

- **default:** `0.475` mL O2/h
- **summary:** Respiration allometric coefficient a in a * W**b.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** respiration_mL_O2_h = 0.475 * dry_tissue_g ** 0.663
- **envelope:** adult dry tissue; approximately 10-13 C, 30 PSU
- **source:** Hamburger et al. (1983), adult Mytilus edulis; synthesis Riisgaard et al. (2025) — https://doi.org/10.1242/bio.062024

### `species.respiration_b`

- **default:** `0.663` -
- **summary:** Respiration allometric exponent b.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** respiration_mL_O2_h = 0.475 * dry_tissue_g ** 0.663
- **envelope:** adult dry tissue; approximately 10-13 C, 30 PSU
- **source:** Hamburger et al. (1983), adult Mytilus edulis; synthesis Riisgaard et al. (2025) — https://doi.org/10.1242/bio.062024

### `species.respiration_ref_temp_c`

- **default:** `12.0` degC
- **summary:** Reference temperature for the respiration Q10 correction.
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** respiration_mL_O2_h = 0.475 * dry_tissue_g ** 0.663
- **envelope:** adult dry tissue; approximately 10-13 C, 30 PSU
- **source:** Hamburger et al. (1983), adult Mytilus edulis; synthesis Riisgaard et al. (2025) — https://doi.org/10.1242/bio.062024

### `species.respiration_q10`

- **default:** `2.0` -
- **summary:** Q10 temperature multiplier for respiration.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** respiration_mL_O2_h = 0.475 * dry_tissue_g ** 0.663
- **envelope:** adult dry tissue; approximately 10-13 C, 30 PSU
- **source:** Hamburger et al. (1983), adult Mytilus edulis; synthesis Riisgaard et al. (2025) — https://doi.org/10.1242/bio.062024

### `species.retention_efficiency`

- **default:** `0.8` -
- **summary:** Retention efficiency for phytoplankton/chlorophyll particles.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Direct editable retention efficiency (no size resolution)
- **envelope:** Natural seston; strong particle-size and temporal variability (peak retention near 30-35 um)
- **source:** Strohmeier et al. (2012) — https://doi.org/10.1016/j.jembe.2011.11.006

### `species.particulate_retention_efficiency`

- **default:** `0.8` -
- **summary:** Retention efficiency for general suspended particulate matter.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Direct editable retention efficiency (no size resolution)
- **envelope:** Natural seston; strong particle-size and temporal variability (peak retention near 30-35 um)
- **source:** Strohmeier et al. (2012) — https://doi.org/10.1016/j.jembe.2011.11.006

### `species.assimilation_efficiency`

- **default:** `0.7` -
- **summary:** Fraction of ingested organic matter assimilated into tissue.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)

### `species.pseudofaeces_fraction`

- **default:** `0.1` -
- **summary:** Fraction of captured material rejected as pseudofaeces.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)

### `species.ammonia_mg_n_g_dw_h`

- **default:** `0.015` mg N/g dw/h
- **summary:** Ammonia (dissolved inorganic N) excretion rate per dry tissue.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Direct editable excretion prior
- **envelope:** reported annual range approximately 0.3-2.2 umol N per g dry tissue per hour; covaried with food and temperature
- **source:** Jansen et al. (2012) — https://doi.org/10.1016/j.jembe.2011.11.009

### `species.activity_fraction`

- **default:** `1.0` -
- **summary:** Baseline fraction of time the animals actively filter.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)

### `species.low_food_threshold_ug_l`

- **default:** `0.7` ug/L
- **summary:** Chlorophyll below which filtration begins to down-regulate.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `species.low_food_transition_ug_l`

- **default:** `0.15` ug/L
- **summary:** Width of the low-food down-regulation transition.
- **range:** strictly > 0
- **Case:** exposed (settable in musselflow_grammar.json)

### `species.salinity_zero_low_psu`

- **default:** `4.0` PSU
- **summary:** Salinity at/below which activity is zero (low end of trapezoid).
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** No universal multiplier; speeds above the envelope are flagged, not silently penalised
- **envelope:** 0.05-1.4 m/s; aggregation strongly changed the response (groups of 20 held clearance where groups of 3 did not)
- **source:** Nielsen & Vismann (2014) — https://doi.org/10.2983/035.033.0214

### `species.salinity_full_low_psu`

- **default:** `15.0` PSU
- **summary:** Salinity above which activity is full (low breakpoint).
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `species.salinity_full_high_psu`

- **default:** `32.0` PSU
- **summary:** Salinity below which activity is full (high breakpoint).
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `species.salinity_zero_high_psu`

- **default:** `42.0` PSU
- **summary:** Salinity at/above which activity is zero (high end of trapezoid).
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `species.valid_flow_max_m_s`

- **default:** `1.4` m/s
- **summary:** Upper current speed of the cited clearance study envelope. Currents above it are constraint-infeasible unless allow_extrapolation is set.
- **range:** strictly > 0
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** No universal multiplier; speeds above the envelope are flagged, not silently penalised
- **envelope:** 0.05-1.4 m/s; aggregation strongly changed the response (groups of 20 held clearance where groups of 3 did not)
- **source:** Nielsen & Vismann (2014) — https://doi.org/10.2983/035.033.0214

## stocking

### `stocking.mode`

- **default:** `'animals'` -
- **summary:** Whether stock is given as animal counts or dry biomass. The current strict case accepts only 'animals'.
- **range:** animals | dry_biomass
- **Case:** exposed (settable in musselflow_grammar.json)

### `stocking.mussels_per_obstacle`

- **default:** `[1000.0]` count
- **summary:** Animal count per obstacle (scalar broadcasts to every obstacle).
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `stocking.dry_tissue_kg_per_obstacle`

- **default:** `[-1.0]` kg
- **summary:** Dry tissue mass per obstacle for the dry_biomass mode. The current strict case does not expose that mode.
- **range:** -1 = unused, else >= 0
- **Case:** dropped (not represented by the current case)

### `stocking.live_wet_g_per_individual`

- **default:** `20.0` g
- **summary:** Live wet mass per individual, for harvest accounting.
- **range:** strictly > 0
- **Case:** exposed (settable in musselflow_grammar.json)

### `stocking.annual_mortality_fraction`

- **default:** `0.1` 1/year
- **summary:** Annual mortality fraction; feeds organic loading and survival.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)

## structure

### `structure.porosity`

- **default:** `[0.7]` -
- **summary:** Open-area fraction per obstacle (solidity S = 1 - porosity).
- **range:** 0.0 .. 0.99
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Cd_solid = 1 + 1.37*S + 0.78*S**2 for a knotless net, with S = 1 - porosity
- **envelope:** Net-panel empirical relation; depends on Reynolds number, mesh geometry, and angle
- **source:** Lader et al. (2009), aquaculture net drag — https://doi.org/10.1016/j.aquaeng.2009.04.003

### `structure.twine_diameter_m`

- **default:** `0.003` m
- **summary:** Net twine diameter; sets the reported twine Reynolds number.
- **range:** >0
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Cd_solid = 1 + 1.37*S + 0.78*S**2 for a knotless net, with S = 1 - porosity
- **envelope:** Net-panel empirical relation; depends on Reynolds number, mesh geometry, and angle
- **source:** Lader et al. (2009), aquaculture net drag — https://doi.org/10.1016/j.aquaeng.2009.04.003

### `structure.drag_multiplier`

- **default:** `1.0` -
- **summary:** Scalar multiplier on the net-drag prior for sensitivity studies.
- **range:** >0
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Cd_solid = 1 + 1.37*S + 0.78*S**2 for a knotless net, with S = 1 - porosity
- **envelope:** Net-panel empirical relation; depends on Reynolds number, mesh geometry, and angle
- **source:** Lader et al. (2009), aquaculture net drag — https://doi.org/10.1016/j.aquaeng.2009.04.003

### `structure.fallback_plan_size_m`

- **default:** `0.1` m
- **summary:** Plan dimension assumed when an obstacle's descriptor is degenerate.
- **range:** >0
- **Case:** exposed (settable in musselflow_grammar.json)

### `structure.fallback_height_m`

- **default:** `1.0` m
- **summary:** Height assumed when an obstacle's descriptor is degenerate.
- **range:** >0
- **Case:** exposed (settable in musselflow_grammar.json)

## hydrodynamics

### `hydrodynamics.wake_spread`

- **default:** `0.12` -
- **summary:** Lateral spread rate of the algebraic Jensen-style wake. Model tuning parameter, not a measured constant.
- **range:** >0
- **Case:** exposed (settable in musselflow_grammar.json)

### `hydrodynamics.food_plume_spread`

- **default:** `0.16` -
- **summary:** Lateral spread rate of the food-depletion plume. Model parameter.
- **range:** >0
- **Case:** exposed (settable in musselflow_grammar.json)

### `hydrodynamics.food_recovery_lengths`

- **default:** `8.0` -
- **summary:** Downstream length scales over which food recovers. Model parameter.
- **range:** >0
- **Case:** exposed (settable in musselflow_grammar.json)

### `hydrodynamics.min_speed_ratio`

- **default:** `0.02` -
- **summary:** Floor on local/free-stream speed ratio to keep the algebra finite.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)

### `hydrodynamics.kinematic_viscosity_m2_s`

- **default:** `1.3e-06` m2/s
- **summary:** Seawater kinematic viscosity for Reynolds numbers.
- **range:** >0
- **Case:** exposed (settable in musselflow_grammar.json)

### `hydrodynamics.water_density_kg_m3`

- **default:** `1020.0` kg/m3
- **summary:** Seawater density for drag and Reynolds numbers.
- **range:** >0
- **Case:** exposed (settable in musselflow_grammar.json)

## scenario

### `scenario.weights`

- **default:** `[1.0]` -
- **summary:** Relative weight per current scenario. In the strict case these are derived from step durations, not set directly.
- **range:** 0.0 .. +inf
- **Case:** derived (computed by the runtime)

### `scenario.duration_h`

- **default:** `[24.0]` h
- **summary:** Duration each current scenario is held while advancing oxygen and sediment stocks.
- **range:** strictly > 0
- **Case:** exposed (settable in musselflow_grammar.json)

### `scenario.repeat_count`

- **default:** `1` -
- **summary:** Number of times the ordered timeline repeats.
- **range:** 1 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

## sediment

### `sediment.settling_velocity_m_s`

- **default:** `0.006` m/s
- **summary:** Biodeposit settling velocity for the settling-distance screen.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Settling-distance screen from depth, transport speed, and settling velocity; first-order sediment box
- **envelope:** Biodeposit production and dispersion vary by mussel age, particle size, current, and site; not a seabed model
- **source:** Callier et al. (2006) — https://doi.org/10.3354/meps322129

### `sediment.in_domain_deposition_fraction`

- **default:** `-1.0` -
- **summary:** Fraction of biodeposits landing inside the domain. -1 = compute from a settling-distance screen instead of a seabed transport model.
- **range:** -1 (automatic) or 0..1
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Settling-distance screen from depth, transport speed, and settling velocity; first-order sediment box
- **envelope:** Biodeposit production and dispersion vary by mussel age, particle size, current, and site; not a seabed model
- **source:** Callier et al. (2006) — https://doi.org/10.3354/meps322129

### `sediment.oxygen_demand_kg_o2_per_kg_organic`

- **default:** `1.0` kg O2/kg
- **summary:** Oxygen demand per kilogram of deposited organic matter.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Settling-distance screen from depth, transport speed, and settling velocity; first-order sediment box
- **envelope:** Biodeposit production and dispersion vary by mussel age, particle size, current, and site; not a seabed model
- **source:** Callier et al. (2006) — https://doi.org/10.3354/meps322129

### `sediment.decay_per_day`

- **default:** `0.05` 1/day
- **summary:** First-order decay of the deposited organic stock.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Settling-distance screen from depth, transport speed, and settling velocity; first-order sediment box
- **envelope:** Biodeposit production and dispersion vary by mussel age, particle size, current, and site; not a seabed model
- **source:** Callier et al. (2006) — https://doi.org/10.3354/meps322129

### `sediment.resuspension_per_day`

- **default:** `0.0` 1/day
- **summary:** First-order resuspension of the deposited organic stock.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Settling-distance screen from depth, transport speed, and settling velocity; first-order sediment box
- **envelope:** Biodeposit production and dispersion vary by mussel age, particle size, current, and site; not a seabed model
- **source:** Callier et al. (2006) — https://doi.org/10.3354/meps322129

### `sediment.mortality_deposition_fraction`

- **default:** `1.0` -
- **summary:** Fraction of mortality organic loss deposited to the seabed.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Settling-distance screen from depth, transport speed, and settling velocity; first-order sediment box
- **envelope:** Biodeposit production and dispersion vary by mussel age, particle size, current, and site; not a seabed model
- **source:** Callier et al. (2006) — https://doi.org/10.3354/meps322129

### `sediment.initial_organic_stock_kg`

- **default:** `0.0` kg
- **summary:** Initial deposited organic stock at the start of a run.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** Settling-distance screen from depth, transport speed, and settling velocity; first-order sediment box
- **envelope:** Biodeposit production and dispersion vary by mussel age, particle size, current, and site; not a seabed model
- **source:** Callier et al. (2006) — https://doi.org/10.3354/meps322129

## harvest

### `harvest.fraction_per_year`

- **default:** `0.8` 1/year
- **summary:** Fraction of standing stock harvested per year.
- **range:** 0.0 .. 1.0
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** 13.7 kg N and 0.9 kg P per harvested wet tonne (Danish); 8.50 +/- 0.59 kg N and 0.95 +/- 0.07 kg P per live tonne (UK) -- composition must stay editable
- **envelope:** Management accounting prior; NOT predicted growth
- **source:** Danish mitigation-mussel values (Taylor et al. 2019); UK rope-culture survey (Mascorda Cabre et al. 2021) — https://doi.org/10.3389/fmars.2019.00698

### `harvest.turnovers_per_year`

- **default:** `1.0` 1/year
- **summary:** Harvest turnovers per year for the accounting scenario.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** 13.7 kg N and 0.9 kg P per harvested wet tonne (Danish); 8.50 +/- 0.59 kg N and 0.95 +/- 0.07 kg P per live tonne (UK) -- composition must stay editable
- **envelope:** Management accounting prior; NOT predicted growth
- **source:** Danish mitigation-mussel values (Taylor et al. 2019); UK rope-culture survey (Mascorda Cabre et al. 2021) — https://doi.org/10.3389/fmars.2019.00698

### `harvest.n_kg_per_t_wet`

- **default:** `13.7` kg N/t
- **summary:** Nitrogen removed per harvested wet tonne (editable composition).
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** 13.7 kg N and 0.9 kg P per harvested wet tonne (Danish); 8.50 +/- 0.59 kg N and 0.95 +/- 0.07 kg P per live tonne (UK) -- composition must stay editable
- **envelope:** Management accounting prior; NOT predicted growth
- **source:** Danish mitigation-mussel values (Taylor et al. 2019); UK rope-culture survey (Mascorda Cabre et al. 2021) — https://doi.org/10.3389/fmars.2019.00698

### `harvest.p_kg_per_t_wet`

- **default:** `0.9` kg P/t
- **summary:** Phosphorus removed per harvested wet tonne (editable composition).
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)
- **equation:** 13.7 kg N and 0.9 kg P per harvested wet tonne (Danish); 8.50 +/- 0.59 kg N and 0.95 +/- 0.07 kg P per live tonne (UK) -- composition must stay editable
- **envelope:** Management accounting prior; NOT predicted growth
- **source:** Danish mitigation-mussel values (Taylor et al. 2019); UK rope-culture survey (Mascorda Cabre et al. 2021) — https://doi.org/10.3389/fmars.2019.00698

## constraint

### `constraint.min_do_mg_l`

- **default:** `-1.0` mg/L
- **summary:** Minimum modeled dissolved oxygen. Negative disables the constraint; no legal threshold is asserted by the program.
- **range:** -1 disables, else >= 0
- **Case:** exposed (settable in musselflow_grammar.json)

### `constraint.min_probe_speed_m_s`

- **default:** `-1.0` m/s
- **summary:** Minimum required probe/flushing speed. Negative disables.
- **range:** -1 disables, else >= 0
- **Case:** exposed (settable in musselflow_grammar.json)

### `constraint.max_biodeposition_kg_m2_day`

- **default:** `-1.0` kg/m2/day
- **summary:** Maximum allowed biodeposition. Negative disables.
- **range:** -1 disables, else >= 0
- **Case:** exposed (settable in musselflow_grammar.json)

### `constraint.min_obstacle_clearance_m`

- **default:** `0.0` m
- **summary:** Minimum spacing between obstacles (collision constraint).
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

## objective

### `objective.target_chlorophyll_capture_g_day`

- **default:** `-1.0` g/day
- **summary:** Absolute chlorophyll-capture target. Negative uses the design's own theoretical clearance capacity for the extraction score.
- **range:** -1 = auto, else >= 0
- **Case:** exposed (settable in musselflow_grammar.json)

### `objective.target_biodeposition_kg_m2_day`

- **default:** `0.05` kg/m2/day
- **summary:** Reference biodeposition used to normalise the low-deposition score.
- **range:** strictly > 0
- **Case:** exposed (settable in musselflow_grammar.json)

### `weight.extraction`

- **default:** `0.3` -
- **summary:** Objective weight: particulate/chlorophyll extraction.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `weight.flushing`

- **default:** `0.2` -
- **summary:** Objective weight: flushing/water exchange.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `weight.food_delivery`

- **default:** `0.15` -
- **summary:** Objective weight: food delivery to obstacles.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `weight.uniformity`

- **default:** `0.1` -
- **summary:** Objective weight: uniformity of food delivery.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `weight.oxygen`

- **default:** `0.15` -
- **summary:** Objective weight: dissolved-oxygen safety proxy.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

### `weight.low_deposition`

- **default:** `0.1` -
- **summary:** Objective weight: low biodeposition proxy.
- **range:** 0.0 .. +inf
- **Case:** exposed (settable in musselflow_grammar.json)

## validation

### `validation.geometry_calibrated`

- **default:** `False` -
- **summary:** Governance flag: geometry descriptors validated against known shapes.
- **Case:** derived (computed by the runtime)

### `validation.hydrodynamics_calibrated`

- **default:** `False` -
- **summary:** Governance flag: wake/drag validated against current measurements.
- **Case:** derived (computed by the runtime)

### `validation.biology_calibrated`

- **default:** `False` -
- **summary:** Governance flag: clearance/capture validated against measurements.
- **Case:** derived (computed by the runtime)

### `validation.oxygen_calibrated`

- **default:** `False` -
- **summary:** Governance flag: oxygen box validated over tidal/low-flow sequences.
- **Case:** derived (computed by the runtime)

### `validation.sediment_calibrated`

- **default:** `False` -
- **summary:** Governance flag: deposition validated against traps and benthic flux.
- **Case:** derived (computed by the runtime)

### `validation.allow_extrapolation`

- **default:** `False` -
- **summary:** If true, use-envelope violations become soft warnings instead of a hard infeasible constraint. A deliberate, declared choice.
- **Case:** exposed (settable in musselflow_grammar.json)
