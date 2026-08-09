# MusselFlow feeding and energy screening model

This document describes the deterministic biological layer used by
`musselflow_bio_optimizer_core.py`. It is a fast, transparent **DEB-lite
screening model**, not a calibrated Dynamic Energy Budget state model, growth
forecast, carrying-capacity assessment, CFD model, or legal evidence.

## Runtime path

For each flow scenario and obstacle, the solver calculates:

1. local speed and food fractions from the reduced wake/plume model;
2. maximum allometric clearance and respiration from mussel dry mass;
3. low-food valve activity and Holling-like ingestion saturation;
4. oxygen-saturation and aggregation-dependent current activity;
5. two-class particle-size retention;
6. TSM-dependent pseudofaeces and food-quality-dependent assimilation;
7. separate assimilated matter, faeces, pseudofaeces and biodeposition;
8. carbon accounting and an energy-balance scope-for-growth proxy.

The organic mass balance is enforced as:

`filtered organic = assimilated organic + faeces + pseudofaeces`

Potential growth is reported only as:

`max(assimilated energy - respiration energy, 0) / tissue energy density`

It does not update individual size, density, mortality, reproduction, harvest,
or self-thinning through time.

## Sources and interpretation

- Maar et al. (2023), DOI `10.1016/j.scitotenv.2023.164168`, motivates the
  process structure: particle classes, saturating ingestion, food quality,
  assimilation, respiration and biodeposit pathways.
- Nielsen & Vismann (2014), DOI `10.2983/035.033.0214`, motivates the editable
  aggregation-dependent current response.
- Kamermans & Saurel (2022), DOI `10.1051/alr/2022001`, motivates explicit
  oxygen-food-temperature interaction. The implemented oxygen ramp is an
  editable screening approximation, not a fitted universal response.
- Strohmeier et al. (2012), DOI `10.1016/j.jembe.2011.11.006`, motivates
  particle-size-dependent retention.

All numerical values in the grammar are **screening priors**. Literature can
justify model structure and plausible ranges; it does not make a site-specific
value true. Field or laboratory measurements must replace the defaults before
scientific interpretation.

## Calibration priorities

The first measurements to calibrate are current speed/direction distributions,
sock porosity and frontal area, mussel dry-mass distribution and aggregation,
clearance versus current, particle-size spectra, chlorophyll and TSM, organic
fraction and food quality, dissolved oxygen through low-flow periods, and
biodeposition/settling observations.

The one-box oxygen and sediment calculations cannot resolve stratification,
near-bed boundary layers, resuspension transport, or spatial hypoxia. Selected
layouts must therefore be checked with time-resolved field data and a suitable
3D hydrodynamic-biogeochemical model.

## Grasshopper outputs

The seven existing component outputs are unchanged. Physical rates are nested
inside `Result.screening_metrics` and per-obstacle values inside
`Result.obstacles`. `Fitness` and the six objective values remain normalized
and dimensionless because Galapagos requires one comparable scalar and rates
with unlike units must not be added directly.

