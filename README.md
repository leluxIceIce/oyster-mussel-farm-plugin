# MusselFlow

MusselFlow is a Rhino 8 / Grasshopper research toolkit for computational exploration of large-scale mussel farm design strategies. It has a solver build up as a deterministic reduced-order screening of current exposure, food delivery, filtration, dissolved-oxygen risk, deposition, and ecological constraints. Copernicus components provide regional site forcing and spatial field exploration.

## Current status

This repository contains the canonical working Python sources selected from the development workspace.
MusselFlow is a comparative screening model. It is not CFD, a calibrated carrying-capacity model, verified nutrient or carbon removal, a legal assessment, or a substitute for site monitoring.

## Pipeline

```mermaid
flowchart LR
  G[“Mesh geometry<br/>centres, orientation, rod dimensions"] --> O["MusselFlow optimizer"]
  C["Ecological grammar<br/>coefficients, objectives, constraints"] --> O
  X["Search context"] --> B["Copernicus browser"]
  B --> R["Fetch request"]
  R --> S["Site data router"]
  S -->|"flow vectors + simulation case"| O
  S --> F["Site field"]
  F --> M["Multivariate exploration"]
  O --> Q["Fitness, feasibility, objectives, result"]
  Q --> GA["Galapagos"]
  O --> V["Detailed flow visualizer"]
```

## Canonical components

- `musselflow_component_gh_sdk.py` — ecological optimizer and Galapagos fitness input.
- `musselflow_gh_sdk.py` — detailed plan-view flow and food visualization.
- `copernicus_search_context_gh_sdk.py` — validated WGS84/time/depth/variable query.
- `copernicus_data_browser_gh_sdk.py` — internal API-called Eto catalogue browser.
- `musselflow_site_data_gh_sdk.py` — Sample regional Copernicus satellite data and compile the field data.
- `musselflow_site_field_gh_sdk.py` — spatial environmental field generation.
- `musselflow_multivariate_gh_sdk.py` — correlation and PCA exploration.
- `musselflow_cluster_genome_gh_sdk.py`, `musselflow_family_decoder_gh_sdk.py`, `musselflow_geometry_gh_sdk.py`, and `kangaroo_frame_loft_gh_sdk.py` — parametric farm geometry tools.

The optimizer requires the colocated runtime sidecars listed in `MUSSELFLOW_SETUP.md`. The flat source layout is deliberate for the current Grasshopper workflow. A compiled cross-platform plug-in can package these modules later.

## Data contracts

- `musselflow_grammar.json` is the single editable ecological model.
- `musselflow_case_schema.json` is its machine-readable validation contract.
- `MUSSELFLOW_JSON_CONTRACTS.md` explains Search Context, Fetch Request, Site Data, Simulation Case, and optimizer Result.

## Tests

From this directory:

```bash
python3 -m unittest \
  test_copernicus_data_browser \
  test_musselflow_bio_optimizer \
  test_musselflow_case_core \
  test_musselflow_grammar_book
```

The initial canonical import passes 85 relevant regression tests. Native Eto-window behavior still requires a Rhino 8 smoke test.

## Machine learning

No reinforcement-learning controller or unvalidated surrogate is active in the fitness path. `musselflow_surrogate_core.py` is a guarded offline residual-model utility reserved for future calibration against CFD or field observations.
