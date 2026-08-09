# MusselFlow JSON contracts

The visible Grasshopper names describe the state of the information. Port
positions remain unchanged; this naming convention does not alter existing
wires or numerical behavior.

## Canonical sequence

```text
SearchContextJson
  -> Copernicus Data Browser
FetchRequestJson
  -> MusselFlow Site Data (+ fetch Button + optional BaseModelJson)
SiteDataJson
  -> Site Field / Multivariate Explorer
SimulationCaseJson
  -> MusselFlow Ecological Optimizer
```

## Meaning of each document

### SearchContextJson

User intent before choosing a data product: WGS84 location, time range,
depths, requested variables, optional layer/pixel preference and keywords.
It contains no catalogue selection and no environmental values.

### FetchRequestJson

SearchContextJson plus the product accepted in the Copernicus Data Browser.
Its `data_state` is `not_fetched`. It is an executable request description,
not evidence that a network request succeeded.

### SiteDataJson

Actual Copernicus responses assembled by MusselFlow Site Data: timestamps,
depth profiles, physical values and units, layer URLs, requested and sampled
coordinates, missing fields, HTTP provenance and scientific limitations. Its
`data_state` is `sampled`.

### BaseModelJson

The site-independent MusselFlow ecological model: mussel species and biomass,
clearance/respiration relations, food partition, structure, oxygen/sediment
terms, objectives, constraints, evidence and validation flags. It does not
represent one fetched day at one coordinate.

### SimulationCaseJson

BaseModelJson patched with the SiteDataJson forcing timeline. It is the
executable optimizer input. It intentionally contains only the environmental
subset needed by the reduced-order solver; the complete profiles and request
provenance remain in SiteDataJson. Flow-vector values stay on the parallel
typed `FlowVectors` wire, while timeline steps reference their indices.

### FieldDataJson and AnalysisDataJson

FieldDataJson contains one spatially sampled variable over a Rhino domain.
AnalysisDataJson contains correlation/PCA results derived from SiteDataJson or
aligned FieldDataJson documents. Neither replaces the source SiteDataJson.

## Why SiteDataJson and SimulationCaseJson stay separate

Combining them would mix source evidence with model assumptions and make it
hard to audit which values came from Copernicus and which came from the mussel
grammar. The two-document boundary also lets one downloaded SiteDataJson be
tested against several species, stocking or objective configurations without
fetching the ocean data again.
