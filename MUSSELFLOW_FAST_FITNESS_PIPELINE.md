# MusselFlow fast-fitness pipeline

```mermaid
flowchart LR
    G["Gene pools / Galapagos genome<br/>centres, UV, rods, orientation"] --> GH["Existing Grasshopper geometry algorithm<br/>mussel socks / hydraulic meshes"]
    D["Optimization domain<br/>closed planar Curve"] --> C
    P["Optional diagnostic probes<br/>Point3d list"] --> C
    F["Tidal/current scenarios<br/>Vector3d list; length = m/s"] --> C
    J["Ecological grammar Panel<br/>musselflow_grammar.json"] --> C
    Q["qualityMode<br/>0 FAST / 1 REFINE"] --> C
    S["speedMode<br/>True = aggressive Galapagos search"] --> C

    subgraph SDK["MusselFlow Grasshopper SDK component"]
        C["musselflow_component_gh_sdk.py<br/>typed ports, Rhino units, validation"]
        B["musselflow_bio_optimizer_gh_sdk.py<br/>SPEED: 16-point descriptors<br/>FAST: 64-point descriptors<br/>REFINE: projected mesh area"]
        CP["musselflow_case_core.py<br/>strict JSON parsing and scenario mapping"]
        EV["musselflow_ecogrammar_core.py<br/>88 parameter checks and broadcasting"]
        N["musselflow_bio_optimizer_core.py<br/>deterministic reduced-order physics"]
        C --> B
        C --> CP
        CP --> EV
        B --> N
        EV --> N
        F --> N
    end

    N --> A["Constraint-dominated aggregate"]
    A --> O1["Fitness<br/>single Galapagos number"]
    A --> O2["Feasible + Constraints"]
    A --> O3["Objectives"]
    A --> O4["Result JSON<br/>future recorder / preview"]
    A --> O5["Status + Report"]
    O1 --> GA["Galapagos selection / mutation"]
    GA --> G
```

## Runtime path

1. Grasshopper generates one candidate geometry set from the current genome.
2. The SDK validates domain, geometry, probes, flow vectors, and Rhino units.
3. The unchanged ecological grammar is parsed and compiled once, then cached.
4. `speedMode = True` uses at most 16 geometry points and 64 domain samples
   for broad Galapagos exploration while retaining all ecological equations.
   With speed mode off, `qualityMode = 0` reduces Rhino geometry to six values:
   `x, y, major, minor, yaw, height`. `qualityMode = 1` retains that stable
   envelope and additionally measures direction-specific projected frontal
   area from a coarse analysis mesh.
5. The NumPy core evaluates wake, food, filtration, oxygen, sediment, objectives,
   and hard constraints.
6. A feasible design receives fitness in `[0.25, 1]`; an infeasible design
   remains below `0.25`, so ecological reward cannot overpower safety limits.
7. Galapagos changes the upstream genes and repeats the pipeline.

The detailed `Result` output is deliberately separate from viewport
visualization. It can later feed the wild-camera recorder and time-lapse
preview without placing animation work inside every Galapagos evaluation.
