# r: numpy
"""
MusselFlow Fitness Assembly
===========================

REDUCED ORDER: Xplore how mussel socks (nets and frames) behave as one coordinated hydraulic field under multiple currents.

Every farm layout becomes a test of ecological swarm behaviour: where water is
redirected, where fresh food reaches the mussels, where filtration accumulates,
and where weak flushing raises deposition or oxygen risk. Rhino geometry and
the ecological grammar are read together across changing current scenarios,
then distilled into one Galapagos fitness and a detailed result. Flow and food
visualizations remain in their dedicated component, keeping the evolutionary
search light enough to explore thousands of arrangements.

Name: MusselFlowOptimizer
Updated: 260808
Author: Felix Berger
Copyright: Apache License 2.0 

Scientific scope
----------------
The core is intentionally reduced-order: fast enough to search a large design
space while remaining clearly labelled ``UNVALIDATED_SCREENING``. Its outputs
are comparative estimates, not CFD, carrying-capacity proof, regulatory
assessment, predicted annual growth, verified nutrient mitigation, or carbon
removal. No trained machine-learning or reinforcement-learning model is active.

Required companion files
------------------------
Keep these internal modules together with this component in the
``Miesmuschel`` folder, or beside the saved Grasshopper definition:

* ``musselflow_case_core.py``
* ``musselflow_ecogrammar_core.py``
* ``musselflow_bio_optimizer_core.py``
* ``musselflow_bio_optimizer_gh_sdk.py``

Inputs
------
run : item / bool
    Evaluate the current candidate. False returns a stable WAITING result.
obstacles : list / Rhino.Geometry.GeometryBase
    Mussel socks, nets, framed surfaces, meshes, Breps, curves, extrusions, or
    SubDs. All list items form one farm evaluation.
domain : item / Rhino.Geometry.Curve
    One closed planar optimization boundary. Rhino document units are converted
    to metres.
probes : list / Rhino.Geometry.Point3d
    Optional diagnostic measurement points. Invalid, duplicate, and
    out-of-domain probes are reported and ignored. With no accepted probes, the
    numerical core uses obstacle centres as its explicit fallback.
flowVectors : list / Rhino.Geometry.Vector3d
    Uniform current scenarios. Planar direction is flow direction and planar
    vector length is speed in m/s. A normal component is reported and discarded.
SimulationCaseJson : item / str
    Complete ecological grammar JSON, normally pasted from
    ``musselflow_grammar.json`` into a Grasshopper Panel.
qualityMode : item / int
    ``0`` FAST for Galapagos; ``1`` REFINE for direction-specific projected
    mesh/Brep frontal area and elite-design re-ranking.
speedMode : item / bool
    True enables the aggressively sampled Galapagos search path. It overrides
    REFINE geometry but retains the complete ecological solver and constraints.

Outputs — seven ports in this order
-----------------------------------
Fitness : item / float
    Scalar in [0,1] for Galapagos maximization.
Feasible : item / bool
    True only when every active hard constraint passes.
Objectives : list / str
    Named normalized objective values.
Constraints : list / str
    Named constraint margins; non-negative values pass.
Result : item / str
    Canonical structured result JSON for downstream Inspector/Preview nodes.
Status : item / str
    WAITING, INVALID_GEOMETRY, INVALID_FLOW, INVALID_CASE, UNSUPPORTED_CASE,
    SOLVER_ERROR, OUT_OF_ENVELOPE, INFEASIBLE, or UNVALIDATED_SCREENING.
Report : list / str
    Timing, units, scenario mapping, assumptions, warnings, and limitations.

Current forcing limitation
--------------------------
The parser validates both timelines and ensembles. This component currently
evaluates ordered timelines only. Ensemble execution already exists in the
numerical core and is the next Grasshopper integration step.
"""

import importlib
import math
import os
import sys
import time

import Grasshopper
import numpy as np
import Rhino
import System


COMPONENT_METADATA = {
    "name": "MusselFlow Ecological Optimizer",
    "nickname": "MusselFlow",
    "description": (
        "Fast reduced-order screening of mussel-farm geometry under a strict "
        "ecological case. Returns one Galapagos fitness and structured "
        "diagnostics; it is not CFD or regulatory evidence."),
}

COMPONENT_BUILD = "2026-08-08b"
OUTPUT_DECIMALS = 6
OUTPUT_SIGNIFICANT_DIGITS = 6

# The ecological case and its compiled numerical configuration stay constant
# throughout a Galapagos run. Cache that immutable work by exact source text,
# counts, and parser source stamp; changing any of them causes a safe miss.
_CASE_COMPILE_CACHE = {}

INPUT_METADATA = (
    ("run", "run", "Evaluate the current candidate. False returns WAITING."),
    ("obstacles", "obstacles", "Mussel geometry supplied as one Grasshopper list."),
    ("domain", "domain", "One closed planar boundary; model units become metres."),
    ("probes", "probes", "Optional Point3d list; obstacle centres are the fallback."),
    ("flowVectors", "flowVectors", "Vector3d list; planar length is speed in m/s."),
    ("SimulationCaseJson", "SimulationCase", "Executable ecological model with site forcing from Site Data."),
    ("qualityMode", "qualityMode", "0 = FAST Galapagos; 1 = REFINE elite designs."),
    ("speedMode", "speedMode", "True = aggressive SPEED search; full biology remains active."),
)

OUTPUT_METADATA = (
    ("Fitness", "Fitness", "Scalar fitness in [0,1] for Galapagos."),
    ("Feasible", "Feasible", "True when every active hard constraint passes."),
    ("Objectives", "Objectives", "Named normalized objective values."),
    ("Constraints", "Constraints", "Named margins; non-negative values pass."),
    ("Result", "Result", "Canonical result JSON for downstream components."),
    ("Status", "Status", "Runtime and scientific-validity status."),
    ("Report", "Report", "Timing, units, assumptions, warnings, and limits."),
)


def apply_component_metadata(component):
    """Assign Grasshopper labels and hover help without changing port topology.

    Port creation/removal is intentionally excluded: mutating parameter topology
    during a Grasshopper solution can invalidate existing wires. The typed SDK
    signature owns input creation; seven output ports are configured once and
    are best preserved in a Grasshopper User Object.
    """
    if component is None:
        return
    component.Name = COMPONENT_METADATA["name"]
    component.NickName = COMPONENT_METADATA["nickname"]
    component.Description = COMPONENT_METADATA["description"]

    for index, (name, nickname, description) in enumerate(INPUT_METADATA):
        if index >= component.Params.Input.Count:
            break
        parameter = component.Params.Input[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description
        if index in (1, 3, 4):
            # Keep farm geometry, probes, and all current scenarios in one
            # RunScript call even when this SDK script is pasted over an older
            # component whose sockets still carry Item access.
            parameter.Access = Grasshopper.Kernel.GH_ParamAccess.list

    for index, (name, nickname, description) in enumerate(OUTPUT_METADATA):
        if index >= component.Params.Output.Count:
            break
        parameter = component.Params.Output[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description


def empty_outputs(status, message):
    """Return the stable seven-output contract for waiting/error states."""
    return (
        0.0,
        False,
        [],
        [],
        "{}",
        status,
        ["BUILD | "+COMPONENT_BUILD, message],
    )


def tree_records(value):
    """Read a DataTree without discarding each item's original GH path."""
    if value is None:
        return []
    try:
        branches = value.Branches
        paths = value.Paths
    except Exception:
        try:
            items = list(value)
        except TypeError:
            items = [value]
        records = []
        def append_item(item):
            if isinstance(item, (list, tuple)):
                for nested in item:
                    append_item(nested)
                return
            try:
                if isinstance(item, System.Collections.IList):
                    for nested in item:
                        append_item(nested)
                    return
            except Exception:
                pass
            if item is not None:
                records.append(("{0}", 0, len(records), item))
        for item in items:
            append_item(item)
        return records

    records = []
    def append_item(path, branch_index, item):
        if isinstance(item, (list, tuple)):
            for nested in item:
                append_item(path, branch_index, nested)
            return
        try:
            if isinstance(item, System.Collections.IList):
                for nested in item:
                    append_item(path, branch_index, nested)
                return
        except Exception:
            pass
        if item is not None:
            records.append((path, branch_index, len(records), item))
    for branch_index, branch in enumerate(branches):
        path = str(paths[branch_index])
        for item in branch:
            append_item(path, branch_index, item)
    return records


def resolve_rhino_value(value):
    """Unwrap GH goo and document GUIDs into RhinoCommon values."""
    if value is None:
        return None
    if hasattr(value, "ScriptVariable"):
        try:
            value = value.ScriptVariable()
        except Exception:
            pass
    if hasattr(value, "Value"):
        try:
            value = value.Value
        except Exception:
            pass
    if isinstance(value, System.Guid):
        document = Rhino.RhinoDoc.ActiveDoc
        rhino_object = (
            document.Objects.FindId(value) if document is not None else None)
        return rhino_object.Geometry if rhino_object is not None else None
    return value


def curve_plane(curve, tolerance):
    """Return a curve plane across RhinoCommon overload variations."""
    try:
        success, plane = curve.TryGetPlane(tolerance)
        if success:
            return plane
    except Exception:
        try:
            success, plane = curve.TryGetPlane()
            if success:
                return plane
        except Exception:
            pass
    return None


def is_finite_point(point):
    return all(math.isfinite(value) for value in (point.X, point.Y, point.Z))


def is_finite_vector(vector):
    return all(math.isfinite(value) for value in (vector.X, vector.Y, vector.Z))


def resolve_point(value):
    """Resolve Point3d, Rhino point geometry, or GH point goo."""
    value = resolve_rhino_value(value)
    if isinstance(value, Rhino.Geometry.Point3d):
        return Rhino.Geometry.Point3d(value)
    if isinstance(value, Rhino.Geometry.Point):
        return Rhino.Geometry.Point3d(value.Location)
    try:
        return Rhino.Geometry.Point3d(value)
    except Exception:
        return None


def resolve_vector(value):
    """Resolve Vector3d or GH vector goo without applying model-unit scaling."""
    value = resolve_rhino_value(value)
    if isinstance(value, Rhino.Geometry.Vector3d):
        return Rhino.Geometry.Vector3d(value)
    try:
        return Rhino.Geometry.Vector3d(value)
    except Exception:
        return None


def supported_obstacle(geometry):
    """Accept the geometry families currently intended for sock/net obstacles."""
    accepted = [
        Rhino.Geometry.Curve,
        Rhino.Geometry.Mesh,
        Rhino.Geometry.Brep,
        Rhino.Geometry.Surface,
        Rhino.Geometry.Extrusion,
    ]
    if hasattr(Rhino.Geometry, "SubD"):
        accepted.append(Rhino.Geometry.SubD)
    return isinstance(geometry, tuple(accepted))


def model_scale_and_tolerance():
    """Return active-document metres/model-unit, tolerance, and warning."""
    document = Rhino.RhinoDoc.ActiveDoc
    if document is None:
        return 1.0, 0.001, (
            "No active Rhino document; geometry is being treated as metres.")
    tolerance = max(float(document.ModelAbsoluteTolerance), 1e-12)
    try:
        scale = Rhino.RhinoMath.UnitScale(
            document.ModelUnitSystem, Rhino.UnitSystem.Meters)
    except Exception:
        return 1.0, tolerance, (
            "Could not read Rhino model units; geometry is being treated "
            "as metres.")
    if not math.isfinite(scale) or scale <= 0:
        return 1.0, tolerance, (
            "Rhino model units are unset; geometry is being treated as metres.")
    return float(scale), tolerance, None


def validate_domain(value, tolerance):
    """Return ``(curve, plane, area_model_units)`` or an error message."""
    curve = resolve_rhino_value(value)
    if not isinstance(curve, Rhino.Geometry.Curve):
        return None, None, None, (
            "domain must resolve to one Rhino.Geometry.Curve.")
    if not curve.IsClosed:
        return None, None, None, "domain Curve must be closed."
    plane = curve_plane(curve, tolerance)
    if plane is None:
        return None, None, None, "domain Curve must be planar."
    try:
        properties = Rhino.Geometry.AreaMassProperties.Compute(curve)
        area = float(properties.Area) if properties is not None else 0.0
    except Exception:
        area = 0.0
    if not math.isfinite(area) or area <= tolerance*tolerance:
        return None, None, None, (
            "domain Curve must enclose measurable positive area.")
    return curve, plane, area, None


def project_folders(component):
    """Return ordered, unique folders that may contain MusselFlow sidecars."""
    candidate_folders = []
    try:
        document = component.OnPingDocument()
        definition_path = document.FilePath if document is not None else ""
        if definition_path:
            candidate_folders.append(os.path.dirname(definition_path))
    except Exception:
        pass
    candidate_folders.extend([
        os.getcwd(),
        "/Users/lelux/Desktop/Datajunk_Graphx_Entre/Python",
    ])
    folders = []
    for folder in candidate_folders:
        if not folder:
            continue
        folder = os.path.abspath(folder)
        if folder not in folders:
            folders.append(folder)
    return folders


def load_project_module(component, module_name):
    """Import a sidecar and reload it only when its source file changed.

    Rhino keeps imported Python modules in ``sys.modules`` across Grasshopper
    solutions. Without the source stamp, replacing a sidecar can leave an old
    parser or solver active until Rhino restarts.
    """
    for folder in project_folders(component):
        module_path = os.path.join(folder, module_name+".py")
        if not os.path.isfile(module_path):
            continue
        if folder not in sys.path:
            sys.path.insert(0, folder)
        importlib.invalidate_caches()
        try:
            source_path = os.path.abspath(module_path)
            source_stamp = os.stat(source_path).st_mtime_ns
            module = sys.modules.get(module_name)
            loaded_path = os.path.abspath(
                getattr(module, "__file__", "")) if module is not None else ""
            if module is not None and loaded_path != source_path:
                del sys.modules[module_name]
                module = None
            if module is None:
                module = importlib.import_module(module_name)
            elif getattr(
                    module, "__musselflow_source_stamp__", None) != source_stamp:
                module = importlib.reload(module)
            module.__musselflow_source_stamp__ = source_stamp
            return module
        except (ImportError, OSError) as exception:
            last_exception = exception
    try:
        module = importlib.import_module(module_name)
        return module
    except ImportError as exception:
        last_exception = exception
    raise ImportError(
        "Could not load %s.py. Put all required MusselFlow files beside the "
        "saved Grasshopper definition. DETAILS | %s"
        % (module_name, last_exception))


def load_case_core(component):
    """Load the Rhino-independent case parser."""
    module = load_project_module(component, "musselflow_case_core")
    if (getattr(module, "SCHEMA_VERSION", None) != "1.0.0" or
            not hasattr(module, "compile_timeline")):
        raise ImportError(
            "OUTDATED SIDECAR | %s | Replace this file with the current "
            "musselflow_case_core.py."
            % getattr(module, "__file__", "unknown path"))
    return module


def load_runtime_modules(component):
    """Load the case parser, numerical core, and Rhino descriptor bridge."""
    # Dependency order matters when source files were edited between solutions.
    load_project_module(component, "musselflow_ecogrammar_core")
    case_core = load_case_core(component)
    module_names = (
        "musselflow_bio_optimizer_core",
        "musselflow_bio_optimizer_gh_sdk",
    )
    modules = []
    for module_name in module_names:
        modules.append(load_project_module(component, module_name))
    return case_core, modules[0], modules[1]


def compile_case_cached(
        case_core, source, obstacle_count, flow_count):
    """Parse and compile an unchanged ecological case only once per run."""
    source = str(source)
    cache_key = (
        source,
        int(obstacle_count),
        int(flow_count),
        getattr(case_core, "__musselflow_source_stamp__", None),
    )
    cached = _CASE_COMPILE_CACHE.get(cache_key)
    if cached is not None:
        return cached, True

    case, warnings = case_core.parse_case(
        source,
        flow_count=flow_count,
        obstacle_count=obstacle_count)
    digest = case_core.case_hash(case)
    config, flow_indices, bridge_warnings = case_core.compile_timeline(
        case,
        obstacle_count=obstacle_count,
        flow_count=flow_count)
    compiled = (
        case, digest, config, flow_indices,
        tuple(warnings), tuple(bridge_warnings))
    # One Grasshopper component normally uses one case. A small bound prevents
    # repeated Panel edits from retaining arbitrarily many large dictionaries.
    if len(_CASE_COMPILE_CACHE) >= 4:
        _CASE_COMPILE_CACHE.clear()
    _CASE_COMPILE_CACHE[cache_key] = compiled
    return compiled, False


def compact_number(value):
    """Keep six meaningful digits without collapsing small values to zero."""
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Result contains a non-finite number.")
    if number == 0.0:
        return 0.0
    decimal_places = (
        OUTPUT_SIGNIFICANT_DIGITS-1-
        int(math.floor(math.log10(abs(number)))))
    rounded = round(number, decimal_places)
    return 0.0 if rounded == 0.0 else rounded


def display_number(value):
    """Format public text compactly without displaying nonzero values as zero."""
    number = compact_number(value)
    if number != 0.0 and abs(number) < 10.0**(-OUTPUT_DECIMALS):
        return "%.6g" % number
    return "%.6f" % number


def json_safe(value):
    """Convert results to finite JSON values rounded for public output."""
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, dict):
        return {
            str(key): json_safe(item)
            for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return value
    if isinstance(value, (int, float)):
        number = float(value) if isinstance(value, float) else value
        if not math.isfinite(float(number)):
            raise ValueError("Result contains a non-finite number.")
        if isinstance(number, float):
            return compact_number(number)
        return number
    return str(value)


def named_values(mapping):
    """Format a stable Grasshopper-friendly list without losing names."""
    return [
        "%s = %s" % (name, display_number(mapping[name]))
        for name in sorted(mapping)]


def make_result_document(
        case, case_digest, result, scenario_items, flow_indices,
        domain_area_m2, obstacle_records, flow_records, total_ms,
        geometry_ms, solver_ms, warnings, fidelity, descriptor_name):
    """Build the canonical downstream Result JSON structure."""
    scenario_documents = []
    for scenario_index, scenario in enumerate(scenario_items):
        source_flow_index = flow_indices[scenario_index]
        scenario_documents.append({
            "id": scenario["id"],
            "order": scenario.get("order"),
            "flow_vector_index": source_flow_index,
            "flow_path": flow_records[source_flow_index][0],
            "duration_h": scenario.get("duration_h"),
            "speed_m_s": result["scenario_speeds_m_s"][scenario_index],
            "deposition_fraction":
                result["scenario_deposition_fraction"][scenario_index],
            "twine_reynolds":
                result["scenario_twine_reynolds"][scenario_index],
            "mass_balance_scale":
                result["scenario_mass_balance_scale"][scenario_index],
            "mean_frontal_fill":
                float(np.mean(
                    result["scenario_frontal_fill"][scenario_index])),
        })

    obstacle_ids = [
        "%s[%d]" % (path, item_index)
        for path, branch_index, item_index, geometry in obstacle_records]
    status = (
        "OUT_OF_ENVELOPE"
        if result["constraint_margins"].get(
            "validation_envelope", 0.0) < 0.0
        else "INFEASIBLE" if not result["feasible"]
        else "UNVALIDATED_SCREENING")

    return {
        "schema_version": "1.0.0",
        "run": {
            "status": status,
            "model_status": result["model_status"],
            "feasible": result["feasible"],
            "fitness": result["fitness"],
            "raw_objective": result["raw_objective"],
            "timing_ms": {
                "geometry": geometry_ms,
                "solver": solver_ms,
                "total": total_ms,
            },
        },
        "case": {
            "case_id": case["case_id"],
            "case_hash": case_digest,
            "forcing_mode": case["forcing"]["mode"],
            "domain_area_m2": domain_area_m2,
        },
        "objectives": result["objectives"],
        "objective_weights": result["objective_weights"],
        "constraints": result["constraint_margins"],
        "scenarios": scenario_documents,
        "obstacles": {
            "ids": obstacle_ids,
            "speed_m_s": result["obstacle_speed_m_s"],
            "speed_ratio": result["obstacle_speed_ratio"],
            "food_fraction": result["obstacle_food_fraction"],
            "removal_fraction": result["obstacle_removal_fraction"],
            "chlorophyll_capture_g_day":
                result["chlorophyll_capture_g_day_by_obstacle"],
            "active_clearance_l_h":
                result["active_clearance_l_h_by_obstacle"],
            "phyto_retention": result["phyto_retention_by_obstacle"],
            "particulate_retention":
                result["particulate_retention_by_obstacle"],
            "ingested_organic_kg_day":
                result["ingested_organic_kg_day_by_obstacle"],
            "assimilated_organic_kg_day":
                result["assimilated_organic_kg_day_by_obstacle"],
            "pseudofaeces_organic_kg_day":
                result["pseudofaeces_organic_kg_day_by_obstacle"],
            "faeces_organic_kg_day":
                result["faeces_organic_kg_day_by_obstacle"],
            "pseudofaeces_fraction":
                result["pseudofaeces_fraction_by_obstacle"],
            "assimilation_efficiency":
                result["assimilation_efficiency_by_obstacle"],
            "oxygen_activity": result["oxygen_activity_by_obstacle"],
            "current_activity": result["current_activity_by_obstacle"],
            "scope_for_growth_kj_day":
                result["scope_for_growth_kj_day_by_obstacle"],
            "potential_growth_g_dw_day":
                result["potential_growth_g_dw_day_by_obstacle"],
        },
        "probes": {
            "source": result["probe_source"],
            "points_xy_m": result["probe_points"],
            "speed_m_s": result["probe_speed_m_s"],
            "speed_ratio": result["probe_speed_ratio"],
            "food_fraction": result["probe_food_fraction"],
            "deflection_deg": result["probe_deflection_deg"],
            "scenario_speed_m_s": result["scenario_probe_speed_m_s"],
            "scenario_food_fraction":
                result["scenario_probe_food_fraction"],
        },
        "screening_metrics": {
            "effective_cleared_m3_day":
                result["effective_cleared_m3_day"],
            "chlorophyll_capture_g_day":
                result["chlorophyll_capture_g_day"],
            "particulate_capture_kg_day":
                result["particulate_capture_kg_day"],
            "assimilated_organic_kg_day":
                result["assimilated_organic_kg_day"],
            "ingested_organic_kg_day":
                result["ingested_organic_kg_day"],
            "pseudofaeces_organic_kg_day":
                result["pseudofaeces_organic_kg_day"],
            "faeces_organic_kg_day": result["faeces_organic_kg_day"],
            "biodeposit_organic_kg_day":
                result["biodeposit_organic_kg_day"],
            "filtered_carbon_kg_day": result["filtered_carbon_kg_day"],
            "assimilated_carbon_kg_day":
                result["assimilated_carbon_kg_day"],
            "phytoplankton_carbon_capture_kg_day":
                result["phytoplankton_carbon_capture_kg_day"],
            "assimilated_energy_kj_day":
                result["assimilated_energy_kj_day"],
            "respiration_energy_kj_day":
                result["respiration_energy_kj_day"],
            "scope_for_growth_kj_day":
                result["scope_for_growth_kj_day"],
            "potential_growth_g_dw_day":
                result["potential_growth_g_dw_day"],
            "minimum_do_mg_l": result["oxygen"]["minimum_mg_l"],
            "final_do_mg_l": result["oxygen"]["final_mg_l"],
            "mussel_respiration_kg_o2_day":
                result["mussel_respiration_kg_o2_day"],
            "ammonia_excretion_kg_n_day":
                result["ammonia_excretion_kg_n_day"],
        },
        "prescribed_harvest_accounting_not_predicted_growth": {
            "wet_t_year": result["harvested_wet_t_year"],
            "n_kg_year": result["harvest_n_kg_year"],
            "p_kg_year": result["harvest_p_kg_year"],
        },
        "uncertainty": {
            "available": False,
            "outside_validation_envelope": status == "OUT_OF_ENVELOPE",
            "surrogate_active": False,
            "reinforcement_learning_active": False,
        },
        "provenance": {
            "case_hash": case_digest,
            "fidelity": fidelity,
            "descriptor": descriptor_name,
            "numerical_core": "musselflow_bio_optimizer_core",
            "scientific_status": "UNVALIDATED_SCREENING",
        },
        "warnings": list(warnings)+list(result["warnings"]),
    }


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def BeforeRunScript(self):
        """Keep component and parameter metadata synchronized before each solve."""
        apply_component_metadata(getattr(self, "Component", None))

    def RunScript(
            self,
            run: bool,
            obstacles: list[Rhino.Geometry.GeometryBase],
            domain: Rhino.Geometry.Curve,
            probes: list[Rhino.Geometry.Point3d],
            flowVectors: list[Rhino.Geometry.Vector3d],
            caseJson: str,
            qualityMode: int,
            speedMode: bool):
        """
        Inputs:
            run: Evaluate the current candidate {item,bool}
            obstacles: Mussel sock/net/frame geometry for one farm
                {list,GeometryBase}
            domain: Closed planar optimization boundary {item,Curve}
            probes: Optional diagnostic measurement locations
                {list,Point3d}
            flowVectors: Current scenarios; planar direction is flow direction
                and vector length is speed in m/s {list,Vector3d}
            SimulationCaseJson: Strict MusselFlow ecological model plus site forcing {item,str}
            qualityMode: 0 FAST or 1 REFINE {item,int}
            speedMode: True for aggressively sampled Galapagos search;
                overrides REFINE geometry but retains full biology {item,bool}

        Outputs, in this exact order:
            Fitness: Scalar value for Galapagos maximization {item,float}
            Feasible: True only when every active hard constraint passes
                {item,bool}
            Objectives: Named normalized objective values {list,str}
            Constraints: Named hard-constraint margins {list,str}
            Result: Canonical structured result JSON {item,str}
            Status: Runtime/model validity status {item,str}
            Report: Timing, assumptions, errors, and warnings {list,str}
        """
        if not run:
            return empty_outputs("WAITING", "Set run to True to evaluate.")
        total_start = time.perf_counter()
        try:
            quality_mode = 0 if qualityMode is None else int(qualityMode)
        except Exception:
            return empty_outputs(
                "INVALID_CASE", "qualityMode must be integer 0 or 1.")
        if quality_mode not in (0, 1):
            return empty_outputs(
                "INVALID_CASE", "qualityMode must be 0 FAST or 1 REFINE.")
        speed_mode = False if speedMode is None else bool(speedMode)
        if speed_mode:
            fidelity = "SPEED"
            descriptor_name = "PCA_OBB_6_SPEED_16"
        else:
            fidelity = "FAST" if quality_mode == 0 else "REFINE"
            descriptor_name = (
                "PCA_OBB_6_FAST"
                if quality_mode == 0 else
                "PCA_OBB_6_PLUS_PROJECTED_FRONTAL_AREA")

        metres_per_model_unit, tolerance, unit_warning = (
            model_scale_and_tolerance())
        domain_curve, plane, area_model, domain_error = validate_domain(
            domain, tolerance)
        if domain_error:
            return empty_outputs("INVALID_GEOMETRY", domain_error)
        area_m2 = area_model*metres_per_model_unit*metres_per_model_unit

        raw_obstacle_records = tree_records(obstacles)
        if not raw_obstacle_records:
            return empty_outputs(
                "INVALID_GEOMETRY",
                "No obstacle items are connected; case validation needs "
                "obstacle_count >= 1.")
        obstacle_records = []
        obstacle_errors = []
        obstacle_outside_centres = []
        for path, branch_index, item_index, value in raw_obstacle_records:
            geometry = resolve_rhino_value(value)
            label = "%s[%d]" % (path, item_index)
            if geometry is None:
                obstacle_errors.append(label+" did not resolve to geometry.")
                continue
            if not supported_obstacle(geometry):
                obstacle_errors.append(
                    "%s is %s; expected Curve, Mesh, Brep, Surface, "
                    "Extrusion, or SubD."
                    % (label, type(geometry).__name__))
                continue
            try:
                bounding_box = geometry.GetBoundingBox(False)
                valid_box = (
                    bounding_box.IsValid and
                    bounding_box.Diagonal.Length > tolerance)
            except Exception:
                valid_box = False
            if not valid_box:
                obstacle_errors.append(
                    label+" has an invalid or degenerate bounding box.")
                continue
            projected_center = plane.ClosestPoint(bounding_box.Center)
            try:
                containment = domain_curve.Contains(
                    projected_center, plane, tolerance)
                if containment == Rhino.Geometry.PointContainment.Outside:
                    obstacle_outside_centres.append(label)
            except Exception:
                pass
            obstacle_records.append(
                (path, branch_index, item_index, geometry))
        if obstacle_errors:
            return (
                0.0, False, [], [], "{}", "INVALID_GEOMETRY",
                ["INVALID OBSTACLE | "+message
                 for message in obstacle_errors])

        raw_flow_records = tree_records(flowVectors)
        if not raw_flow_records:
            return empty_outputs(
                "INVALID_FLOW",
                "No flowVectors are connected; each scenario must reference "
                "a vector whose planar length will later represent m/s.")
        flow_records = []
        flow_errors = []
        vertical_flow_count = 0
        flow_speeds = []
        for path, branch_index, item_index, value in raw_flow_records:
            vector = resolve_vector(value)
            label = "%s[%d]" % (path, item_index)
            if vector is None or not is_finite_vector(vector):
                flow_errors.append(label+" is not a finite Vector3d.")
                continue
            x = float(vector*plane.XAxis)
            y = float(vector*plane.YAxis)
            normal = float(vector*plane.Normal)
            planar_speed = math.sqrt(x*x+y*y)
            if planar_speed <= 1e-12:
                flow_errors.append(
                    label+" has zero planar speed or is perpendicular to domain.")
                continue
            if abs(normal) > max(1e-9, planar_speed*1e-6):
                vertical_flow_count += 1
            flow_records.append(
                (path, branch_index, item_index, (x, y)))
            flow_speeds.append(planar_speed)
        if flow_errors:
            return (
                0.0, False, [], [], "{}", "INVALID_FLOW",
                ["INVALID FLOW | "+message for message in flow_errors])

        probe_records = []
        rejected_probes = []
        outside_probes = []
        duplicate_probes = []
        tolerance_m = tolerance*metres_per_model_unit
        accepted_probe_xy_m = []
        for path, branch_index, item_index, value in tree_records(probes):
            point = resolve_point(value)
            label = "%s[%d]" % (path, item_index)
            if point is None or not is_finite_point(point):
                rejected_probes.append(label)
                continue
            projected = plane.ClosestPoint(point)
            delta = projected-plane.Origin
            xy_m = (
                float(delta*plane.XAxis)*metres_per_model_unit,
                float(delta*plane.YAxis)*metres_per_model_unit)
            if any(
                    (xy_m[0]-other[0])**2+(xy_m[1]-other[1])**2
                    <= tolerance_m*tolerance_m
                    for other in accepted_probe_xy_m):
                duplicate_probes.append(label)
                continue
            try:
                containment = domain_curve.Contains(
                    projected, plane, tolerance)
            except Exception:
                containment = Rhino.Geometry.PointContainment.Unset
            if containment == Rhino.Geometry.PointContainment.Outside:
                outside_probes.append(label)
                continue
            accepted_probe_xy_m.append(xy_m)
            probe_records.append(
                (path, branch_index, item_index, point))

        if caseJson is None or not str(caseJson).strip():
            return empty_outputs(
                "INVALID_CASE",
                "SimulationCaseJson is empty. Connect Site Data.SimulationCaseJson, or the base model for a manual case.")

        setup_ms = (time.perf_counter()-total_start)*1000.0
        case_start = time.perf_counter()
        try:
            case_core, optimizer_core, geometry_bridge = (
                load_runtime_modules(getattr(self, "Component", None)))
        except Exception as exception:
            return empty_outputs("INVALID_CASE", str(exception))

        try:
            compiled_case, case_cache_hit = compile_case_cached(
                case_core,
                caseJson,
                len(obstacle_records),
                len(flow_records))
            (case, digest, config, flow_indices,
             warnings, bridge_warnings) = compiled_case
        except case_core.UnsupportedCaseError as exception:
            return empty_outputs(
                "UNSUPPORTED_CASE", "UNSUPPORTED CASE | %s" % exception)
        except case_core.CaseValidationError as exception:
            report = [
                "INVALID CASE | optimizer received %d flattened flow vectors."
                % len(flow_records)]
            report.extend(
                "INVALID CASE | "+error for error in exception.errors)
            return (
                0.0, False, [], [], "{}", "INVALID_CASE", report)
        except Exception as exception:
            return empty_outputs(
                "INVALID_CASE", "CASE PARSER/COMPILER ERROR | %s" % exception)
        case_ms = (time.perf_counter()-case_start)*1000.0

        forcing = case["forcing"]
        scenario_key = "steps" if forcing["mode"] == "timeline" else "states"
        scenarios = forcing[scenario_key]

        geometry_start = time.perf_counter()
        try:
            domain_world = geometry_bridge.curve_samples(
                domain_curve, 64 if speed_mode else 200)
            domain_polygon = np.asarray([
                geometry_bridge.plane_coordinates(
                    point, plane, metres_per_model_unit)[:2]
                for point in domain_world], dtype=float)
            if (len(domain_polygon) < 3 or
                    optimizer_core.polygon_area(domain_polygon) <= 1e-12):
                return empty_outputs(
                    "INVALID_GEOMETRY",
                    "domain sampling did not produce a positive polygon.")

            geometries = [record[3] for record in obstacle_records]
            descriptors, centres, plan_clouds, inside_fractions, rejected_geo = (
                geometry_bridge.build_descriptors(
                    geometries, plane, metres_per_model_unit,
                    domain_polygon, tolerance_m, config,
                    quality_mode=quality_mode,
                    speed_mode=speed_mode))
            if rejected_geo or len(descriptors) != len(geometries):
                return empty_outputs(
                    "INVALID_GEOMETRY",
                    "Descriptor generation rejected obstacle indices: %s"
                    % ", ".join(str(index) for index in rejected_geo))
            boundary_value = (
                min(inside_fractions) if inside_fractions else 0.0)
            collision_value = geometry_bridge.collision_score(
                descriptors,
                config["constraint.min_obstacle_clearance_m"])
            probe_array = (
                np.asarray(accepted_probe_xy_m, dtype=float).reshape((-1, 2))
                if accepted_probe_xy_m
                else np.empty((0, 2), dtype=float))
            all_flow_xy = [record[3] for record in flow_records]
            selected_flow_xy = np.asarray(
                [all_flow_xy[index] for index in flow_indices], dtype=float)
            hydraulic_profiles = (
                None
                if speed_mode or quality_mode == 0 else
                geometry_bridge.build_hydraulic_profiles(
                    geometries, plane, metres_per_model_unit,
                    descriptors, selected_flow_xy))
            refined_mesh_count = (
                0 if hydraulic_profiles is None or not hydraulic_profiles
                else int(np.count_nonzero(
                    hydraulic_profiles[0]["mesh_supported"])))
        except Exception as exception:
            return empty_outputs(
                "INVALID_GEOMETRY",
                "DESCRIPTOR ERROR | %s" % exception)
        geometry_ms = (time.perf_counter()-geometry_start)*1000.0

        solver_start = time.perf_counter()
        try:
            result = optimizer_core.evaluate_layout(
                descriptors,
                domain_polygon,
                probe_array,
                selected_flow_xy,
                config,
                boundary_score=boundary_value,
                collision_score=collision_value,
                hydraulic_profiles=hydraulic_profiles)
        except Exception as exception:
            return empty_outputs("SOLVER_ERROR", "SOLVER ERROR | %s" % exception)
        solver_ms = (time.perf_counter()-solver_start)*1000.0
        pre_result_ms = (time.perf_counter()-total_start)*1000.0

        combined_warnings = list(warnings)+list(bridge_warnings)
        result_start = time.perf_counter()
        try:
            result_document = make_result_document(
                case, digest, result, scenarios, flow_indices,
                area_m2, obstacle_records, flow_records,
                pre_result_ms, geometry_ms, solver_ms, combined_warnings,
                fidelity, descriptor_name)
            canonical_result = case_core.canonical_json(
                json_safe(result_document))
        except Exception as exception:
            return empty_outputs(
                "SOLVER_ERROR", "RESULT SERIALIZATION ERROR | %s" % exception)
        result_ms = (time.perf_counter()-result_start)*1000.0
        total_ms = (time.perf_counter()-total_start)*1000.0

        status = result_document["run"]["status"]
        document = Rhino.RhinoDoc.ActiveDoc
        unit_name = (
            str(document.ModelUnitSystem)
            if document is not None else "Unset")
        domain_span_x = float(np.ptp(domain_polygon[:, 0]))
        domain_span_y = float(np.ptp(domain_polygon[:, 1]))
        flow_summary = " | ".join(
            "%s %.6f m/s" % (
                scenario["id"],
                result["scenario_speeds_m_s"][scenario_index])
            for scenario_index, scenario in enumerate(scenarios))
        objective_summary = " | ".join(
            "%s %s" % (
                name, display_number(result["objectives"][name]))
            for name in sorted(result["objectives"]))
        constraint_summary = " | ".join(
            "%s %s" % (
                name, display_number(result["constraint_margins"][name]))
            for name in sorted(result["constraint_margins"]))
        report = [
            "MUSSELFLOW | build %s | %d obstacles | %.3f m2 domain | "
            "%d probes (%s) | %d timeline scenarios | %.3f ms"
            % (COMPONENT_BUILD, len(obstacle_records), area_m2,
               len(result["probe_points"]), result["probe_source"],
               len(scenarios), total_ms),
            "TIMING | setup %.3f ms | case %.3f ms (%s) | geometry %.3f ms "
            "| solver %.3f ms | result %.3f ms"
            % (setup_ms, case_ms,
               "cache" if case_cache_hit else "compiled",
               geometry_ms, solver_ms, result_ms),
            "FIDELITY | %s | %s"
            % (fidelity, descriptor_name),
            "UNITS | Rhino %s | 1 model unit = %.6f m"
            % (unit_name, metres_per_model_unit),
            "DOMAIN | local bounds %.6f x %.6f m | area %.6f m2"
            % (domain_span_x, domain_span_y, area_m2),
            "FLOW | "+flow_summary,
            "FITNESS %s | raw %s | feasible %s | %s"
            % (display_number(result["fitness"]),
               display_number(result["raw_objective"]),
               result["feasible"], status),
            "OBJECTIVES | "+objective_summary,
            "CONSTRAINTS | "+constraint_summary,
            "CASE | %s | %s | hash %s"
            % (case["case_id"], forcing["mode"], digest),
            "RUNTIME | parser schema %s | %s"
            % (case_core.SCHEMA_VERSION,
               getattr(case_core, "__file__", "embedded/unknown")),
            "SCENARIOS | %d | source vector indices %s"
            % (len(scenarios), ", ".join(
                str(index) for index in flow_indices)),
            "SCREENING METRICS | cleared %s m3/day | chlorophyll %s "
            "g/day | particulate %s kg/day | minimum DO %s mg/L"
            % (display_number(result["effective_cleared_m3_day"]),
               display_number(result["chlorophyll_capture_g_day"]),
               display_number(result["particulate_capture_kg_day"]),
               display_number(result["oxygen"]["minimum_mg_l"])),
            "FEEDING PARTITION | ingested %s | assimilated %s | faeces %s "
            "| pseudofaeces %s kg organic/day"
            % (display_number(result["ingested_organic_kg_day"]),
               display_number(result["assimilated_organic_kg_day"]),
               display_number(result["faeces_organic_kg_day"]),
               display_number(result["pseudofaeces_organic_kg_day"])),
            "ENERGY PROXY | scope for growth %s kJ/day | potential growth "
            "%s g dry tissue/day (uncalibrated)"
            % (display_number(result["scope_for_growth_kj_day"]),
               display_number(result["potential_growth_g_dw_day"])),
            "SCIENTIFIC LIMIT | reduced-order descriptor/physics model; "
            "not CFD, site validation, carrying capacity, or legal evidence.",
            "HARVEST LIMIT | N/P harvest values in Result are prescribed "
            "standing-stock accounting, not predicted growth or verified "
            "eutrophication removal.",
            "ML STATUS | No trained surrogate or reinforcement-learning "
            "controller is active.",
        ]
        if unit_warning:
            report.append("UNIT WARNING | "+unit_warning)
        if obstacle_outside_centres:
            report.append(
                "GEOMETRY WARNING | %d obstacle bounding-box centre(s) are "
                "outside the domain: %s. Full boundary containment will be "
                "checked by the descriptor stage."
                % (len(obstacle_outside_centres),
                   ", ".join(obstacle_outside_centres)))
        if rejected_probes:
            report.append(
                "PROBE WARNING | Rejected non-point inputs: %s"
                % ", ".join(rejected_probes))
        if outside_probes:
            report.append(
                "PROBE WARNING | Ignored out-of-domain probes: %s"
                % ", ".join(outside_probes))
        if duplicate_probes:
            report.append(
                "PROBE WARNING | Removed duplicate probes: %s"
                % ", ".join(duplicate_probes))
        if vertical_flow_count:
            report.append(
                "FLOW WARNING | %d vector(s) have a normal component; only "
                "their domain-plane projection will be evaluated."
                % vertical_flow_count)
        if (not speed_mode and quality_mode == 1 and
                refined_mesh_count < len(geometries)):
            report.append(
                "FIDELITY WARNING | REFINE obtained projected mesh area for "
                "%d/%d obstacle(s); unsupported curves or geometry use their "
                "FAST projected envelope."
                % (refined_mesh_count, len(geometries)))
        if speed_mode and quality_mode == 1:
            report.append(
                "FIDELITY NOTE | speedMode overrides REFINE geometry for this "
                "evaluation; the complete ecological solver remains active.")
        report.extend(
            "CASE/BRIDGE WARNING | "+warning
            for warning in combined_warnings)
        report.extend(
            "MODEL WARNING | "+warning for warning in result["warnings"])
        return (
            compact_number(result["fitness"]),
            bool(result["feasible"]),
            named_values(result["objectives"]),
            named_values(result["constraint_margins"]),
            canonical_result,
            status,
            report,
        )
