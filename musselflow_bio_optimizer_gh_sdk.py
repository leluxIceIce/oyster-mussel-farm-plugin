# r: numpy
"""
MusselFlow Bio Optimizer - Rhino 8 Grasshopper Python 3 SDK mode.

This is the fast, non-animated Galapagos component.  Geometry is converted to
metres, several current vectors are evaluated as weighted scenarios, and a
strict ecological grammar supplies all biological/environmental coefficients.

Place ``musselflow_ecogrammar_core.py`` and
``musselflow_bio_optimizer_core.py`` beside the saved Grasshopper definition,
or in this script's development folder.

Inputs (the RunScript signature creates the type hints and access modes):
    run: Evaluate this design {item,bool}
    obstacles: Completed sock/star geometry; all branches are one farm
               {tree,GeometryBase}
    domain: Closed planar farm boundary {item,Curve}
    probes: Optional measurement points; obstacle centres are the explicit
            fallback when this tree is empty {tree,Point3d}
    flowVectors: Uniform-current scenarios.  Each vector's LENGTH is speed in
                 m/s and its direction is projected into the domain plane
                 {tree,Vector3d}
    grammar: Panel lines using ``key = value`` {tree,string}

Outputs - create/rename in this exact order:
    Fitness, RawObjective, Feasible, ObjectiveNames, ObjectiveValues,
    ConstraintNames, ConstraintMargins, ProbePoints, ProbeSpeed,
    ProbeSpeedRatio, ProbeFood, ProbeDeflection, ScenarioProbeSpeed,
    ScenarioProbeFood, ObstacleCenters, ObstacleSpeed, ObstacleFood,
    ObstacleRemoval, CaptureByObstacle, CaptureChl, CapturePM, ClearedWater,
    Biodeposit, MinDO, FinalDO, DOSeries, RespirationO2, ExcretionN,
    HarvestWet, HarvestN, HarvestP, ModelStatus, Report

Units:
    speed m/s; chlorophyll capture g/day; particulate and biodeposit kg/day;
    cleared water m3/day; DO mg/L; oxygen kg O2/day; excretion kg N/day;
    harvest wet t/year and N/P kg/year.

The node is an unvalidated reduced-order screening model until its calibration
flags are supported by site/flume/CFD evidence.
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

rc = Rhino


# =============================================================== conversion
def to_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    try:
        return list(value)
    except TypeError:
        return [value]


def flatten_tree(value):
    """Flatten a typed Grasshopper DataTree, with list fallback for old nodes."""
    if value is None:
        return []
    try:
        branches = value.Branches
    except Exception:
        return to_list(value)
    flattened = []
    for branch in branches:
        flattened.extend(list(branch))
    return flattened


def resolve_geometry(value):
    """Unwrap GH goo and legacy Guid inputs without requiring rhinoscriptsyntax."""
    if value is None:
        return None
    if hasattr(value, "Value"):
        try:
            value = value.Value
        except Exception:
            pass
    if isinstance(value, System.Guid):
        document = rc.RhinoDoc.ActiveDoc
        obj = document.Objects.FindId(value) if document is not None else None
        return obj.Geometry if obj is not None else None
    return value


def load_numeric_modules(component):
    """Load the testable numerical modules without hard-wiring them into Rhino."""
    module_names = (
        "musselflow_ecogrammar_core",
        "musselflow_bio_optimizer_core",
    )
    try:
        return tuple(importlib.import_module(name) for name in module_names)
    except ImportError:
        pass

    candidate_folders = []
    try:
        gh_document = component.OnPingDocument()
        definition_path = gh_document.FilePath if gh_document is not None else ""
        if definition_path:
            candidate_folders.append(os.path.dirname(definition_path))
    except Exception:
        pass
    candidate_folders.extend([
        os.getcwd(),
        "/Users/lelux/Desktop/Datajunk_Graphx_Entre/Python",
    ])
    for folder in candidate_folders:
        if not folder or folder in sys.path:
            continue
        if all(os.path.isfile(os.path.join(folder, name+".py"))
               for name in module_names):
            sys.path.insert(0, folder)
            importlib.invalidate_caches()
            try:
                return tuple(
                    importlib.import_module(name) for name in module_names)
            except ImportError:
                continue
    raise ImportError(
        "Could not load MusselFlow cores. Put musselflow_ecogrammar_core.py "
        "and musselflow_bio_optimizer_core.py beside the saved .gh file.")


def curve_plane(curve, tolerance):
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


def curve_samples(curve, count=160):
    try:
        success, polyline = curve.TryGetPolyline()
        if success and polyline.Count >= 3:
            limit = polyline.Count-1 if polyline.IsClosed else polyline.Count
            return [rc.Geometry.Point3d(polyline[index])
                    for index in range(limit)]
    except Exception:
        pass
    parameters = curve.DivideByCount(max(4, int(count)), False)
    if parameters:
        return [curve.PointAt(parameter) for parameter in parameters]
    return [curve.PointAtStart, curve.PointAtEnd]


def plane_coordinates(point, plane, metres_per_model_unit):
    delta = point-plane.Origin
    return (
        (delta*plane.XAxis)*metres_per_model_unit,
        (delta*plane.YAxis)*metres_per_model_unit,
        (delta*plane.Normal)*metres_per_model_unit,
    )


def world_point(x_m, y_m, z_m, plane, metres_per_model_unit):
    inverse = 1.0/metres_per_model_unit
    return (
        plane.Origin +
        plane.XAxis*(x_m*inverse) +
        plane.YAxis*(y_m*inverse) +
        plane.Normal*(z_m*inverse))


# ================================================================ geometry
def _sample_surface(surface, count_u=5, count_v=5):
    points = []
    try:
        domain_u = surface.Domain(0)
        domain_v = surface.Domain(1)
        for iu in range(count_u):
            u = domain_u.ParameterAt(
                float(iu)/float(max(count_u-1, 1)))
            for iv in range(count_v):
                v = domain_v.ParameterAt(
                    float(iv)/float(max(count_v-1, 1)))
                points.append(surface.PointAt(u, v))
    except Exception:
        pass
    return points


def _sample_brep(brep, maximum_vertices=160, maximum_edges=16):
    """Sample Brep boundaries without evaluating every face interior.

    The optimizer reduces every obstacle to plan dimensions, orientation and
    height. Vertices plus representative edge samples preserve those extents;
    repeatedly calling ``IsPointOnFace`` was the dominant Rhino runtime for
    detailed sock Breps while adding no information used by the descriptor.
    """
    vertex_count = brep.Vertices.Count
    if vertex_count >= 4:
        if vertex_count <= maximum_vertices:
            return [
                brep.Vertices[index].Location
                for index in range(vertex_count)]
        vertex_indices = sorted(set(
            int(index) for index in np.linspace(
                0, vertex_count-1, maximum_vertices)))
        return [
            brep.Vertices[index].Location
            for index in vertex_indices]

    # Smooth primitives such as cylinders may expose too few vertices to
    # preserve their plan extents. Only those exceptional Breps need edge
    # evaluation.
    points = [
        brep.Vertices[index].Location
        for index in range(vertex_count)]
    edge_count = brep.Edges.Count
    if edge_count <= maximum_edges:
        edge_indices = range(edge_count)
    else:
        edge_indices = sorted(set(
            int(index) for index in np.linspace(
                0, edge_count-1, maximum_edges)))
    for edge_index in edge_indices:
        points.extend(curve_samples(brep.Edges[edge_index], 4))
    if len(points) < 4:
        box = brep.GetBoundingBox(False)
        if box.IsValid:
            points.extend(list(box.GetCorners()))
    return points


def geometry_points(
        value, maximum=64, quality_mode=0, speed_mode=False):
    """Return a deterministic point cloud for fast or refined descriptors."""
    geometry = resolve_geometry(value)
    G = rc.Geometry
    if geometry is None:
        return []
    if isinstance(geometry, G.Point3d):
        return [geometry]
    if isinstance(geometry, G.Point):
        return [geometry.Location]
    if isinstance(geometry, G.Curve):
        sample_count = 8 if speed_mode else (24 if quality_mode else 12)
        return curve_samples(geometry, sample_count)
    if isinstance(geometry, G.Mesh):
        count = geometry.Vertices.Count
        if count <= maximum:
            return [G.Point3d(geometry.Vertices[index])
                    for index in range(count)]
        indices = np.linspace(0, count-1, maximum, dtype=int)
        points = [G.Point3d(geometry.Vertices[int(index)])
                  for index in indices]
        box = geometry.GetBoundingBox(False)
        if box.IsValid:
            points.extend(list(box.GetCorners()))
        return points
    if isinstance(geometry, G.Extrusion):
        geometry = geometry.ToBrep()
    elif isinstance(geometry, G.Surface):
        count = 3 if speed_mode else (7 if quality_mode else 4)
        return _sample_surface(geometry, count, count)
    elif hasattr(G, "SubD") and isinstance(geometry, G.SubD):
        try:
            geometry = geometry.ToBrep()
        except Exception:
            pass
    if isinstance(geometry, G.Brep):
        points = _sample_brep(
            geometry,
            maximum_vertices=(
                16 if speed_mode else 320 if quality_mode else maximum),
            maximum_edges=(
                6 if speed_mode else 48 if quality_mode else 12))
        if len(points) > maximum:
            indices = np.linspace(0, len(points)-1, maximum, dtype=int)
            points = [points[int(index)] for index in indices]
        return points
    try:
        box = geometry.GetBoundingBox(False)
        return list(box.GetCorners()) if box.IsValid else []
    except Exception:
        return []


def convex_hull(points):
    unique = sorted(set((float(point[0]), float(point[1]))
                        for point in points))
    if len(unique) <= 2:
        return unique

    def cross(origin, first, second):
        return (
            (first[0]-origin[0])*(second[1]-origin[1]) -
            (first[1]-origin[1])*(second[0]-origin[0]))

    lower = []
    for point in unique:
        while (len(lower) >= 2 and
               cross(lower[-2], lower[-1], point) <= 0.0):
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(unique):
        while (len(upper) >= 2 and
               cross(upper[-2], upper[-1], point) <= 0.0):
            upper.pop()
        upper.append(point)
    return lower[:-1]+upper[:-1]


def point_segment_distance_squared(point, first, second):
    px, py = point
    ax, ay = first
    bx, by = second
    dx, dy = bx-ax, by-ay
    length_squared = dx*dx+dy*dy
    if length_squared <= 1e-24:
        return (px-ax)**2+(py-ay)**2
    parameter = max(
        0.0, min(1.0, ((px-ax)*dx+(py-ay)*dy)/length_squared))
    qx = ax+parameter*dx
    qy = ay+parameter*dy
    return (px-qx)**2+(py-qy)**2


def point_in_polygon(point, polygon, tolerance):
    count = len(polygon)
    if count < 3:
        return False
    tolerance_squared = tolerance*tolerance
    for index in range(count):
        if point_segment_distance_squared(
                point, polygon[index], polygon[(index+1) % count]
                ) <= tolerance_squared:
            return True
    x, y = point
    inside = False
    previous = count-1
    for current in range(count):
        xi, yi = polygon[current]
        xj, yj = polygon[previous]
        if (yi > y) != (yj > y):
            crossing_x = (xj-xi)*(y-yi)/(yj-yi)+xi
            if x < crossing_x:
                inside = not inside
        previous = current
    return inside


def descriptor_from_points(points, plane, metres_per_model_unit,
                           fallback_plan_m, fallback_height_m):
    coordinates = np.asarray([
        plane_coordinates(point, plane, metres_per_model_unit)
        for point in points], dtype=float)
    if coordinates.ndim != 2 or len(coordinates) == 0:
        return None
    plan = coordinates[:, :2]
    centre_mean = np.mean(plan, axis=0)
    centred = plan-centre_mean
    if len(plan) >= 2 and np.any(np.abs(centred) > 1e-15):
        covariance = centred.T@centred/max(len(plan)-1, 1)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        major_axis = eigenvectors[:, int(np.argmax(eigenvalues))]
    else:
        major_axis = np.array([1.0, 0.0])
    if major_axis[0] < 0.0 or (
            abs(major_axis[0]) < 1e-15 and major_axis[1] < 0.0):
        major_axis *= -1.0
    minor_axis = np.array([-major_axis[1], major_axis[0]])
    major_values = plan@major_axis
    minor_values = plan@minor_axis
    major_min, major_max = np.min(major_values), np.max(major_values)
    minor_min, minor_max = np.min(minor_values), np.max(minor_values)
    major_length = max(float(major_max-major_min), fallback_plan_m)
    minor_length = max(float(minor_max-minor_min), fallback_plan_m)
    centre_plan = (
        major_axis*(0.5*(major_min+major_max)) +
        minor_axis*(0.5*(minor_min+minor_max)))
    if minor_length > major_length:
        major_length, minor_length = minor_length, major_length
        major_axis = minor_axis
    yaw = math.atan2(major_axis[1], major_axis[0])
    z_min = float(np.min(coordinates[:, 2]))
    z_max = float(np.max(coordinates[:, 2]))
    height = max(z_max-z_min, fallback_height_m)
    centre_z = 0.5*(z_min+z_max)
    descriptor = np.array([
        centre_plan[0], centre_plan[1], major_length, minor_length,
        yaw, height], dtype=float)
    return descriptor, centre_z, plan


def build_descriptors(geometries, plane, metres_per_model_unit,
                      domain_polygon, tolerance_m, config, quality_mode=0,
                      speed_mode=False):
    descriptors = []
    centres = []
    plan_clouds = []
    inside_fractions = []
    rejected = []
    for source_index, geometry in enumerate(geometries):
        points = geometry_points(
            geometry,
            maximum=(
                16 if speed_mode else 320 if quality_mode else 64),
            quality_mode=quality_mode,
            speed_mode=speed_mode)
        result = descriptor_from_points(
            points, plane, metres_per_model_unit,
            config["structure.fallback_plan_size_m"],
            config["structure.fallback_height_m"])
        if result is None:
            rejected.append(source_index)
            continue
        descriptor, centre_z, plan = result
        flags = [
            point_in_polygon(point, domain_polygon, tolerance_m)
            for point in plan]
        inside_fractions.append(
            sum(1 for flag in flags if flag)/float(max(len(flags), 1)))
        descriptors.append(descriptor)
        centres.append(world_point(
            descriptor[0], descriptor[1], centre_z, plane,
            metres_per_model_unit))
        plan_clouds.append(plan)
    return (
        np.asarray(descriptors, dtype=float),
        centres,
        plan_clouds,
        inside_fractions,
        rejected,
    )


def _hydraulic_mesh(value):
    """Return one analysis mesh or ``None`` without modifying source geometry."""
    geometry = resolve_geometry(value)
    G = rc.Geometry
    if geometry is None:
        return None
    if isinstance(geometry, G.Mesh):
        return geometry
    if isinstance(geometry, G.Extrusion):
        geometry = geometry.ToBrep()
    elif isinstance(geometry, G.Surface):
        geometry = geometry.ToBrep()
    elif hasattr(G, "SubD") and isinstance(geometry, G.SubD):
        try:
            geometry = geometry.ToBrep()
        except Exception:
            return None
    if not isinstance(geometry, G.Brep):
        return None
    try:
        meshes = G.Mesh.CreateFromBrep(
            geometry, G.MeshingParameters.FastRenderMesh)
    except Exception:
        return None
    if not meshes:
        return None
    joined = G.Mesh()
    for mesh in meshes:
        if mesh is not None:
            joined.Append(mesh)
    return joined if joined.Vertices.Count and joined.Faces.Count else None


def _triangle_projected_area(first, second, third, flow):
    area_vector = 0.5*np.cross(second-first, third-first)
    return abs(float(np.dot(area_vector, flow)))


def _projected_mesh_area(
        mesh, plane, metres_per_model_unit, flow_xy, maximum_faces=6000):
    """Approximate flow-normal projected area in square metres.

    Closed meshes expose front and back surfaces, so their summed projection is
    halved. Open meshes are treated as one hydraulic sheet. Concave overlap is
    capped later by the obstacle's projected envelope.
    """
    face_count = mesh.Faces.Count
    if face_count <= 0:
        return None
    if face_count <= maximum_faces:
        face_indices = list(range(face_count))
        sample_scale = 1.0
    else:
        face_indices = sorted(set(
            int(index) for index in np.linspace(
                0, face_count-1, maximum_faces)))
        sample_scale = face_count/float(len(face_indices))

    vertices = mesh.Vertices
    flow_length = math.sqrt(
        float(flow_xy[0])*float(flow_xy[0]) +
        float(flow_xy[1])*float(flow_xy[1]))
    if flow_length <= 1e-15:
        return None
    flow = np.array([
        float(flow_xy[0])/flow_length,
        float(flow_xy[1])/flow_length,
        0.0])
    coordinate_cache = {}

    def coordinate(index):
        cached = coordinate_cache.get(index)
        if cached is not None:
            return cached
        point = rc.Geometry.Point3d(vertices[index])
        value = np.asarray(
            plane_coordinates(point, plane, metres_per_model_unit),
            dtype=float)
        coordinate_cache[index] = value
        return value

    total = 0.0
    for face_index in face_indices:
        face = mesh.Faces[face_index]
        a = coordinate(face.A)
        b = coordinate(face.B)
        c = coordinate(face.C)
        total += _triangle_projected_area(a, b, c, flow)
        if face.IsQuad:
            d = coordinate(face.D)
            total += _triangle_projected_area(a, c, d, flow)
    total *= sample_scale
    try:
        if mesh.IsClosed:
            total *= 0.5
    except Exception:
        pass
    return total if math.isfinite(total) and total > 0.0 else None


def build_hydraulic_profiles(
        geometries, plane, metres_per_model_unit, descriptors, flow_vectors):
    """Build direction-specific frontal areas for refined fidelity.

    The returned list has one profile per flow scenario. Each profile retains
    the outer wake width but scales blockage and intercepted flux by the
    actual projected mesh area. Explicit textile porosity remains a separate
    ecological input.
    """
    descriptors = np.asarray(descriptors, dtype=float)
    flow_vectors = np.asarray(flow_vectors, dtype=float).reshape((-1, 2))
    meshes = [_hydraulic_mesh(geometry) for geometry in geometries]
    profiles = []
    for flow_vector in flow_vectors:
        flow, speed = np.asarray(flow_vector, dtype=float), float(
            np.linalg.norm(flow_vector))
        if speed > 1e-15:
            flow = flow/speed
        cross = np.array([-flow[1], flow[0]]) if speed > 1e-15 else np.array(
            [0.0, 1.0])
        yaw = descriptors[:, 4]
        major_axis = np.column_stack((np.cos(yaw), np.sin(yaw)))
        minor_axis = np.column_stack((-np.sin(yaw), np.cos(yaw)))
        cross_width = (
            descriptors[:, 2]*np.abs(major_axis@cross) +
            descriptors[:, 3]*np.abs(minor_axis@cross))
        envelope_area = np.maximum(
            cross_width*descriptors[:, 5], 1e-12)
        frontal_area = envelope_area.copy()
        mesh_supported = np.zeros(len(descriptors), dtype=bool)
        for obstacle_index, mesh in enumerate(meshes):
            if mesh is None:
                continue
            projected = _projected_mesh_area(
                mesh, plane, metres_per_model_unit, flow_vector)
            if projected is None:
                continue
            frontal_area[obstacle_index] = min(
                max(projected, 0.01*envelope_area[obstacle_index]),
                envelope_area[obstacle_index])
            mesh_supported[obstacle_index] = True
        profiles.append({
            "frontal_area_m2": frontal_area,
            "frontal_fill": np.clip(
                frontal_area/envelope_area, 0.01, 1.0),
            "mesh_supported": mesh_supported,
        })
    return profiles


def rectangle_corners(descriptor, clearance):
    centre = descriptor[:2]
    major = 0.5*descriptor[2]+clearance
    minor = 0.5*descriptor[3]+clearance
    axis = np.array([math.cos(descriptor[4]), math.sin(descriptor[4])])
    cross = np.array([-axis[1], axis[0]])
    return np.asarray([
        centre-major*axis-minor*cross,
        centre+major*axis-minor*cross,
        centre+major*axis+minor*cross,
        centre-major*axis+minor*cross,
    ])


def oriented_overlap_severity(first, second):
    axes = []
    for rectangle in (first, second):
        for index in (0, 1):
            edge = rectangle[(index+1) % 4]-rectangle[index]
            length = np.linalg.norm(edge)
            if length > 1e-15:
                axes.append(np.array([-edge[1], edge[0]])/length)
    penetration = []
    for axis in axes:
        first_projection = first@axis
        second_projection = second@axis
        overlap = (
            min(np.max(first_projection), np.max(second_projection)) -
            max(np.min(first_projection), np.min(second_projection)))
        if overlap <= 0.0:
            return 0.0
        penetration.append(overlap)
    scale = min(
        np.linalg.norm(first[1]-first[0]),
        np.linalg.norm(first[2]-first[1]),
        np.linalg.norm(second[1]-second[0]),
        np.linalg.norm(second[2]-second[1]))
    return float(np.clip(min(penetration)/max(scale, 1e-12), 0.0, 1.0))


def collision_score(descriptors, clearance):
    if len(descriptors) < 2:
        return 1.0
    rectangles = [
        rectangle_corners(descriptor, 0.5*clearance)
        for descriptor in descriptors]
    maximum_severity = 0.0
    for first_index in range(len(rectangles)-1):
        for second_index in range(first_index+1, len(rectangles)):
            maximum_severity = max(
                maximum_severity,
                oriented_overlap_severity(
                    rectangles[first_index], rectangles[second_index]))
            if maximum_severity >= 1.0:
                return 0.0
    return 1.0-maximum_severity


def unique_probes(points, plane, metres_per_model_unit,
                  domain_polygon, tolerance_m):
    accepted_world = []
    accepted_xy = []
    rejected = 0
    duplicate = 0
    quantisation = max(tolerance_m, 1e-9)
    seen = set()
    for value in points:
        point = resolve_geometry(value)
        if isinstance(point, rc.Geometry.Point):
            point = point.Location
        if not isinstance(point, rc.Geometry.Point3d):
            continue
        x, y, unused_z = plane_coordinates(
            point, plane, metres_per_model_unit)
        if not point_in_polygon((x, y), domain_polygon, tolerance_m):
            rejected += 1
            continue
        key = (int(round(x/quantisation)), int(round(y/quantisation)))
        if key in seen:
            duplicate += 1
            continue
        seen.add(key)
        accepted_world.append(point)
        accepted_xy.append((x, y))
    return accepted_world, accepted_xy, rejected, duplicate


def scenario_tree(values):
    tree = Grasshopper.DataTree[object]()
    for scenario_index, row in enumerate(np.asarray(values)):
        path = Grasshopper.Kernel.Data.GH_Path(scenario_index)
        for value in np.ravel(row):
            tree.Add(float(value), path)
    return tree


def empty_outputs(message):
    return (
        0.0, 0.0, False, [], [], [], [], [], [], [], [], [],
        Grasshopper.DataTree[object](), Grasshopper.DataTree[object](),
        [], [], [], [], [], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, [],
        0.0, 0.0, 0.0, 0.0, 0.0, "WAITING", [message])


# =========================================================== SDK component
class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def RunScript(
            self,
            run: bool,
            obstacles: Grasshopper.DataTree[
                Rhino.Geometry.GeometryBase],
            domain: Rhino.Geometry.Curve,
            probes: Grasshopper.DataTree[Rhino.Geometry.Point3d],
            flowVectors: Grasshopper.DataTree[
                Rhino.Geometry.Vector3d],
            grammar: Grasshopper.DataTree[str]):

        if not run:
            return empty_outputs("Waiting: set run to True.")
        total_start = time.perf_counter()

        try:
            grammar_core, optimizer_core = load_numeric_modules(
                getattr(self, "Component", None))
        except Exception as exception:
            return empty_outputs(str(exception))

        domain_curve = resolve_geometry(domain)
        if (domain_curve is None or
                not isinstance(domain_curve, rc.Geometry.Curve) or
                not domain_curve.IsClosed):
            return empty_outputs("domain must be one closed planar Curve.")

        document = rc.RhinoDoc.ActiveDoc
        tolerance_model = (
            document.ModelAbsoluteTolerance if document is not None else 0.001)
        plane = curve_plane(domain_curve, tolerance_model)
        if plane is None:
            return empty_outputs("domain Curve must be planar.")

        metres_per_model_unit = 1.0
        unit_warning = None
        if document is not None:
            try:
                metres_per_model_unit = rc.RhinoMath.UnitScale(
                    document.ModelUnitSystem, rc.UnitSystem.Meters)
            except Exception:
                unit_warning = (
                    "Could not read Rhino model units; geometry was treated as metres.")
        if not np.isfinite(metres_per_model_unit) or metres_per_model_unit <= 0:
            metres_per_model_unit = 1.0
            unit_warning = (
                "Rhino model units are unset; geometry was treated as metres.")
        tolerance_m = max(
            tolerance_model*metres_per_model_unit, 1e-9)

        domain_world = curve_samples(domain_curve, 200)
        domain_polygon = np.asarray([
            plane_coordinates(point, plane, metres_per_model_unit)[:2]
            for point in domain_world], dtype=float)
        if (len(domain_polygon) < 3 or
                optimizer_core.polygon_area(domain_polygon) <= 1e-12):
            return empty_outputs("domain has no measurable planar area.")

        geometry_values = [
            resolve_geometry(value) for value in flatten_tree(obstacles)]
        geometry_values = [
            value for value in geometry_values if value is not None]
        if not geometry_values:
            return empty_outputs(
                "No obstacles resolved. Input must use Tree access with "
                "Rhino.Geometry.GeometryBase type hint.")

        flow_values = []
        vertical_flow_count = 0
        normal_only_flow_count = 0
        for vector in flatten_tree(flowVectors):
            if vector is None:
                continue
            try:
                x = float(vector*plane.XAxis)
                y = float(vector*plane.YAxis)
                normal = float(vector*plane.Normal)
            except Exception:
                continue
            vector_length = math.sqrt(x*x+y*y+normal*normal)
            planar_length = math.sqrt(x*x+y*y)
            if vector_length > 1e-12 and planar_length <= 1e-12:
                normal_only_flow_count += 1
                continue
            if abs(normal) > 1e-9:
                vertical_flow_count += 1
            flow_values.append((x, y))
        if normal_only_flow_count:
            return empty_outputs(
                "FLOW VECTOR ERROR | %d vector(s) are perpendicular to the "
                "domain plane, so plan-view speed is zero. Use a vector in "
                "the domain plane (Unit X or Unit Y for an XY domain) and "
                "make its length the current speed in m/s, for example "
                "Unit X * 0.15." % normal_only_flow_count)
        if not flow_values:
            return empty_outputs(
                "flowVectors needs at least one Vector3d; vector length is m/s.")

        grammar_start = time.perf_counter()
        config, grammar_warnings, grammar_errors, notes = (
            grammar_core.parse_grammar(
                flatten_tree(grammar),
                obstacle_count=len(geometry_values),
                scenario_count=len(flow_values)))
        grammar_ms = (time.perf_counter()-grammar_start)*1000.0
        if grammar_errors:
            if any("flow scenarios" in error for error in grammar_errors):
                return empty_outputs(
                    "FLOW/GRAMMAR COUNT ERROR | %d flow vector(s) connected. "
                    "In the grammar use scalar defaults "
                    "'scenario.weights = 1' and "
                    "'scenario.duration_h = 6', or supply exactly one comma "
                    "value per flow vector. DETAILS | %s"
                    % (len(flow_values), " | ".join(grammar_errors)))
            return empty_outputs(
                "GRAMMAR ERROR | "+" | ".join(grammar_errors))

        geometry_start = time.perf_counter()
        descriptors, centres, plan_clouds, inside_fractions, rejected_geo = (
            build_descriptors(
                geometry_values, plane, metres_per_model_unit,
                domain_polygon, tolerance_m, config))
        if len(descriptors) == 0:
            return empty_outputs(
                "No obstacle produced a usable 3D/plan descriptor.")
        if len(descriptors) != len(geometry_values):
            return empty_outputs(
                "Some obstacles were rejected (%s). Fix geometry rather than "
                "silently changing stocking-list correspondence."
                % ", ".join(str(index) for index in rejected_geo))

        boundary_value = min(inside_fractions) if inside_fractions else 0.0
        collision_value = collision_score(
            descriptors, config["constraint.min_obstacle_clearance_m"])

        supplied_probes = flatten_tree(probes)
        accepted_world, accepted_xy, rejected_probes, duplicate_probes = (
            unique_probes(
                supplied_probes, plane, metres_per_model_unit,
                domain_polygon, tolerance_m))
        if supplied_probes and not accepted_xy:
            return empty_outputs(
                "Probe input was supplied, but no unique probe lies in domain.")
        probe_array = (
            np.asarray(accepted_xy, dtype=float).reshape((-1, 2))
            if accepted_xy else np.empty((0, 2), dtype=float))
        geometry_ms = (time.perf_counter()-geometry_start)*1000.0

        solver_start = time.perf_counter()
        try:
            result = optimizer_core.evaluate_layout(
                descriptors,
                domain_polygon,
                probe_array,
                np.asarray(flow_values, dtype=float),
                config,
                boundary_score=boundary_value,
                collision_score=collision_value)
        except Exception as exception:
            return empty_outputs("SOLVER ERROR | %s" % exception)
        solver_ms = (time.perf_counter()-solver_start)*1000.0

        output_probe_points = (
            accepted_world if result["probe_source"] == "user_probes"
            else centres)
        objective_names = list(result["objectives"].keys())
        objective_values = [
            float(result["objectives"][name]) for name in objective_names]
        constraint_names = list(result["constraint_margins"].keys())
        constraint_margins = [
            float(result["constraint_margins"][name])
            for name in constraint_names]
        total_ms = (time.perf_counter()-total_start)*1000.0

        report = [
            "BIO OPTIMIZER | %d obstacles | %d probes (%s) | %d currents | "
            "%.3f ms total"
            % (len(descriptors), len(output_probe_points),
               result["probe_source"], len(flow_values), total_ms),
            "TIMING | grammar %.3f ms | geometry %.3f ms | solver %.3f ms"
            % (grammar_ms, geometry_ms, solver_ms),
            "Fitness %.4f | raw objective %.4f | feasible %s | %s"
            % (result["fitness"], result["raw_objective"],
               result["feasible"], result["model_status"]),
            "Constraint domination: feasible fitness = 0.25 + 0.75*objective; "
            "infeasible fitness <= 0.25.",
            "Water cleared %.3f m3/day | chlorophyll captured %.3f g/day | "
            "particulate captured %.3f kg/day"
            % (result["effective_cleared_m3_day"],
               result["chlorophyll_capture_g_day"],
               result["particulate_capture_kg_day"]),
            "Organic biodeposit %.3f kg/day | minimum DO %.3f mg/L | "
            "mussel respiration %.4f kg O2/day"
            % (result["biodeposit_organic_kg_day"],
               result["oxygen"]["minimum_mg_l"],
               result["mussel_respiration_kg_o2_day"]),
            "Management scenario (not predicted growth): %.3f wet t/year | "
            "%.3f kg N/year | %.3f kg P/year"
            % (result["harvested_wet_t_year"],
               result["harvest_n_kg_year"],
               result["harvest_p_kg_year"]),
            "Flow vectors are uniform scenarios in m/s, not a spatial vector "
            "field. Weights: %s | durations h: %s"
            % (", ".join("%.3f" % value
                         for value in result["scenario_weights"]),
               ", ".join("%.3f" % value
                         for value in result["scenario_durations_h"])),
            "Calibration fraction %.0f%%. Physical-unit outputs remain "
            "screening estimates until relevant calibration flags are evidenced."
            % (100.0*result["calibration_fraction"]),
        ]
        if unit_warning:
            report.append(unit_warning)
        if vertical_flow_count:
            report.append(
                "%d flow vector(s) had a normal component; only the planar "
                "projection was evaluated." % vertical_flow_count)
        if rejected_probes:
            report.append(
                "%d out-of-domain probe(s) were ignored." % rejected_probes)
        if duplicate_probes:
            report.append(
                "%d duplicate probe(s) were removed to prevent accidental "
                "objective reweighting." % duplicate_probes)
        report.extend("GRAMMAR | "+warning for warning in grammar_warnings)
        report.extend("MODEL | "+warning for warning in result["warnings"])
        if notes:
            report.append(
                "%d note/source grammar entries preserved but not executed."
                % len(notes))

        return (
            result["fitness"],
            result["raw_objective"],
            result["feasible"],
            objective_names,
            objective_values,
            constraint_names,
            constraint_margins,
            output_probe_points,
            [float(value) for value in result["probe_speed_m_s"]],
            [float(value) for value in result["probe_speed_ratio"]],
            [float(value) for value in result["probe_food_fraction"]],
            [float(value) for value in result["probe_deflection_deg"]],
            scenario_tree(result["scenario_probe_speed_m_s"]),
            scenario_tree(result["scenario_probe_food_fraction"]),
            centres,
            [float(value) for value in result["obstacle_speed_m_s"]],
            [float(value) for value in result["obstacle_food_fraction"]],
            [float(value) for value in result["obstacle_removal_fraction"]],
            [float(value) for value in
             result["chlorophyll_capture_g_day_by_obstacle"]],
            result["chlorophyll_capture_g_day"],
            result["particulate_capture_kg_day"],
            result["effective_cleared_m3_day"],
            result["biodeposit_organic_kg_day"],
            result["oxygen"]["minimum_mg_l"],
            result["oxygen"]["final_mg_l"],
            [float(value) for value in result["oxygen"]["series_mg_l"]],
            result["mussel_respiration_kg_o2_day"],
            result["ammonia_excretion_kg_n_day"],
            result["harvested_wet_t_year"],
            result["harvest_n_kg_year"],
            result["harvest_p_kg_year"],
            result["model_status"],
            report)
