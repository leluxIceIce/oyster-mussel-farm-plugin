"""
Prepare real mussel-net Meshes/Surfaces for the MusselFlow solvers.

This Rhino 8 Grasshopper Python SDK-mode component intersects every incoming
Mesh, Brep, Surface, Extrusion or closed Curve with one analysis plane. It returns
closed footprint Curves for the flow solvers and filled planar Breps so the exact
hydraulic obstacle seen by the 2D model is visible immediately in Rhino.

No Point inputs are used: this component is for topology-driven net geometry.

    Inputs:
        run: Activates the component {item,bool}
        sockGeometry: Meshes, Breps, Surfaces, Extrusions or closed Curves
                      {list,geometry}
        sectionPlane: Plane at the water depth to analyse; WorldXY if empty
                      {item,plane}
        flowDir: Water-flow direction projected into sectionPlane {item,vector}

    Outputs (in this exact return order):
        Footprints: Closed planar Curves used as obstacles by both flow solvers
        SectionSurfaces: Filled planar Breps for shaded Rhino preview
        Centers: Area centroids corresponding to the footprints
        Areas: Section areas in model units squared
        FrontalWidths: Projected widths perpendicular to flow
        StreamLengths: Projected lengths parallel to flow
        SourceIndex: Index of the source sock geometry for each footprint
        Report: Per-object validation and summary

This is a 2D sectioning stage. It does not yet simulate water around the complete
3D sock surface. Use multiple planes later for a 2.5D stack, or CFD for full 3D.
"""

import System
import Rhino
import Grasshopper
import math

rc = Rhino


def unitizedInPlane(vector, plane):
    """Project a vector into a plane and return a valid unit flow direction."""
    if vector is None or not vector.IsValid or vector.Length < 1e-12:
        vector = plane.XAxis
    normal = plane.Normal
    dot = (vector.X*normal.X + vector.Y*normal.Y + vector.Z*normal.Z)
    projected = rc.Geometry.Vector3d(
        vector.X-dot*normal.X,
        vector.Y-dot*normal.Y,
        vector.Z-dot*normal.Z)
    if projected.Length < 1e-12:
        projected = plane.XAxis
    projected.Unitize()
    return projected


def closedSections(geometry, plane, tolerance):
    """Return closed section curves from one supported Rhino geometry object."""
    G = rc.Geometry
    curves = []

    if isinstance(geometry, G.Mesh):
        polylines = G.Intersect.Intersection.MeshPlane(geometry, plane)
        if polylines:
            curves.extend(G.PolylineCurve(polyline) for polyline in polylines)

    elif isinstance(geometry, G.Brep):
        success, sectionCurves, points = G.Intersect.Intersection.BrepPlane(
            geometry, plane, tolerance)
        if success and sectionCurves:
            curves.extend(sectionCurves)

    elif isinstance(geometry, G.Surface):
        success, sectionCurves, points = G.Intersect.Intersection.BrepPlane(
            geometry.ToBrep(), plane, tolerance)
        if success and sectionCurves:
            curves.extend(sectionCurves)

    elif isinstance(geometry, G.Extrusion):
        success, sectionCurves, points = G.Intersect.Intersection.BrepPlane(
            geometry.ToBrep(), plane, tolerance)
        if success and sectionCurves:
            curves.extend(sectionCurves)

    elif isinstance(geometry, G.SubD):
        brep = geometry.ToBrep()
        if brep:
            success, sectionCurves, points = G.Intersect.Intersection.BrepPlane(
                brep, plane, tolerance)
            if success and sectionCurves:
                curves.extend(sectionCurves)

    elif isinstance(geometry, G.Curve):
        duplicate = geometry.DuplicateCurve()
        if duplicate and duplicate.IsClosed:
            curves.append(duplicate)

    if not curves:
        return []

    joined = G.Curve.JoinCurves(curves, tolerance)
    candidates = list(joined) if joined else curves
    output = []
    for curve in candidates:
        if curve and curve.IsValid and curve.IsClosed:
            output.append(curve)
    return output


def footprintMetrics(curve, plane, flowDirection):
    """Calculate centroid, area and projected hydraulic dimensions."""
    properties = rc.Geometry.AreaMassProperties.Compute(curve)
    if properties is None:
        return None

    center = properties.Centroid
    area = float(properties.Area)
    parameters = curve.DivideByCount(64, True)
    points = ([curve.PointAt(parameter) for parameter in parameters]
              if parameters else [curve.PointAtStart])

    crossDirection = rc.Geometry.Vector3d.CrossProduct(
        plane.Normal, flowDirection)
    if crossDirection.Length < 1e-12:
        crossDirection = plane.YAxis
    crossDirection.Unitize()

    along = []
    across = []
    for point in points:
        delta = point-center
        along.append(delta.X*flowDirection.X +
                     delta.Y*flowDirection.Y +
                     delta.Z*flowDirection.Z)
        across.append(delta.X*crossDirection.X +
                      delta.Y*crossDirection.Y +
                      delta.Z*crossDirection.Z)

    streamLength = max(along)-min(along) if along else 0.0
    frontalWidth = max(across)-min(across) if across else 0.0
    return center, area, float(frontalWidth), float(streamLength)


def planarSurface(curve, tolerance):
    """Create one filled planar Brep for immediate shaded GH/Rhino preview."""
    breps = rc.Geometry.Brep.CreatePlanarBreps(curve, tolerance)
    if breps and len(breps):
        return breps[0]
    return None


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def RunScript(
            self,
            run: bool,
            sockGeometry: System.Collections.Generic.List[
                Rhino.Geometry.GeometryBase],
            sectionPlane: Rhino.Geometry.Plane,
            flowDir: Rhino.Geometry.Vector3d):

        empty = ([], [], [], [], [], [], [], ["Waiting for geometry."])
        if not run or not sockGeometry:
            return empty

        plane = (sectionPlane
                 if sectionPlane is not None and sectionPlane.IsValid
                 else rc.Geometry.Plane.WorldXY)
        direction = unitizedInPlane(flowDir, plane)
        document = rc.RhinoDoc.ActiveDoc
        tolerance = document.ModelAbsoluteTolerance if document else 0.001

        footprints = []
        sectionSurfaces = []
        centers = []
        areas = []
        frontalWidths = []
        streamLengths = []
        sourceIndices = []
        report = []

        for sourceIndex, geometry in enumerate(sockGeometry):
            if geometry is None:
                report.append("Geometry %d: empty input." % sourceIndex)
                continue

            sections = closedSections(geometry, plane, tolerance)
            if not sections:
                report.append(
                    "Geometry %d: no closed intersection at section plane."
                    % sourceIndex)
                continue

            accepted = 0
            for curve in sections:
                metrics = footprintMetrics(curve, plane, direction)
                if metrics is None:
                    report.append(
                        "Geometry %d: closed curve has no measurable area."
                        % sourceIndex)
                    continue

                center, area, frontalWidth, streamLength = metrics
                surface = planarSurface(curve, tolerance)
                footprints.append(curve)
                sectionSurfaces.append(surface)
                centers.append(center)
                areas.append(area)
                frontalWidths.append(frontalWidth)
                streamLengths.append(streamLength)
                sourceIndices.append(sourceIndex)
                accepted += 1

            report.append(
                "Geometry %d: %d valid hydraulic footprint(s)."
                % (sourceIndex, accepted))

        report.insert(
            0,
            "SECTION | %d source geometries -> %d closed footprints | Z %.3f"
            % (len(sockGeometry), len(footprints), plane.Origin.Z))
        if not footprints:
            report.append(
                "Move sectionPlane so it cuts the sock meshes/surfaces.")

        return (
            footprints,
            sectionSurfaces,
            centers,
            areas,
            frontalWidths,
            streamLengths,
            sourceIndices,
            report)
