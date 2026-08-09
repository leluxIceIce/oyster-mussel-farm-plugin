"""
Build a Kangaroo-ready membrane through ordered, transformed frame curves.

Important:
    Create all frames from ONE base closed curve using Move/Rotate transforms.
    Feed those transformed curves directly into this component, in stack order.
    Do not Rebuild them and do not change their seams upstream.

    Inputs:
        run: Build the mesh {item,bool}
        frames: Closed frame curves in physical stack order {list,curve}
        aroundCount: Requested vertices around each frame {item,int}
        betweenCount: Free rings inserted between frame pairs {item,int}

    Outputs (this exact order):
        Mesh: Connected quad mesh
        AnchorPoints: Mesh vertices located on the rigid frames
        TargetPoints: Matching target points for Kangaroo Anchor
        FramePolylines: Sampled frame rings
        Report: Validation and topology summary

Kangaroo:
    Mesh -> EdgeLength (LengthFactor 1.0)
    AnchorPoints -> Anchor P
    TargetPoints -> Anchor T
"""

import System
import Rhino
import Grasshopper

rc = Rhino


def polyline_vertices(curve):
    """Return unique vertices when curve is a closed polyline."""
    try:
        success, polyline = curve.TryGetPolyline()
        if not success or polyline.Count < 4:
            return None
        count = polyline.Count - 1 if polyline.IsClosed else polyline.Count
        return [rc.Geometry.Point3d(polyline[i]) for i in range(count)]
    except Exception:
        return None


def sample_polyline(vertices, requested):
    """Sample every edge equally while preserving seam and direction."""
    edge_count = len(vertices)
    per_edge = max(1, int(round(requested / float(edge_count))))
    count = edge_count * per_edge
    points = []
    for edge in range(edge_count):
        a = vertices[edge]
        b = vertices[(edge + 1) % edge_count]
        for step in range(per_edge):
            t = step / float(per_edge)
            points.append(rc.Geometry.Point3d(
                a.X + (b.X-a.X)*t,
                a.Y + (b.Y-a.Y)*t,
                a.Z + (b.Z-a.Z)*t))
    return points, count


def sample_curve(curve, count):
    """Fallback for smooth closed curves; preserves existing parameterization."""
    points = []
    domain = curve.Domain
    for index in range(count):
        t = domain.ParameterAt(index / float(count))
        points.append(curve.PointAt(t))
    return points


def interpolate_ring(first, second, parameter):
    return [
        rc.Geometry.Point3d(
            a.X + (b.X-a.X)*parameter,
            a.Y + (b.Y-a.Y)*parameter,
            a.Z + (b.Z-a.Z)*parameter)
        for a, b in zip(first, second)
    ]


def closed_polyline(points):
    return rc.Geometry.PolylineCurve(points + [points[0]])


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def RunScript(
            self,
            run: bool,
            frames: System.Collections.Generic.List[Rhino.Geometry.Curve],
            aroundCount: int,
            betweenCount: int):

        empty = (None, [], [], [], ["Waiting for frames."])
        if not run or not frames:
            return empty

        curves = [
            curve.DuplicateCurve() for curve in frames
            if curve is not None and curve.IsValid and curve.IsClosed
        ]
        if len(curves) < 2:
            return (None, [], [], [], [
                "Need at least two valid closed curves."])

        requested = max(4, int(aroundCount) if aroundCount else 32)
        between = max(
            0, int(betweenCount) if betweenCount is not None else 3)

        vertex_sets = [polyline_vertices(curve) for curve in curves]
        all_polylines = all(vertices is not None for vertices in vertex_sets)

        if all_polylines:
            edge_counts = [len(vertices) for vertices in vertex_sets]
            if len(set(edge_counts)) != 1:
                return (None, [], [], [], [
                    "All polyline frames must have the same number of edges.",
                    "Use transformed copies of one base frame."])
            rings = []
            actual = None
            for vertices in vertex_sets:
                ring, actual = sample_polyline(vertices, requested)
                rings.append(ring)
        else:
            actual = requested
            rings = [sample_curve(curve, actual) for curve in curves]

        # No sorting, seam shifting, reversing, or nearest-point matching:
        # ring[i][j] is the same material knot on every transformed frame.
        rows = [rings[0]]
        frame_rows = [0]
        for frame_index in range(len(rings)-1):
            first = rings[frame_index]
            second = rings[frame_index+1]
            for step in range(1, between+1):
                rows.append(interpolate_ring(
                    first, second, step / float(between+1)))
            rows.append(second)
            frame_rows.append(len(rows)-1)

        mesh = rc.Geometry.Mesh()
        for row in rows:
            for point in row:
                mesh.Vertices.Add(point)

        for row_index in range(len(rows)-1):
            a0 = row_index * actual
            b0 = (row_index+1) * actual
            for column in range(actual):
                following = (column+1) % actual
                mesh.Faces.AddFace(
                    int(a0+column),
                    int(a0+following),
                    int(b0+following),
                    int(b0+column))

        mesh.Normals.ComputeNormals()
        mesh.Compact()

        anchor_points = []
        for row_index in frame_rows:
            anchor_points.extend(rows[row_index])
        target_points = [
            rc.Geometry.Point3d(point) for point in anchor_points]
        frame_polylines = [closed_polyline(ring) for ring in rings]
        report = [
            "FRAME MEMBRANE | %d ordered frames | %d vertices/frame"
            % (len(rings), actual),
            "Mesh: %d vertices | %d quads | %d anchored vertices"
            % (mesh.Vertices.Count, mesh.Faces.Count, len(anchor_points)),
            "Correspondence is preserved directly from input curve "
            "parameterization.",
            "FramePolylines preserve each input curve's seam and direction.",
            "Bypass upstream Rebuild and Seam components."
        ]

        return (
            mesh,
            anchor_points,
            target_points,
            frame_polylines,
            report)
