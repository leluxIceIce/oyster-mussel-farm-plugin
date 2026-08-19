"""
MusselFlow Viewport Annotations — Rhino 8 Grasshopper Python SDK component.

Draws up to four construction-style dimensions around a closed planar boundary,
plus a horizontal colour legend, site coordinates, and date range pinned to the
Rhino viewport. Nothing is baked into the Rhino document.

SDK setup
---------
Create a Rhino 8 Python 3 component, choose ``Convert To GH_ScriptInstance``,
replace its generated source with this complete file, and remove every output.
The annotated RunScript signature creates these inputs in this exact order:

    visible, domain, latitude, longitude, startTime, endTime,
    color, textSize, font, sideCount, gradientColors

Keep both the component preview and ``visible`` enabled to see the annotation.
"""

import math

import Grasshopper
import Rhino
import System.Drawing


COMPONENT_METADATA = {
    "name": "MusselFlow Viewport Annotations",
    "nickname": "Site Annotations",
    "description": (
        "Viewport-only boundary dimensions and georeferenced site HUD. "
        "No Rhino objects are baked."
    ),
}

INPUT_METADATA = (
    ("visible", "visible", "True draws the dimensions and viewport HUD."),
    ("domain", "domain", "Closed planar rectangular boundary curve to dimension."),
    ("latitude", "latitude", "WGS84 latitude in decimal degrees."),
    ("longitude", "longitude", "WGS84 longitude in decimal degrees."),
    ("startTime", "startTime", "First displayed UTC ISO timestamp."),
    ("endTime", "endTime", "Last displayed UTC ISO timestamp."),
    ("color", "color", "Colour swatch for all lines, arrows, and text."),
    ("textSize", "textSize", "Viewport text height in pixels; default 28."),
    ("font", "font", "Optional installed font face; default Arial."),
    ("sideCount", "sideCount", "Number of boundary dimensions to draw, clamped to 0-4."),
    ("gradientColors", "gradient", "Ordered colour list for the horizontal low-to-high legend."),
)


def apply_component_metadata(component):
    if component is None:
        return
    component.Name = COMPONENT_METADATA["name"]
    component.NickName = COMPONENT_METADATA["nickname"]
    component.Description = COMPONENT_METADATA["description"]
    component.Message = "Viewport only"
    for index, (name, nickname, description) in enumerate(INPUT_METADATA):
        if index >= component.Params.Input.Count:
            break
        parameter = component.Params.Input[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description
        if index == 10:
            parameter.Access = Grasshopper.Kernel.GH_ParamAccess.list


def warning(component, message):
    if component is not None:
        component.AddRuntimeMessage(
            Grasshopper.Kernel.GH_RuntimeMessageLevel.Warning, message)


def tolerance():
    document = Rhino.RhinoDoc.ActiveDoc
    if document is None:
        return 0.001
    return max(float(document.ModelAbsoluteTolerance), 1e-9)


def metres_per_model_unit():
    document = Rhino.RhinoDoc.ActiveDoc
    if document is None:
        return None
    try:
        scale = Rhino.RhinoMath.UnitScale(
            document.ModelUnitSystem, Rhino.UnitSystem.Meters)
        return float(scale) if math.isfinite(scale) and scale > 0.0 else None
    except Exception:
        return None


def finite_number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def clean_points(points, tol):
    cleaned = []
    for point in points:
        point = Rhino.Geometry.Point3d(point)
        if not cleaned or point.DistanceTo(cleaned[-1]) > tol:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0].DistanceTo(cleaned[-1]) <= tol:
        cleaned.pop()

    changed = True
    while changed and len(cleaned) > 4:
        changed = False
        for index in range(len(cleaned)):
            previous = cleaned[index-1]
            current = cleaned[index]
            following = cleaned[(index+1) % len(cleaned)]
            before = current-previous
            after = following-current
            if before.Length <= tol or after.Length <= tol:
                cleaned.pop(index)
                changed = True
                break
            before.Unitize()
            after.Unitize()
            if before*after > 0.999999:
                cleaned.pop(index)
                changed = True
                break
    return cleaned


def curve_plane(curve):
    try:
        success, plane = curve.TryGetPlane()
        if success:
            return plane
    except Exception:
        pass
    raise ValueError("domain must be planar.")


def boundary_corners(curve, tol):
    """Return four ordered corners, using the curve plane as a safe fallback."""
    try:
        success, polyline = curve.TryGetPolyline()
        if success:
            points = clean_points(list(polyline), tol)
            if len(points) == 4:
                return points
    except Exception:
        pass

    try:
        segments = list(curve.DuplicateSegments())
        points = clean_points([segment.PointAtStart for segment in segments], tol)
        if len(points) == 4:
            return points
    except Exception:
        pass

    plane = curve_plane(curve)
    parameters = curve.DivideByCount(64, True)
    if not parameters:
        raise ValueError("The domain curve could not be sampled.")
    uv = []
    for parameter in parameters:
        point = curve.PointAt(parameter)
        delta = point-plane.Origin
        uv.append((delta*plane.XAxis, delta*plane.YAxis))
    u0 = min(item[0] for item in uv)
    u1 = max(item[0] for item in uv)
    v0 = min(item[1] for item in uv)
    v1 = max(item[1] for item in uv)
    if min(u1-u0, v1-v0) <= tol:
        raise ValueError("The domain has a degenerate planar bounding box.")
    return [
        plane.PointAt(u0, v0), plane.PointAt(u1, v0),
        plane.PointAt(u1, v1), plane.PointAt(u0, v1),
    ]


def average_point(points):
    count = float(len(points))
    return Rhino.Geometry.Point3d(
        sum(point.X for point in points)/count,
        sum(point.Y for point in points)/count,
        sum(point.Z for point in points)/count,
    )


def length_label(model_length, metres_per_unit):
    if metres_per_unit is None:
        return "%.3g units" % model_length
    metres = model_length*metres_per_unit
    if abs(metres) >= 1000.0:
        kilometres = metres/1000.0
        return ("%.0f km" % kilometres if abs(kilometres-round(kilometres)) < 1e-6
                else "%.2f km" % kilometres)
    if abs(metres) >= 10.0:
        return "%.0f m" % metres
    if abs(metres) >= 1.0:
        return "%.2f m" % metres
    return "%.0f mm" % (metres*1000.0)


def dimension_geometry(corners, plane, tol, side_count):
    centre = average_point(corners)
    side_lengths = [
        corners[(index+1) % 4].DistanceTo(corners[index])
        for index in range(4)
    ]
    shortest = min(length for length in side_lengths if length > tol)
    margin = max(8.0*tol, 0.065*shortest)
    extension_overhang = 0.25*margin
    arrow_length = 0.32*margin
    arrow_width = 0.16*margin
    unit_scale = metres_per_model_unit()
    lines = []
    labels = []
    clip_points = list(corners)

    for index in range(side_count):
        start = corners[index]
        end = corners[(index+1) % 4]
        tangent = end-start
        length = tangent.Length
        if length <= tol:
            continue
        tangent.Unitize()
        outward = Rhino.Geometry.Vector3d.CrossProduct(plane.Normal, tangent)
        if not outward.Unitize():
            continue
        midpoint = Rhino.Geometry.Point3d(
            0.5*(start.X+end.X), 0.5*(start.Y+end.Y), 0.5*(start.Z+end.Z))
        if outward*(centre-midpoint) > 0.0:
            outward.Reverse()

        dim_start = start+outward*margin
        dim_end = end+outward*margin
        ext_start_0 = start+outward*(0.18*margin)
        ext_start_1 = end+outward*(0.18*margin)
        ext_end_0 = dim_start+outward*extension_overhang
        ext_end_1 = dim_end+outward*extension_overhang
        lines.extend([
            Rhino.Geometry.Line(dim_start, dim_end),
            Rhino.Geometry.Line(ext_start_0, ext_end_0),
            Rhino.Geometry.Line(ext_start_1, ext_end_1),
        ])

        for tip, inward in ((dim_start, tangent), (dim_end, -tangent)):
            lines.append(Rhino.Geometry.Line(
                tip, tip+inward*arrow_length+outward*arrow_width))
            lines.append(Rhino.Geometry.Line(
                tip, tip+inward*arrow_length-outward*arrow_width))

        text_point = midpoint+outward*(margin+0.58*margin)
        labels.append({
            "text": length_label(length, unit_scale),
            "point": text_point,
            "tangent": tangent,
            "outward": outward,
            "line_start": dim_start,
            "line_end": dim_end,
        })
        clip_points.extend([ext_end_0, ext_end_1, text_point])

    box = Rhino.Geometry.BoundingBox(clip_points)
    box.Inflate(max(tol, 0.25*margin))
    return lines, labels, box


def normalized_color(value):
    if isinstance(value, System.Drawing.Color) and not value.IsEmpty:
        return value
    return System.Drawing.Color.FromArgb(255, 59, 48)


def normalized_text_size(value):
    try:
        return max(8, min(160, int(round(float(value)))))
    except (TypeError, ValueError):
        return 28


def normalized_font(value):
    text = "" if value is None else str(value).strip()
    return text or "Arial"


def normalized_side_count(value):
    if value is None:
        return 4
    try:
        return max(0, min(4, int(round(float(value)))))
    except (TypeError, ValueError):
        return 4


def normalized_gradient(value):
    if value is None:
        return []
    try:
        values = list(value)
    except TypeError:
        values = [value]
    return [item for item in values
            if isinstance(item, System.Drawing.Color) and not item.IsEmpty]


def gradient_color(colors, parameter):
    if len(colors) == 1:
        return colors[0]
    position = max(0.0, min(1.0, parameter))*(len(colors)-1)
    index = min(len(colors)-2, int(math.floor(position)))
    local = position-index
    first = colors[index]
    second = colors[index+1]
    channel = lambda a, b: int(round(a+(b-a)*local))
    return System.Drawing.Color.FromArgb(
        channel(first.A, second.A),
        channel(first.R, second.R),
        channel(first.G, second.G),
        channel(first.B, second.B))


def draw_hud(display, viewport, draw_state):
    size = viewport.Size
    x = 0.785*float(size.Width)
    rows = (
        (draw_state["start_time"], 0.060),
        (draw_state["end_time"], 0.102),
        (draw_state["longitude"], 0.800),
        (draw_state["latitude"], 0.840),
    )
    for text, y_ratio in rows:
        if text:
            display.Draw2dText(
                text, draw_state["color"],
                Rhino.Geometry.Point2d(x, y_ratio*float(size.Height)),
                False, draw_state["text_size"], draw_state["font"])


def draw_gradient_legend(display, viewport, draw_state):
    colors = draw_state["gradient_colors"]
    if not colors:
        return
    size = viewport.Size
    x0 = 0.785*float(size.Width)
    x1 = 0.950*float(size.Width)
    y0 = 0.905*float(size.Height)
    bar_height = max(8.0, 0.55*draw_state["text_size"])
    segments = max(32, min(160, int(round(x1-x0))))
    for index in range(segments):
        parameter = index/float(max(1, segments-1))
        x = x0+(x1-x0)*parameter
        display.Draw2dLine(
            System.Drawing.PointF(x, y0),
            System.Drawing.PointF(x, y0+bar_height),
            gradient_color(colors, parameter), 2.0)

    border = draw_state["color"]
    display.Draw2dLine(
        System.Drawing.PointF(x0, y0), System.Drawing.PointF(x1, y0), border, 1.0)
    display.Draw2dLine(
        System.Drawing.PointF(x1, y0),
        System.Drawing.PointF(x1, y0+bar_height), border, 1.0)
    display.Draw2dLine(
        System.Drawing.PointF(x1, y0+bar_height),
        System.Drawing.PointF(x0, y0+bar_height), border, 1.0)
    display.Draw2dLine(
        System.Drawing.PointF(x0, y0+bar_height),
        System.Drawing.PointF(x0, y0), border, 1.0)

    label_y = y0+bar_height+1.15*draw_state["text_size"]
    display.Draw2dText(
        "low", border, Rhino.Geometry.Point2d(x0, label_y), False,
        draw_state["text_size"], draw_state["font"])
    display.Draw2dText(
        "high", border, Rhino.Geometry.Point2d(x1, label_y), True,
        draw_state["text_size"], draw_state["font"])


def draw_dimension_label(display, viewport, label, draw_state):
    point = label["point"]
    tangent = Rhino.Geometry.Vector3d(label["tangent"])
    outward = Rhino.Geometry.Vector3d(label["outward"])
    screen_start = viewport.WorldToClient(label["line_start"])
    screen_end = viewport.WorldToClient(label["line_end"])
    dx = screen_end.X-screen_start.X
    dy = screen_end.Y-screen_start.Y

    # Draw3dText is mirrored when its plane normal faces away from the camera.
    # Flip only the Y axis first so the text plane always faces the viewer.
    text_normal = Rhino.Geometry.Vector3d.CrossProduct(tangent, outward)
    to_camera = viewport.CameraLocation-point
    if text_normal*to_camera < 0.0:
        outward.Reverse()

    # Then rotate both axes together when needed, preserving the camera-facing
    # normal while keeping the baseline readable from left to right/upwards.
    if dx < -1e-6 or (abs(dx) <= 1e-6 and dy > 0.0):
        tangent.Reverse()
        outward.Reverse()
    text_plane = Rhino.Geometry.Plane(point, tangent, outward)
    try:
        success, pixels_per_unit = viewport.GetWorldToScreenScale(point)
        if not success or pixels_per_unit <= 1e-9:
            raise ValueError("No world-to-screen scale.")
        world_height = draw_state["text_size"]/float(pixels_per_unit)
        display.Draw3dText(
            label["text"], draw_state["color"], text_plane, world_height,
            draw_state["font"], False, False,
            Rhino.DocObjects.TextHorizontalAlignment.Center,
            Rhino.DocObjects.TextVerticalAlignment.Middle)
    except Exception:
        display.Draw2dText(
            label["text"], draw_state["color"],
            viewport.WorldToClient(point), True,
            draw_state["text_size"], draw_state["font"])


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def BeforeRunScript(self):
        apply_component_metadata(getattr(self, "Component", None))

    def RunScript(
            self,
            visible: bool,
            domain: Rhino.Geometry.Curve,
            latitude: float,
            longitude: float,
            startTime: str,
            endTime: str,
            color: System.Drawing.Color,
            textSize: int,
            font: str,
            sideCount: int,
            gradientColors: list[System.Drawing.Color]):
        """Prepare transient viewport dimensions and site metadata."""
        self._draw = None
        self._draw_error_reported = False
        if not bool(visible):
            return
        if not isinstance(domain, Rhino.Geometry.Curve):
            warning(self.Component, "Connect a closed planar boundary curve to domain.")
            return
        if not domain.IsClosed:
            warning(self.Component, "domain must be a closed planar curve.")
            return

        try:
            tol = tolerance()
            plane = curve_plane(domain)
            corners = boundary_corners(domain, tol)
            if len(corners) != 4:
                raise ValueError("domain did not resolve to four boundary sides.")
            side_count = normalized_side_count(sideCount)
            lines, labels, box = dimension_geometry(
                corners, plane, tol, side_count)
            lat = finite_number(latitude)
            lon = finite_number(longitude)
            self._draw = {
                "lines": lines,
                "labels": labels,
                "box": box,
                "color": normalized_color(color),
                "text_size": normalized_text_size(textSize),
                "font": normalized_font(font),
                "side_count": side_count,
                "gradient_colors": normalized_gradient(gradientColors),
                "start_time": "" if startTime is None else str(startTime),
                "end_time": "" if endTime is None else str(endTime),
                "latitude": "" if lat is None else "%.5f lat" % lat,
                "longitude": "" if lon is None else "%.5f lon" % lon,
            }
        except Exception as exception:
            warning(self.Component, "Viewport annotation: %s" % exception)

    @property
    def ClippingBox(self):
        draw_state = getattr(self, "_draw", None)
        if draw_state:
            return draw_state["box"]
        return Rhino.Geometry.BoundingBox.Empty

    def DrawViewportWires(self, args):
        draw_state = getattr(self, "_draw", None)
        if not draw_state:
            return
        try:
            display = args.Display
            viewport = getattr(args, "Viewport", None) or display.Viewport
            thickness = max(1, min(5, int(round(draw_state["text_size"]/14.0))))
            for line in draw_state["lines"]:
                display.DrawLine(line, draw_state["color"], thickness)
            for label in draw_state["labels"]:
                draw_dimension_label(display, viewport, label, draw_state)
            draw_hud(display, viewport, draw_state)
            draw_gradient_legend(display, viewport, draw_state)
        except Exception as exception:
            if not getattr(self, "_draw_error_reported", False):
                self._draw_error_reported = True
                warning(self.Component, "Viewport draw failed: %s" % exception)

    def AfterRunScript(self):
        pass
