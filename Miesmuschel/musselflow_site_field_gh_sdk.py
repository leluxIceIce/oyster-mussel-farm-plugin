"""
MusselFlow Copernicus Field
===========================

Self-contained Rhino 8 Grasshopper Python 3 SDK-mode component. It turns one
MusselFlow SiteData JSON document into a spatial Copernicus field sampled inside
a closed planar Rhino domain. Raw values remain in physical source units while
the preview uses normalized colour and circle size only for visual comparison.

This component uses Python's standard library, RhinoCommon, Grasshopper, and
System.Drawing included with Rhino. It needs no external executable, desktop
helper, Python package, account credential, or operating-system-specific path.

Important data distinction
--------------------------
The component samples processed Copernicus model or ocean-colour products. It
does not download raw Sentinel multispectral bands and does not reproduce SNAP
atmospheric correction, masking, custom indices, or a multispectral data cube.

Inputs
------
fetch : item / bool
    Connect a Button. True fetches missing samples. Identical URLs are cached in
    memory, so a held Toggle does not repeatedly contact Copernicus.
SiteDataJson : item / str
    SiteDataJson output from MusselFlow Site Data. It supplies the coordinate,
    regional layer catalogue, and available timestamps.
domain : item / Rhino.Geometry.Curve
    Closed planar Rhino curve. Its bounding-box centre is georeferenced to the
    requested SiteData coordinate. Rhino document units are converted to metres.
variable : item / str
    current_speed, temperature, salinity, oxygen, chlorophyll,
    satellite_chlorophyll, phytoplankton_carbon, nitrate, phosphate, tsm,
    or turbidity.
timeIndex : item / str
    Frame index or ISO UTC timestamp from SiteData.Times. Numeric indices are
    clamped; timestamps select the exact or nearest available frame.
depth : item / float
    Positive sampling depth in metres. Surface products ignore this value.
resolution : item / int
    Cells along the domain's longer side, clamped to 4-24. Twelve is recommended.
sizePower : item / float
    Hotspot radius exponent. 1 is linear; 2 strongly emphasizes high values.
colors : list / System.Drawing.Color
    Optional ordered colour stops from low to high. Leave empty for the default
    blue-cyan-yellow-red-purple sequence.
northVector : item / Rhino.Geometry.Vector3d
    Geographic north in the Rhino domain plane. Defaults to World +Y.
placementPoint : item / Rhino.Geometry.Point3d
    Optional Rhino point at which the centre of the visualized tile is placed.
    This translates preview geometry only; it never changes the SiteData
    latitude/longitude used for Copernicus sampling.

Outputs — create nine ports once in this exact order
----------------------------------------------------
FieldMesh : item / Rhino.Geometry.Mesh
    Disconnected coloured cells for every valid sample inside the domain.
HotspotMesh : item / Rhino.Geometry.Mesh
    Filled circles whose radius follows normalized_value ** sizePower.
Circles : list / Rhino.Geometry.Curve
    Circle outlines matching HotspotMesh, useful for downstream geometry.
Points : list / Rhino.Geometry.Point3d
    Rhino sample points corresponding one-to-one with Values and Normalized.
Values : list / float
    Unnormalized physical source values.
Normalized : list / float
    Display-only values from 0 to 1.
FieldDataJson : item / str
    Canonical ML-ready JSON with requested/source coordinates and raw values.
Legend : list / str
    Variable, source units, raw range, display transform, time, and depth.
Report : list / str
    Timing, sample counts, API/cache state, resolution warnings, and limits.

SDK setup
---------
Create a Rhino 8 Python 3 component, convert the default component with
``Convert To GH_ScriptInstance``, then replace its generated text with this
complete file. The RunScript annotations create and type-hint the eleven inputs.
Add nine output sockets once in the exact order above. BeforeRunScript applies
all names and human-readable hover tooltips.
"""

import concurrent.futures
import datetime
import hashlib
import json
import math
import time
import urllib.parse
import urllib.request

import Grasshopper
import Rhino
import System.Drawing


COMPONENT_BUILD = "2026-08-10a"
WMTS_ENDPOINT = "https://wmts.marine.copernicus.eu/teroWmts/"
TILE_LEVEL = 10
HTTP_TIMEOUT_S = 25
MAX_WORKERS = 8
MAX_SAMPLES = 576

VARIABLES = {
    "current_speed": {
        "layer": "current", "units": "m/s", "kind": "vector_speed",
        "daily": False, "surface": False,
    },
    "temperature": {
        "layer": "temperature", "units": "degC", "kind": "scalar",
        "daily": False, "surface": False,
    },
    "salinity": {
        "layer": "salinity", "units": "PSU", "kind": "scalar",
        "daily": False, "surface": False,
    },
    "oxygen": {
        "layer": "oxygen", "units": "mg/L", "kind": "oxygen",
        "daily": True, "surface": False,
    },
    "chlorophyll": {
        "layer": "chlorophyll", "units": "ug/L", "kind": "scalar",
        "daily": True, "surface": False,
    },
    "satellite_chlorophyll": {
        "layer": "satellite_chlorophyll", "units": "ug/L",
        "kind": "scalar", "daily": True, "surface": True,
    },
    "phytoplankton_carbon": {
        "layer": "phytoplankton_carbon", "units": "mmol/m3",
        "kind": "scalar", "daily": True, "surface": False,
    },
    "nitrate": {
        "layer": "nitrate", "units": "mmol/m3", "kind": "scalar",
        "daily": True, "surface": False,
    },
    "phosphate": {
        "layer": "phosphate", "units": "mmol/m3", "kind": "scalar",
        "daily": True, "surface": False,
    },
    "tsm": {
        "layer": "tsm", "units": "mg/L", "kind": "scalar",
        "daily": True, "surface": True,
    },
    "turbidity": {
        "layer": "turbidity", "units": "FNU", "kind": "scalar",
        "daily": True, "surface": True,
    },
}

ALIASES = {
    "speed": "current_speed",
    "current": "current_speed",
    "flow": "current_speed",
    "current_speed_m_s": "current_speed",
    "temp": "temperature",
    "thetao": "temperature",
    "temperature_c": "temperature",
    "so": "salinity",
    "salinity_psu": "salinity",
    "o2": "oxygen",
    "do": "oxygen",
    "dissolved_oxygen": "oxygen",
    "dissolved_oxygen_mg_l": "oxygen",
    "chl": "chlorophyll",
    "chlorophyll_a": "chlorophyll",
    "chlorophyll_a_ug_l": "chlorophyll",
    "satellite_chl": "satellite_chlorophyll",
    "ocean_colour_chlorophyll": "satellite_chlorophyll",
    "ocean_color_chlorophyll": "satellite_chlorophyll",
    "phyc": "phytoplankton_carbon",
    "no3": "nitrate",
    "po4": "phosphate",
    "spm": "tsm",
    "suspended_matter": "tsm",
    "tur": "turbidity",
}

COMPONENT_METADATA = {
    "name": "MusselFlow Copernicus Field",
    "nickname": "SiteField",
    "description": (
        "Sample one Copernicus variable over a georeferenced Rhino domain and "
        "output coloured hotspot geometry plus ML-ready physical values."),
}

INPUT_METADATA = (
    ("fetch", "fetch", "Button: fetch missing field samples; repeated URLs use memory cache."),
    ("SiteDataJson", "SiteData", "Actual sampled SiteDataJson from MusselFlow Site Data."),
    ("domain", "domain", "Closed planar curve defining the sampled tile size and shape."),
    ("variable", "variable", "Physical field name, e.g. chlorophyll, oxygen, nitrate, or current_speed."),
    ("timeIndex", "timeIndex", "Frame index or ISO UTC timestamp; timestamps select the exact or nearest SiteData frame."),
    ("depth", "depth", "Positive depth in metres; ignored by surface satellite products."),
    ("resolution", "resolution", "Grid cells along the longer domain side, 4-24; twelve recommended."),
    ("sizePower", "sizePower", "Hotspot radius exponent; 2 emphasizes high values."),
    ("colors", "colors", "Optional ordered System.Drawing colour stops, low to high."),
    ("northVector", "northVector", "Geographic north in the Rhino domain plane; defaults to World +Y."),
    ("placementPoint", "placementPoint", "Optional Rhino point for the displayed tile centre; geospatial sampling is unchanged."),
)

OUTPUT_METADATA = (
    ("FieldMesh", "FieldMesh", "Vertex-coloured grid cells for valid physical samples."),
    ("HotspotMesh", "HotspotMesh", "Filled value-scaled hotspot circles."),
    ("Circles", "Circles", "Hotspot circle outline curves."),
    ("Points", "Points", "Rhino sample points matching Values and Normalized."),
    ("Values", "Values", "Raw physical values in the Legend units."),
    ("Normalized", "Normalized", "Display-only values in the interval 0-1."),
    ("FieldDataJson", "FieldData", "Canonical spatial field records for export or later ML."),
    ("Legend", "Legend", "Variable, units, range, timestamp, and display mapping."),
    ("Report", "Report", "Request/cache state, coverage, warnings, and scientific scope."),
)

_HTTP_CACHE = {}


def apply_component_metadata(component):
    if component is None:
        return
    component.Name = COMPONENT_METADATA["name"]
    component.NickName = COMPONENT_METADATA["nickname"]
    component.Description = COMPONENT_METADATA["description"]
    component.Message = "Copernicus field"
    for index, (name, nickname, description) in enumerate(INPUT_METADATA):
        if index >= component.Params.Input.Count:
            break
        parameter = component.Params.Input[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description
        if index == 8:
            parameter.Access = Grasshopper.Kernel.GH_ParamAccess.list
    for index, (name, nickname, description) in enumerate(OUTPUT_METADATA):
        if index >= component.Params.Output.Count:
            break
        parameter = component.Params.Output[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description


def empty_outputs(status, message):
    return (None, None, [], [], [], [], "{}", [], [
        "MUSSELFLOW SITE FIELD | build "+COMPONENT_BUILD+" | "+status,
        message,
    ])


def finite_number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def rounded(value, digits=6):
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, float):
        return round(value, digits)
    if isinstance(value, list):
        return [rounded(item, digits) for item in value]
    if isinstance(value, dict):
        return {key: rounded(item, digits) for key, item in value.items()}
    return value


def canonical_json(value):
    return json.dumps(
        rounded(value), sort_keys=True, separators=(",", ":"),
        ensure_ascii=False, allow_nan=False)


def variable_name(value):
    text = "chlorophyll" if value is None else str(value).strip().lower()
    text = text.replace("-", "_").replace(" ", "_")
    return ALIASES.get(text, text)


def parse_site_data(source):
    try:
        document = json.loads(str(source))
    except Exception as exception:
        raise ValueError("siteData is invalid JSON: %s" % exception)
    if not isinstance(document, dict):
        raise ValueError("siteData root must be a JSON object.")
    if not isinstance(document.get("layers"), dict):
        raise ValueError("siteData is missing its regional layer catalogue.")
    frames = document.get("frames")
    if not isinstance(frames, list) or not frames:
        raise ValueError("siteData contains no valid frames.")
    requested = document.get("requested") or {}
    latitude = finite_number(requested.get("latitude"))
    longitude = finite_number(requested.get("longitude"))
    if latitude is None or longitude is None:
        raise ValueError("siteData is missing requested latitude/longitude.")
    return document, latitude, longitude, frames


def timestamp_text(value):
    """Return a stable UTC-like string for Python, .NET, or Grasshopper values."""
    if value is None:
        return ""
    try:
        return value.ToUniversalTime().ToString("yyyy-MM-ddTHH:mm:ss'Z'")
    except Exception:
        return str(value).strip()


def parse_utc_timestamp(value):
    text = timestamp_text(value)
    if not text:
        return None
    candidate = text[:-1]+"+00:00" if text.endswith("Z") else text
    try:
        parsed = datetime.datetime.fromisoformat(candidate)
    except (TypeError, ValueError):
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=datetime.timezone.utc)
    return parsed.astimezone(datetime.timezone.utc)


def select_frame(frames, selector):
    """Resolve a zero-based index or timestamp to one SiteData frame."""
    timestamps = [str(frame.get("time_utc") or "").strip() for frame in frames]
    if any(not value for value in timestamps):
        raise ValueError("one or more SiteData frames have no time_utc value.")

    text = timestamp_text(selector)
    if not text:
        return 0, timestamps[0], "default frame 0"

    try:
        numeric = float(text)
    except (TypeError, ValueError):
        numeric = None
    if numeric is not None and math.isfinite(numeric):
        requested = int(numeric)
        index = max(0, min(len(frames)-1, requested))
        note = "frame index %d" % requested
        if index != requested:
            note += " clamped to %d" % index
        return index, timestamps[index], note

    for index, timestamp in enumerate(timestamps):
        if text == timestamp:
            return index, timestamp, "exact timestamp"

    requested_time = parse_utc_timestamp(text)
    candidates = [
        (index, parse_utc_timestamp(timestamp))
        for index, timestamp in enumerate(timestamps)
    ]
    candidates = [
        (index, value) for index, value in candidates if value is not None
    ]
    if requested_time is None or not candidates:
        raise ValueError(
            "timeIndex must be a zero-based frame index or ISO timestamp "
            "present in SiteData.Times; received '%s'." % text)

    index, _ = min(
        candidates,
        key=lambda item: abs((item[1]-requested_time).total_seconds()))
    return index, timestamps[index], "nearest timestamp to %s" % text


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
            return [Rhino.Geometry.Point3d(polyline[index])
                    for index in range(limit)]
    except Exception:
        pass
    parameters = curve.DivideByCount(max(8, int(count)), False)
    if parameters:
        return [curve.PointAt(parameter) for parameter in parameters]
    return []


def model_scale_and_tolerance():
    document = Rhino.RhinoDoc.ActiveDoc
    if document is None:
        return 1.0, 0.001, "No active Rhino document; model units treated as metres."
    tolerance = max(float(document.ModelAbsoluteTolerance), 1e-12)
    try:
        scale = Rhino.RhinoMath.UnitScale(
            document.ModelUnitSystem, Rhino.UnitSystem.Meters)
    except Exception:
        return 1.0, tolerance, "Could not read Rhino units; model units treated as metres."
    if not math.isfinite(scale) or scale <= 0.0:
        return 1.0, tolerance, "Rhino units are unset; model units treated as metres."
    return float(scale), tolerance, None


def plane_uv(point, plane):
    delta = point-plane.Origin
    return delta*plane.XAxis, delta*plane.YAxis


def projected_north(value, plane):
    try:
        north = Rhino.Geometry.Vector3d(value)
    except Exception:
        north = Rhino.Geometry.Vector3d.YAxis
    north -= plane.Normal*(north*plane.Normal)
    if not north.IsValid or north.Length <= 1e-9:
        north = Rhino.Geometry.Vector3d(plane.YAxis)
    north.Unitize()
    east = Rhino.Geometry.Vector3d.CrossProduct(north, plane.Normal)
    east.Unitize()
    return east, north


def grid_samples(curve, plane, tolerance, resolution, metres_per_unit,
                 anchor_latitude, anchor_longitude, north_vector):
    boundary = curve_samples(curve)
    if len(boundary) < 3:
        raise ValueError("domain could not be sampled as a closed planar curve.")
    uv = [plane_uv(point, plane) for point in boundary]
    u_min = min(item[0] for item in uv)
    u_max = max(item[0] for item in uv)
    v_min = min(item[1] for item in uv)
    v_max = max(item[1] for item in uv)
    width = u_max-u_min
    height = v_max-v_min
    if width <= tolerance or height <= tolerance:
        raise ValueError("domain has a degenerate planar bounding box.")
    longest = max(width, height)
    nx = max(2, int(round(resolution*width/longest)))
    ny = max(2, int(round(resolution*height/longest)))
    while nx*ny > MAX_SAMPLES:
        if nx >= ny:
            nx -= 1
        else:
            ny -= 1
    du = width/float(nx)
    dv = height/float(ny)
    centre_u = 0.5*(u_min+u_max)
    centre_v = 0.5*(v_min+v_max)
    anchor = plane.PointAt(centre_u, centre_v)
    east_axis, north_axis = projected_north(north_vector, plane)
    longitude_scale = max(
        1e-9, 111320.0*math.cos(math.radians(anchor_latitude)))
    samples = []
    for row in range(ny):
        v = v_min+(row+0.5)*dv
        for column in range(nx):
            u = u_min+(column+0.5)*du
            point = plane.PointAt(u, v)
            containment = curve.Contains(point, plane, tolerance)
            if containment == Rhino.Geometry.PointContainment.Outside:
                continue
            delta = point-anchor
            east_m = (delta*east_axis)*metres_per_unit
            north_m = (delta*north_axis)*metres_per_unit
            samples.append({
                "row": row,
                "column": column,
                "u": u,
                "v": v,
                "point": point,
                "latitude": anchor_latitude+north_m/111320.0,
                "longitude": anchor_longitude+east_m/longitude_scale,
            })
    return (
        samples, nx, ny, du, dv,
        width*metres_per_unit, height*metres_per_unit, anchor)


def wmts_position(latitude, longitude, level=TILE_LEVEL):
    matrix_width = 2**(level+1)
    matrix_height = 2**level
    x = (longitude+180.0)/360.0*matrix_width
    y = (90.0-latitude)/180.0*matrix_height
    column = min(matrix_width-1, max(0, int(math.floor(x))))
    row = min(matrix_height-1, max(0, int(math.floor(y))))
    pixel_x = min(255, max(0, int(math.floor((x-column)*256.0))))
    pixel_y = min(255, max(0, int(math.floor((y-row)*256.0))))
    return row, column, pixel_x, pixel_y


def feature_url(layer, latitude, longitude, timestamp, depth=None):
    row, column, pixel_x, pixel_y = wmts_position(latitude, longitude)
    query = {
        "service": "WMTS",
        "request": "GetFeatureInfo",
        "layer": layer,
        "tilematrixset": "EPSG:4326",
        "tilematrix": str(TILE_LEVEL),
        "tilerow": str(row),
        "tilecol": str(column),
        "i": str(pixel_x),
        "j": str(pixel_y),
        "INFOFORMAT": "application/json",
        "time": timestamp,
    }
    if depth is not None:
        query["elevation"] = str(-abs(float(depth)))
    return WMTS_ENDPOINT+"?"+urllib.parse.urlencode(query)


def http_json(url):
    request = urllib.request.Request(
        url, headers={"User-Agent": "MusselFlow-Rhino/"+COMPONENT_BUILD})
    with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_S) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_urls(urls, allow_network):
    unique = sorted(set(urls))
    missing = [url for url in unique if url not in _HTTP_CACHE]
    errors = {}
    if missing and allow_network:
        workers = min(MAX_WORKERS, len(missing))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(http_json, url): url for url in missing}
            for future in concurrent.futures.as_completed(futures):
                url = futures[future]
                try:
                    _HTTP_CACHE[url] = future.result()
                except Exception as exception:
                    errors[url] = "%s: %s" % (
                        type(exception).__name__, str(exception))
    documents = {url: _HTTP_CACHE[url] for url in unique if url in _HTTP_CACHE}
    return documents, errors, len(missing), len(unique)-len(missing)


def feature_properties(document):
    try:
        features = document.get("features") or []
        return features[0].get("properties") or {}
    except (AttributeError, IndexError):
        return {}


def field_value(properties, kind):
    if kind == "vector_speed":
        east = finite_number(properties.get("component1Value"))
        north = finite_number(properties.get("component2Value"))
        return None if east is None or north is None else math.hypot(east, north)
    value = finite_number(properties.get("value"))
    if value is not None and kind == "oxygen":
        return value*0.032
    return value


def daily_time(timestamp):
    return str(timestamp)[:10]+"T00:00:00Z"


def default_colors():
    color = System.Drawing.Color
    return [
        color.FromArgb(25, 65, 180),
        color.FromArgb(0, 210, 230),
        color.FromArgb(250, 235, 40),
        color.FromArgb(245, 70, 30),
        color.FromArgb(105, 25, 145),
    ]


def color_stops(values):
    result = []
    if values is not None:
        try:
            candidates = list(values)
        except TypeError:
            candidates = [values]
        for value in candidates:
            if isinstance(value, System.Drawing.Color):
                result.append(value)
    return result if len(result) >= 2 else default_colors()


def interpolate_color(stops, value):
    t = min(1.0, max(0.0, float(value)))
    position = t*(len(stops)-1)
    index = min(len(stops)-2, int(math.floor(position)))
    local = position-index
    a, b = stops[index], stops[index+1]
    channel = lambda x, y: int(round(x+(y-x)*local))
    return System.Drawing.Color.FromArgb(
        channel(a.A, b.A), channel(a.R, b.R),
        channel(a.G, b.G), channel(a.B, b.B))


def add_colored_quad(mesh, plane, u, v, half_u, half_v, color):
    start = mesh.Vertices.Count
    points = (
        plane.PointAt(u-half_u, v-half_v),
        plane.PointAt(u+half_u, v-half_v),
        plane.PointAt(u+half_u, v+half_v),
        plane.PointAt(u-half_u, v+half_v),
    )
    for point in points:
        mesh.Vertices.Add(point)
        mesh.VertexColors.Add(color.R, color.G, color.B)
    mesh.Faces.AddFace(int(start), int(start+1), int(start+2), int(start+3))


def add_colored_disc(mesh, plane, centre, radius, color, segments=16):
    start = mesh.Vertices.Count
    mesh.Vertices.Add(centre)
    mesh.VertexColors.Add(color.R, color.G, color.B)
    local_plane = Rhino.Geometry.Plane(centre, plane.XAxis, plane.YAxis)
    for index in range(segments):
        angle = 2.0*math.pi*index/float(segments)
        point = local_plane.PointAt(radius*math.cos(angle), radius*math.sin(angle))
        mesh.Vertices.Add(point)
        mesh.VertexColors.Add(color.R, color.G, color.B)
    for index in range(segments):
        a = start+1+index
        b = start+1+((index+1) % segments)
        mesh.Faces.AddFace(int(start), int(a), int(b))


def make_geometry(valid_samples, plane, du, dv, size_power, stops):
    field_mesh = Rhino.Geometry.Mesh()
    hotspot_mesh = Rhino.Geometry.Mesh()
    circles = []
    points = []
    values = []
    normalized = []
    for sample in valid_samples:
        value = sample["value"]
        display = sample["normalized"]
        color = interpolate_color(stops, display)
        add_colored_quad(
            field_mesh, plane, sample["u"], sample["v"],
            0.48*du, 0.48*dv, color)
        radius = min(du, dv)*(0.06+0.40*(display**size_power))
        display_point = sample.get("display_point", sample["point"])
        centre = display_point+plane.Normal*(0.005*min(du, dv))
        add_colored_disc(hotspot_mesh, plane, centre, radius, color)
        circle_plane = Rhino.Geometry.Plane(centre, plane.XAxis, plane.YAxis)
        circles.append(Rhino.Geometry.Circle(circle_plane, radius).ToNurbsCurve())
        points.append(display_point)
        values.append(value)
        normalized.append(display)
    for mesh in (field_mesh, hotspot_mesh):
        mesh.Normals.ComputeNormals()
        mesh.Compact()
    return field_mesh, hotspot_mesh, circles, points, values, normalized


def normalize_samples(samples):
    values = [sample["value"] for sample in samples if sample["value"] is not None]
    if not values:
        return [], None, None
    low, high = min(values), max(values)
    span = high-low
    valid = []
    for sample in samples:
        value = sample["value"]
        if value is None:
            continue
        item = dict(sample)
        item["normalized"] = 0.5 if span <= 1e-15 else (value-low)/span
        valid.append(item)
    return valid, low, high


def source_cell_count(samples):
    cells = set()
    for sample in samples:
        lat = sample.get("source_latitude")
        lon = sample.get("source_longitude")
        if lat is not None and lon is not None:
            cells.add((round(lat, 6), round(lon, 6)))
    return len(cells)


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def BeforeRunScript(self):
        apply_component_metadata(getattr(self, "Component", None))

    def RunScript(
            self,
            fetch: bool,
            siteData: str,
            domain: Rhino.Geometry.Curve,
            variable: str,
            timeIndex: str,
            depth: float,
            resolution: int,
            sizePower: float,
            colors: list[System.Drawing.Color],
            northVector: Rhino.Geometry.Vector3d,
            placementPoint: Rhino.Geometry.Point3d):
        """Sample and preview one physical Copernicus field."""
        started = time.perf_counter()
        try:
            site, anchor_lat, anchor_lon, frames = parse_site_data(siteData)
            selected_variable = variable_name(variable)
            if selected_variable not in VARIABLES:
                raise ValueError(
                    "unknown variable '%s'. Choose: %s" % (
                        selected_variable, ", ".join(sorted(VARIABLES))))
            specification = VARIABLES[selected_variable]
            layer = site["layers"].get(specification["layer"])
            if not layer:
                raise ValueError(
                    "siteData regional catalogue has no '%s' layer."
                    % specification["layer"])
            index, frame_timestamp, time_selection = select_frame(
                frames, timeIndex)
            timestamp = (
                daily_time(frame_timestamp)
                if specification["daily"] else frame_timestamp)
            sample_depth = abs(finite_number(depth) or 0.0)
            count = 12 if resolution is None else int(resolution)
            count = max(4, min(24, count))
            power = finite_number(sizePower)
            power = 2.0 if power is None else min(4.0, max(0.25, power))
            if domain is None or not isinstance(domain, Rhino.Geometry.Curve):
                raise ValueError("domain must be one closed planar Rhino curve.")
            if not domain.IsClosed:
                raise ValueError("domain curve must be closed.")
            metres_per_unit, tolerance, unit_warning = model_scale_and_tolerance()
            plane = curve_plane(domain, tolerance)
            if plane is None:
                raise ValueError("domain curve must be planar.")
            samples, nx, ny, du, dv, width_m, height_m, domain_anchor = grid_samples(
                domain, plane, tolerance, count, metres_per_unit,
                anchor_lat, anchor_lon, northVector)
            try:
                display_anchor = Rhino.Geometry.Point3d(placementPoint)
                if not display_anchor.IsValid:
                    raise ValueError()
            except Exception:
                display_anchor = Rhino.Geometry.Point3d(domain_anchor)
            translation = display_anchor-domain_anchor
            display_plane = Rhino.Geometry.Plane(
                plane.Origin+translation, plane.XAxis, plane.YAxis)
            for sample in samples:
                sample["display_point"] = sample["point"]+translation
        except Exception as exception:
            return empty_outputs("INVALID_INPUT", str(exception))

        urls = []
        for sample in samples:
            query_depth = None if specification["surface"] else sample_depth
            url = feature_url(
                layer, sample["latitude"], sample["longitude"],
                timestamp, query_depth)
            sample["url"] = url
            urls.append(url)

        documents, errors, missing_before, cached_before = fetch_urls(urls, bool(fetch))
        if missing_before and not fetch:
            return empty_outputs(
                "WAITING",
                "Press the fetch Button once. %d unique API samples are not cached."
                % missing_before)

        for sample in samples:
            properties = feature_properties(documents.get(sample["url"], {}))
            sample["value"] = field_value(properties, specification["kind"])
            sample["source_latitude"] = finite_number(properties.get("lat"))
            sample["source_longitude"] = finite_number(properties.get("lon"))
            sample.pop("url", None)

        valid_samples, low, high = normalize_samples(samples)
        if not valid_samples:
            detail = ""
            if errors:
                detail = " First HTTP error: "+next(iter(errors.values()))
            return empty_outputs(
                "NO_DATA",
                "No valid %s values were returned for this domain/time/depth.%s"
                % (selected_variable, detail))

        stops = color_stops(colors)
        field_mesh, hotspot_mesh, circles, points, values, normalized = (
            make_geometry(valid_samples, display_plane, du, dv, power, stops))
        unique_cells = source_cell_count(valid_samples)
        region = site.get("regional_product") or {}
        field_document = {
            "schema": "musselflow.site_field.1.0.0",
            "build": COMPONENT_BUILD,
            "variable": selected_variable,
            "units": specification["units"],
            "time_utc": timestamp,
            "source_frame_index": index,
            "source_frame_time_utc": frame_timestamp,
            "depth_m": None if specification["surface"] else sample_depth,
            "regional_product": region,
            "anchor": {"latitude": anchor_lat, "longitude": anchor_lon},
            "placement_point_rhino": {
                "x": display_anchor.X,
                "y": display_anchor.Y,
                "z": display_anchor.Z,
            },
            "domain_size_m": {"width": width_m, "height": height_m},
            "grid": {
                "columns": nx,
                "rows": ny,
                "inside_samples": len(samples),
                "valid_samples": len(valid_samples),
                "unique_source_cells": unique_cells,
            },
            "display": {
                "minimum": low,
                "maximum": high,
                "normalization": "linear min-max",
                "hotspot_radius_power": power,
            },
            "records": [
                {
                    "sample_id": "%d:%d" % (
                        sample["row"], sample["column"]),
                    "row": sample["row"],
                    "column": sample["column"],
                    "latitude": sample["latitude"],
                    "longitude": sample["longitude"],
                    "source_latitude": sample["source_latitude"],
                    "source_longitude": sample["source_longitude"],
                    "value": sample["value"],
                    "normalized": sample["normalized"],
                    "rhino_point": {
                        "x": sample["display_point"].X,
                        "y": sample["display_point"].Y,
                        "z": sample["display_point"].Z,
                    },
                }
                for sample in valid_samples
            ],
            "scientific_status": "COPERNICUS_DERIVED_FIELD_NOT_SITE_VALIDATED",
        }
        legend = [
            "VARIABLE | %s" % selected_variable,
            "UNITS | %s" % specification["units"],
            "RAW RANGE | %.6g to %.6g" % (low, high),
            "DISPLAY | linear colour; hotspot radius = normalized^%.3g" % power,
            "TIME | %s" % timestamp,
            "DEPTH | %s" % (
                "surface product" if specification["surface"]
                else "%.6g m" % sample_depth),
        ]
        report = [
            "MUSSELFLOW SITE FIELD | build %s | %s | %d/%d valid | %.3fs"
            % (COMPONENT_BUILD, selected_variable, len(valid_samples),
               len(samples), time.perf_counter()-started),
            "REGION | %s" % (region.get("name") or "from SiteData catalogue"),
            "DOMAIN | %.3f x %.3f m | grid %d x %d | %d inside"
            % (width_m, height_m, nx, ny, len(samples)),
            "PLACEMENT | Rhino %.3f, %.3f, %.3f | SiteData WGS84 unchanged"
            % (display_anchor.X, display_anchor.Y, display_anchor.Z),
            "API | %d uncached before run | %d cached | %d HTTP failures"
            % (missing_before, cached_before, len(errors)),
            "SOURCE CELLS | %d distinct Copernicus cells represented"
            % unique_cells,
            "TIME SELECTION | %s | frame %d | %s"
            % (time_selection, index, frame_timestamp),
        ]
        if unit_warning:
            report.append("UNIT WARNING | "+unit_warning)
        if unique_cells <= 1 and len(valid_samples) > 1:
            report.append(
                "RESOLUTION WARNING | the Rhino domain is smaller than this "
                "product's spatial support; all samples resolve to one source cell.")
        if max(width_m, height_m) > 100000.0:
            report.append(
                "GEOREFERENCE WARNING | domain exceeds 100 km; the local tangent "
                "coordinate approximation should be replaced by a map projection.")
        if errors:
            report.append(
                "HTTP WARNING | failed samples remain holes; no values were invented.")
        report.extend([
            "DATA LIMIT | processed model/ocean-colour product, not raw multispectral bands.",
            "VALIDATION LIMIT | regional fields require product-QC and local observations.",
        ])
        return (
            field_mesh, hotspot_mesh, circles, points, values, normalized,
            canonical_json(field_document), legend, report)

    def AfterRunScript(self):
        pass
