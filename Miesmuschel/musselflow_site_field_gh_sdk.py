"""
MusselFlow Copernicus Field — Georeferenced Environmental Preview
=================================================================

Turns one Copernicus field into a georeferenced Rhino tile for rapid site
exploration. It preserves the WGS84 anchor and physical source values, then draws
coloured cells and value-scaled circles directly in the document.

The component converts the complete Rhino boundary into one WGS84 bounding box,
downloads the native numerical Copernicus grid once when the optional Copernicus
Toolbox is already available, caches that patch, and maps any Rhino preview
resolution onto it locally. It deliberately does not install the Toolbox into
Rhino's shared Python environment: its large dependency tree can invalidate every
Python component if package installation fails. If the Toolbox is unavailable,
the component reports the missing bulk backend and does not issue per-cell WMTS
requests.

``exportFieldData`` controls only whether the canonical FieldDataJson document is
serialized for the PLS/ML pipeline; it never changes sampling, source values, or
preview geometry.

Name: MusselFlow Copernicus Field
Updated: 260818
Author: Felix Berger
Copyright: Apache License 2.0

Inputs:
fetch : item / bool
    Connect a Button. True fetches one missing boundary patch. Identical patch
    requests are cached in memory, so a held Toggle does not contact Copernicus
    again.
SiteDataJson : item / str
    SiteDataJson output from MusselFlow Site Data. It supplies the coordinate,
    regional layer catalogue, and available timestamps.
domain : item / Rhino.Geometry.Curve
    Closed planar Rhino curve. Its bounding-box centre is georeferenced to the
    requested SiteData coordinate. Rhino document units are converted to metres.
variable : item / str
    current_speed, temperature, salinity, oxygen, chlorophyll,
    satellite_chlorophyll, phytoplankton_carbon, nitrate, phosphate, tsm,
    or turbidity. ``chlorophyll`` is the depth-aware 7 km regional model;
    ``satellite_chlorophyll`` is the surface-only 300 m Sentinel-3 OLCI
    open-water field.
timeIndex : item / object
    Prefer an integer frame index: 0 selects the first SiteData frame, 1 the
    second, and so on. For safe direct wiring, an ISO timestamp from
    SiteData.Times is also accepted and resolved to its corresponding index.
    Missing daily satellite observations search backward by up to 30 days.
depth : item / float
    Positive sampling depth in metres. Surface products ignore this value.
resolution : item / int
    Cells along the longer side, clamped to 4-128 in every output mode. A
    square domain at 100 creates 10,000 locally represented Rhino samples;
    changing this value never creates additional Copernicus requests.
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
exportFieldData : item / bool
    Controls only the FieldDataJson output. False returns "{}"; True serializes
    the canonical PLS/ML records. Sampling, values, resolution, and geometry stay
    identical.

Data distinction:
    The component samples processed Copernicus model or ocean-colour products.
    It does not download raw Sentinel multispectral bands and does not reproduce
    SNAP atmospheric correction, masking, custom indices, or a multispectral
    data cube.

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
    Canonical ML/PLS-ready records when exportFieldData is True; otherwise "{}".
Legend : list / str
    Variable, source units, raw range, display transform, time, and depth.
Report : list / str
    Timing, sample counts, API/cache state, resolution warnings, limits, and a
    chronological BACKTEST TRACE block. Copy that complete block when reporting
    a failed or suspicious field run.

SDK setup
---------
Create a Rhino 8 Python 3 component, convert the default component with
``Convert To GH_ScriptInstance``, then replace its generated text with this
complete file. The RunScript annotations create and type-hint the twelve inputs.
Add nine output sockets once in the exact order above. BeforeRunScript applies
all names and human-readable hover tooltips.

Optional native-patch backend
-----------------------------
The component never auto-installs ``copernicusmarine``. A compiled plugin may
bundle it as an internal dependency. During raw-script development it may be
installed into a dedicated Rhino Python environment, never the shared default
environment. If it is unavailable, the component stops before any field request
instead of falling back to hundreds of point requests.

The Toolbox bulk service requires a free Copernicus Marine account and a
one-time ``copernicusmarine.login`` for this isolated runtime. The component
checks for credentials before requesting data and returns ``AUTH_REQUIRED``
instead of opening an interactive prompt inside Grasshopper.

Backtest contract
-----------------
Every solution records the georeferenced boundary signature and WGS84 bounding
box, bulk ``read_dataframe`` arguments, cache key, returned native-grid shape,
local nearest-cell mapping distances, and a final PASS/FAIL verdict. The trace
counts component-level SDK calls; the SDK may perform its own transport details
internally. MusselFlow itself performs zero per-cell GetFeatureInfo requests.
"""

import datetime
import hashlib
import json
import math
import os
import sys
import time
import traceback

import Grasshopper
import Rhino
import System.Drawing


COMPONENT_BUILD = "2026-08-18f"
MAX_PREVIEW_RESOLUTION = 128
MAX_SAMPLES = MAX_PREVIEW_RESOLUTION*MAX_PREVIEW_RESOLUTION
HIGH_RESOLUTION_WARNING = 100
SATELLITE_LOOKBACK_DAYS = 30

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
        "source": "Northwest Shelf regional biogeochemical model",
        "nominal_resolution_m": 7000.0,
    },
    "satellite_chlorophyll": {
        "layer": "satellite_chlorophyll", "units": "ug/L",
        "kind": "scalar", "daily": True, "surface": True,
        "source": "Sentinel-3A/B OLCI open-ocean colour",
        "nominal_resolution_m": 300.0,
        "lookback_days": SATELLITE_LOOKBACK_DAYS,
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
        "lookback_days": SATELLITE_LOOKBACK_DAYS,
    },
    "turbidity": {
        "layer": "turbidity", "units": "FNU", "kind": "scalar",
        "daily": True, "surface": True,
        "lookback_days": SATELLITE_LOOKBACK_DAYS,
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
    "surface_chlorophyll": "satellite_chlorophyll",
    "chlorophyll_surface": "satellite_chlorophyll",
    "chlorophyll_depth_model": "chlorophyll",
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
        "Sample a Copernicus field over the WGS84 extent derived from the "
        "input Rhino domain, then map physical values into coloured geometry "
        "and optional PLS-ready records."),
}

INPUT_METADATA = (
    ("fetch", "fetch", "Button: fetch one missing boundary patch; identical patch requests use memory cache."),
    ("SiteDataJson", "SiteData", "Actual sampled SiteDataJson from MusselFlow Site Data."),
    ("domain", "domain", "Closed planar curve defining the sampled tile size and shape."),
    ("variable", "variable", "Field name. For chlorophyll: satellite_chlorophyll = 300 m surface Sentinel-3 OLCI; chlorophyll = depth-aware regional model."),
    ("timeIndex", "timeIndex", "Integer frame index or exact SiteData timestamp. Missing daily satellite data searches backward up to 30 days."),
    ("depth", "depth", "Positive depth in metres; ignored by surface satellite products."),
    ("resolution", "resolution", "Local Rhino preview cells along the longer side, clamped to 4-128. A square grid at 100 creates 10,000 local samples but no additional Copernicus requests."),
    ("sizePower", "sizePower", "Hotspot radius exponent; 2 emphasizes high values."),
    ("colors", "colors", "Optional ordered System.Drawing colour stops, low to high."),
    ("northVector", "northVector", "Geographic north in the Rhino domain plane; defaults to World +Y."),
    ("placementPoint", "placementPoint", "Optional Rhino point for the displayed tile centre; geospatial sampling is unchanged."),
    ("exportFieldData", "exportFieldData", "Controls only FieldDataJson serialization. It never changes sampling, values, resolution, or geometry."),
)

OUTPUT_METADATA = (
    ("FieldMesh", "FieldMesh", "Vertex-coloured grid cells for valid physical samples."),
    ("HotspotMesh", "HotspotMesh", "Filled value-scaled hotspot circles."),
    ("Circles", "Circles", "Hotspot circle outline curves."),
    ("Points", "Points", "Rhino sample points matching Values and Normalized."),
    ("Values", "Values", "Raw physical values in the Legend units."),
    ("Normalized", "Normalized", "Display-only values in the interval 0-1."),
    ("FieldDataJson", "FieldData", "Canonical records only when exportFieldData is True; otherwise {}."),
    ("Legend", "Legend", "Variable, units, range, timestamp, and display mapping."),
    ("Report", "Report", "Request/cache state, coverage, warnings, scientific scope, and chronological bulk-patch backtest trace."),
)

_PATCH_CACHE = {}


class FieldTrace(object):
    """Chronological, copyable diagnostics for one component solution."""

    def __init__(self):
        self.started = time.perf_counter()
        self.events = []
        self.sequence = 0

    @staticmethod
    def _text(value):
        try:
            return str(value).replace("\r", " ").replace("\n", " / ")
        except Exception:
            return "<unprintable diagnostic value>"

    def add(self, event, detail=""):
        self.sequence += 1
        elapsed = time.perf_counter()-self.started
        line = "TRACE %04d | +%.3fs | %s" % (
            self.sequence, elapsed, self._text(event))
        if detail:
            line += " | "+self._text(detail)
        self.events.append(line)

    def add_exception(self, stage, exception):
        self.add(
            str(stage)+" ERROR",
            "%s: %s" % (type(exception).__name__, self._text(exception)))
        formatted = traceback.format_exc(limit=12).strip()
        if formatted and formatted != "NoneType: None":
            for line in formatted.splitlines():
                if line.strip():
                    self.add("PYTHON", line.strip())

    def lines(self):
        return ([
            "BACKTEST TRACE BEGIN | copy every line through BACKTEST TRACE END",
        ]+list(self.events)+[
            "BACKTEST TRACE END | %d events | %.3fs elapsed"
            % (self.sequence, time.perf_counter()-self.started),
        ])


class CopernicusAuthenticationRequired(RuntimeError):
    pass


def apply_component_metadata(component):
    if component is None:
        return
    component.Name = COMPONENT_METADATA["name"]
    component.NickName = COMPONENT_METADATA["nickname"]
    component.Description = COMPONENT_METADATA["description"]
    component.Message = "Native patch"
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


def empty_outputs(status, message, trace=None):
    report = [
        "MUSSELFLOW SITE FIELD | build "+COMPONENT_BUILD+" | "+status,
        message,
    ]
    if trace is not None:
        report.extend(trace.lines())
    return (None, None, [], [], [], [], "{}", [], report)


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


def normalize_timestamp(value):
    """Normalize an ISO timestamp for exact frame matching."""
    text = "" if value is None else str(value).strip()
    if text.endswith("Z"):
        text = text[:-1]+"+00:00"
    return text


def resolve_frame(value, frames):
    """Return (index, selection mode) from an integer or ISO timestamp."""
    if value is None or not str(value).strip():
        return 0, "default index"

    # Keep integer indices as the canonical internal representation. Decimal
    # strings are accepted because Grasshopper Panels commonly provide text.
    text = str(value).strip()
    try:
        number = float(text)
        if math.isfinite(number) and number.is_integer():
            index = max(0, min(len(frames)-1, int(number)))
            return index, "frame index"
    except (TypeError, ValueError):
        pass

    requested = normalize_timestamp(text)
    for index, frame in enumerate(frames):
        candidate = normalize_timestamp(frame.get("time_utc"))
        if candidate == requested:
            return index, "exact timestamp"
    raise ValueError(
        "timeIndex must be an integer frame index or an exact timestamp from "
        "SiteData.Times. Received: %s" % text)


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

    # Georeference the actual boundary, not merely the preview cell centres.
    # Therefore changing `resolution` cannot change the downloaded source patch.
    boundary_geo = []
    for point in boundary:
        delta = point-anchor
        east_m = (delta*east_axis)*metres_per_unit
        north_m = (delta*north_axis)*metres_per_unit
        boundary_geo.append((
            anchor_latitude+north_m/111320.0,
            anchor_longitude+east_m/longitude_scale))
    field_bbox = {
        "minimum_latitude": min(item[0] for item in boundary_geo),
        "maximum_latitude": max(item[0] for item in boundary_geo),
        "minimum_longitude": min(item[1] for item in boundary_geo),
        "maximum_longitude": max(item[1] for item in boundary_geo),
    }
    boundary_signature = hashlib.sha256(";".join(
        "%.9f,%.9f,%.9f" % (point.X, point.Y, point.Z)
        for point in boundary).encode("utf-8")).hexdigest()[:16]
    domain_request = {
        "source": "input_closed_rhino_curve",
        "boundary_sample_count": len(boundary),
        "boundary_signature": boundary_signature,
        "rhino_plane_bounds": {
            "minimum_u": u_min,
            "maximum_u": u_max,
            "minimum_v": v_min,
            "maximum_v": v_max,
        },
        "metres_per_model_unit": metres_per_unit,
        "geographic_anchor_mode": "SiteData_coordinate_at_domain_bbox_centre",
        "geographic_anchor": {
            "latitude": anchor_latitude,
            "longitude": anchor_longitude,
        },
        "wgs84_bbox": dict(field_bbox),
    }
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
        width*metres_per_unit, height*metres_per_unit, anchor, field_bbox,
        domain_request)


def layer_contract(layer, kind):
    """Translate PRODUCT/DATASET/WMTS_VARIABLE into a Toolbox request."""
    parts = [part for part in str(layer or "").split("/") if part]
    if len(parts) < 3:
        raise ValueError("Copernicus layer is not PRODUCT/DATASET/VARIABLE.")
    wmts_dataset_id = parts[-2]
    dataset_id = wmts_dataset_id
    version_hint = None
    base, separator, suffix = wmts_dataset_id.rpartition("_")
    if (separator and len(suffix) == 6 and suffix.isdigit() and
            suffix.startswith("20")):
        # WMTS layer identifiers append the catalogue version (for example
        # `_202511`). The Toolbox expects the stable dataset ID and version as
        # separate arguments: dataset_id=... and dataset_version="202511".
        dataset_id = base
        version_hint = suffix
    wmts_variable = parts[-1]
    if kind == "vector_speed":
        variables = ["uo", "vo"]
    else:
        variables = [wmts_variable]
    return dataset_id, version_hint, variables


def patch_time_window(requested_timestamp, lookback_days):
    end = str(requested_timestamp)
    if int(lookback_days or 0) <= 0:
        return end, end
    return previous_daily_time(end, int(lookback_days)), end


def patch_cache_key(dataset_id, version_hint, variables, field_bbox,
                    start_time, end_time, depth):
    return (
        dataset_id,
        version_hint,
        tuple(variables),
        round(field_bbox["minimum_longitude"], 8),
        round(field_bbox["maximum_longitude"], 8),
        round(field_bbox["minimum_latitude"], 8),
        round(field_bbox["maximum_latitude"], 8),
        str(start_time),
        str(end_time),
        None if depth is None else round(float(depth), 4),
    )


def dataframe_column(frame, candidates):
    lookup = {str(column).lower(): column for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def iso_time(value):
    if value is None:
        return ""
    try:
        text = value.isoformat()
    except Exception:
        text = str(value)
    text = text.replace("+00:00", "Z")
    if text.endswith(".000000000"):
        text = text[:-10]
    return text


def timestamp_number(value):
    text = normalize_timestamp(value)
    if not text:
        return None
    try:
        return datetime.datetime.fromisoformat(text).timestamp()
    except Exception:
        try:
            return datetime.datetime.strptime(text[:19], "%Y-%m-%dT%H:%M:%S").timestamp()
        except Exception:
            return None


def dataframe_records(frame, variables, kind, requested_timestamp, depth):
    """Extract one latest valid native grid from a Toolbox DataFrame."""
    if frame is None or len(frame.index) == 0:
        return [], None
    table = frame.reset_index()
    latitude_column = dataframe_column(table, ("latitude", "lat"))
    longitude_column = dataframe_column(table, ("longitude", "lon"))
    time_column = dataframe_column(table, ("time", "datetime", "date"))
    depth_column = dataframe_column(table, ("depth", "elevation"))
    if latitude_column is None or longitude_column is None:
        raise ValueError("native patch has no latitude/longitude coordinates.")

    variable_columns = []
    for variable in variables:
        column = dataframe_column(table, (variable,))
        if column is None:
            raise ValueError("native patch is missing variable '%s'." % variable)
        variable_columns.append(column)

    rows = []
    for _, row in table.iterrows():
        latitude = finite_number(row[latitude_column])
        longitude = finite_number(row[longitude_column])
        if latitude is None or longitude is None:
            continue
        if kind == "vector_speed":
            east = finite_number(row[variable_columns[0]])
            north = finite_number(row[variable_columns[1]])
            value = None if east is None or north is None else math.hypot(east, north)
        else:
            value = finite_number(row[variable_columns[0]])
            if value is not None and kind == "oxygen":
                value *= 0.032
        if value is None:
            continue
        row_time = iso_time(row[time_column]) if time_column is not None else str(requested_timestamp)
        row_depth = finite_number(row[depth_column]) if depth_column is not None else None
        rows.append({
            "source_latitude": latitude,
            "source_longitude": longitude,
            "source_time": row_time,
            "source_depth": row_depth,
            "value": value,
        })
    if not rows:
        return [], None

    # A lookback patch may contain several dates. Use the newest available
    # observation not later than the requested frame.
    requested_number = timestamp_number(requested_timestamp)
    available_times = sorted(set(row["source_time"] for row in rows))
    eligible = [item for item in available_times
                if timestamp_number(item) is not None and
                (requested_number is None or timestamp_number(item) <= requested_number+1.0)]
    selected_time = eligible[-1] if eligible else available_times[-1]
    rows = [row for row in rows if row["source_time"] == selected_time]

    # A single requested depth uses the nearest native layer only.
    if depth is not None:
        source_depths = sorted(set(
            row["source_depth"] for row in rows if row["source_depth"] is not None))
        if source_depths:
            selected_depth = min(source_depths, key=lambda item: abs(item-float(depth)))
            rows = [row for row in rows
                    if row["source_depth"] is None or row["source_depth"] == selected_depth]

    for row in rows:
        row["source_cell_id"] = "%.8f:%.8f" % (
            row["source_latitude"], row["source_longitude"])
    return rows, selected_time


def load_copernicus_toolbox():
    """Import and configure the bulk backend after Script_Instance exists."""
    runtime_root = os.path.abspath(os.path.join(
        os.path.dirname(os.__file__), os.pardir, os.pardir))
    dependency_path = os.environ.get("MUSSELFLOW_COPERNICUS_PATH")
    if not dependency_path:
        dependency_path = os.path.join(
            runtime_root, "site-envs", "musselflow-copernicus")
    if os.path.isdir(dependency_path) and dependency_path not in sys.path:
        sys.path.insert(0, dependency_path)
    runtime_info = {
        "host_executable": sys.executable,
        "worker_executable": None,
        "dask_scheduler": "synchronous",
        "dependency_path": dependency_path,
    }
    # Embedded CPython reports Rhinoceros as sys.executable. Without an
    # explicit worker executable, multiprocessing.resource_tracker asks macOS
    # to open its `-c` payload as a Rhino document. Point any unavoidable
    # helper process at Rhino's actual CPython executable instead.
    worker_candidates = [
        os.path.join(runtime_root, "python3.9"),
        os.path.join(runtime_root, "python3"),
        os.path.join(runtime_root, "python.exe"),
    ]
    worker_executable = next(
        (candidate for candidate in worker_candidates
         if os.path.isfile(candidate)), None)
    try:
        multiprocessing = __import__("multiprocessing")
        if worker_executable:
            multiprocessing.set_executable(worker_executable)
            runtime_info["worker_executable"] = worker_executable
    except Exception as exception:
        runtime_info["multiprocessing_warning"] = "%s: %s" % (
            type(exception).__name__, str(exception))
    # read_dataframe is an in-memory subset operation. Force Dask to remain in
    # the Rhino process; this component does not need process workers.
    os.environ["DASK_SCHEDULER"] = "synchronous"
    try:
        dask = __import__("dask")
        dask.config.set(scheduler="synchronous")
        toolbox = __import__("copernicusmarine")
        if not callable(getattr(toolbox, "read_dataframe", None)):
            raise AttributeError(
                "copernicusmarine has no callable read_dataframe API")
        return toolbox, None, runtime_info
    except BaseException as exception:
        return None, "%s: %s | expected isolated path: %s" % (
            type(exception).__name__, str(exception), dependency_path), runtime_info


def copernicus_credentials_source():
    """Return a non-secret credential-source label or an actionable error."""
    username = os.environ.get("COPERNICUSMARINE_SERVICE_USERNAME")
    password = os.environ.get("COPERNICUSMARINE_SERVICE_PASSWORD")
    if username and password:
        return "environment variables", None
    if username or password:
        return None, (
            "Both COPERNICUSMARINE_SERVICE_USERNAME and "
            "COPERNICUSMARINE_SERVICE_PASSWORD must be set together.")

    home = os.path.expanduser("~")
    configured_root = os.environ.get(
        "COPERNICUSMARINE_CREDENTIALS_DIRECTORY")
    credentials_root = configured_root or home
    native_file = os.path.join(
        credentials_root, ".copernicusmarine",
        ".copernicusmarine-credentials")
    if os.path.isfile(native_file):
        return "Copernicus Marine credentials file", None
    motu_file = os.path.join(home, "motuclient", "motuclient-python.ini")
    if os.path.isfile(motu_file):
        return "legacy motuclient credentials file", None
    netrc_file = os.path.join(
        home, "_netrc" if os.name == "nt" else ".netrc")
    if os.path.isfile(netrc_file):
        try:
            netrc_module = __import__("netrc")
            configuration = netrc_module.netrc(netrc_file)
            for host in (
                    "auth.marine.copernicus.eu", "default_host",
                    "nrt.cmems-du.eu", "my.cmems-du.eu"):
                authentication = configuration.authenticators(host)
                if (authentication and authentication[0] and
                        authentication[2]):
                    return "netrc credentials", None
        except Exception:
            pass
    return None, (
        "No Copernicus Marine Toolbox credentials were found. Run "
        "copernicusmarine.login once with the isolated Rhino Python runtime, "
        "then restart Rhino.")


def fetch_native_patch(toolbox, layer, kind, field_bbox, requested_timestamp,
                       lookback_days, depth, allow_network, trace=None):
    dataset_id, version_hint, variables = layer_contract(layer, kind)
    start_time, end_time = patch_time_window(requested_timestamp, lookback_days)
    key = patch_cache_key(
        dataset_id, version_hint, variables, field_bbox,
        start_time, end_time, depth)
    cache_id = hashlib.sha256(repr(key).encode("utf-8")).hexdigest()[:16]
    if trace is not None:
        trace.add(
            "PATCH PLAN",
            "dataset=%s | wmts_version_hint=%s | variables=%s | "
            "time=%s..%s | depth=%s | cache_key=%s"
            % (dataset_id, version_hint or "none", ",".join(variables),
               start_time, end_time,
               "surface" if depth is None else "%.6g m" % float(depth),
               cache_id))
        trace.add(
            "PATCH BBOX",
            "west=%.8f | south=%.8f | east=%.8f | north=%.8f | selection=outside"
            % (field_bbox["minimum_longitude"],
               field_bbox["minimum_latitude"],
               field_bbox["maximum_longitude"],
               field_bbox["maximum_latitude"]))
    if key in _PATCH_CACHE:
        cached = dict(_PATCH_CACHE[key])
        cached["cache_hit"] = True
        cached["sdk_call_count"] = 0
        cached["sdk_seconds"] = 0.0
        cached["cache_id"] = cache_id
        if trace is not None:
            trace.add(
                "PATCH CACHE HIT",
                "cache_key=%s | records=%d | component SDK calls this run=0"
                % (cache_id, len(cached.get("records") or [])))
        return cached
    if not allow_network:
        if trace is not None:
            trace.add(
                "PATCH CACHE MISS",
                "cache_key=%s | fetch=False | network call suppressed" % cache_id)
        raise RuntimeError("native numerical patch is not cached; press fetch once.")
    if toolbox is None:
        raise RuntimeError(
            "optional Copernicus Toolbox bulk backend is unavailable.")
    credential_source, credential_error = copernicus_credentials_source()
    if credential_source is None:
        if trace is not None:
            trace.add("AUTH REQUIRED", credential_error)
        raise CopernicusAuthenticationRequired(credential_error)
    if trace is not None:
        trace.add("AUTH READY", "source="+credential_source)

    arguments = {
        "dataset_id": dataset_id,
        "variables": variables,
        "minimum_longitude": field_bbox["minimum_longitude"],
        "maximum_longitude": field_bbox["maximum_longitude"],
        "minimum_latitude": field_bbox["minimum_latitude"],
        "maximum_latitude": field_bbox["maximum_latitude"],
        "start_datetime": start_time,
        "end_datetime": end_time,
        "coordinates_selection_method": "outside",
        "disable_progress_bar": True,
    }
    if version_hint:
        # Preserve the exact WMTS/catalogue version without corrupting the
        # Toolbox dataset ID. This is the supported read_dataframe signature.
        arguments["dataset_version"] = version_hint
    if depth is not None:
        arguments["minimum_depth"] = float(depth)
        arguments["maximum_depth"] = float(depth)

    if trace is not None:
        trace.add(
            "PATCH CACHE MISS",
            "cache_key=%s | cache_entries_before=%d" % (
                cache_id, len(_PATCH_CACHE)))
        trace.add(
            "BULK SDK CALL START",
            "copernicusmarine.read_dataframe | component_call=1 | args=%s"
            % canonical_json(arguments))
    sdk_started = time.perf_counter()
    try:
        frame = toolbox.read_dataframe(**arguments)
    except Exception as exception:
        if trace is not None:
            trace.add(
                "BULK SDK CALL FAILED",
                "%.3fs | %s: %s" % (
                    time.perf_counter()-sdk_started,
                    type(exception).__name__, str(exception)))
        raise
    sdk_seconds = time.perf_counter()-sdk_started
    try:
        frame_rows = len(frame.index)
        frame_columns = [str(column) for column in frame.columns]
        index_names = [str(name) for name in frame.index.names]
    except Exception:
        frame_rows = -1
        frame_columns = []
        index_names = []
    if trace is not None:
        trace.add(
            "BULK SDK CALL END",
            "%.3fs | rows=%d | columns=%s | index=%s"
            % (sdk_seconds, frame_rows, ",".join(frame_columns),
               ",".join(index_names)))
    records, selected_time = dataframe_records(
        frame, variables, kind, requested_timestamp, depth)
    if not records:
        raise ValueError("native numerical patch contains no valid values.")
    latitudes = [record["source_latitude"] for record in records]
    longitudes = [record["source_longitude"] for record in records]
    values = [record["value"] for record in records]
    source_depths = sorted(set(
        record["source_depth"] for record in records
        if record.get("source_depth") is not None))
    if trace is not None:
        trace.add(
            "PATCH EXTRACTED",
            "records=%d | selected_time=%s | source_depths=%s | "
            "lon=%.8f..%.8f | lat=%.8f..%.8f | value=%.8g..%.8g"
            % (len(records), selected_time,
               ",".join("%.6g" % item for item in source_depths) or "surface",
               min(longitudes), max(longitudes),
               min(latitudes), max(latitudes), min(values), max(values)))
    result = {
        "dataset_id": dataset_id,
        "dataset_version": version_hint,
        "wmts_version_hint": version_hint,
        "variables": variables,
        "records": records,
        "selected_time": selected_time,
        "requested_bbox": dict(field_bbox),
        "cache_hit": False,
        "cache_id": cache_id,
        "sdk_call_count": 1,
        "sdk_seconds": sdk_seconds,
        "frame_row_count": frame_rows,
        "frame_columns": frame_columns,
        "source_depths": source_depths,
    }
    _PATCH_CACHE[key] = dict(result)
    if trace is not None:
        trace.add(
            "PATCH CACHE STORE",
            "cache_key=%s | cache_entries_after=%d" % (
                cache_id, len(_PATCH_CACHE)))
    return result


def apply_native_patch(samples, records, anchor_latitude):
    """Nearest-cell resampling is local; it performs no HTTP operations."""
    longitude_weight = max(1e-9, math.cos(math.radians(anchor_latitude)))
    distances_m = []
    source_cells = set()
    for sample in samples:
        latitude = sample["latitude"]
        longitude = sample["longitude"]
        nearest = min(records, key=lambda record: (
            (record["source_latitude"]-latitude)**2+
            ((record["source_longitude"]-longitude)*longitude_weight)**2))
        distance_degrees = math.sqrt(
            (nearest["source_latitude"]-latitude)**2+
            ((nearest["source_longitude"]-longitude)*longitude_weight)**2)
        distances_m.append(distance_degrees*111320.0)
        sample["value"] = nearest["value"]
        sample["source_latitude"] = nearest["source_latitude"]
        sample["source_longitude"] = nearest["source_longitude"]
        sample["source_time"] = nearest.get("source_time")
        sample["source_depth"] = nearest.get("source_depth")
        sample["source_cell_id"] = nearest.get("source_cell_id")
        source_cells.add(sample["source_cell_id"])
    return {
        "preview_samples": len(samples),
        "patch_records": len(records),
        "unique_source_cells": len(source_cells),
        "mean_nearest_distance_m": (
            sum(distances_m)/float(len(distances_m)) if distances_m else 0.0),
        "maximum_nearest_distance_m": max(distances_m) if distances_m else 0.0,
        "per_cell_http_calls": 0,
    }


def ocean_colour_layer(layers, timestamp):
    """Choose OLCI MY for archives and NRT only for the recent rolling window."""
    archive = layers.get("satellite_chlorophyll")
    recent = layers.get("satellite_chlorophyll_nrt")
    if not archive:
        return None
    if not recent:
        return archive
    try:
        requested = datetime.datetime.strptime(
            str(timestamp)[:10], "%Y-%m-%d").date()
        age_days = (datetime.datetime.utcnow().date()-requested).days
    except Exception:
        return archive
    return recent if -2 <= age_days <= 35 else archive


def daily_time(timestamp):
    return str(timestamp)[:10]+"T00:00:00Z"


def previous_daily_time(timestamp, days_back):
    """Return a UTC midnight timestamp a whole number of days earlier."""
    try:
        day = datetime.datetime.strptime(str(timestamp)[:10], "%Y-%m-%d")
    except (TypeError, ValueError):
        raise ValueError("cannot derive a daily fallback from time: %s" % timestamp)
    day -= datetime.timedelta(days=int(days_back))
    return day.strftime("%Y-%m-%dT00:00:00Z")


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
            SiteDataJson: str,
            domain: Rhino.Geometry.Curve,
            variable: str,
            timeIndex: object,
            depth: float,
            resolution: int,
            sizePower: float,
            colors: list[System.Drawing.Color],
            northVector: Rhino.Geometry.Vector3d,
            placementPoint: Rhino.Geometry.Point3d,
            exportFieldData: bool):
        """Preview one field quickly; optionally export canonical analysis records."""
        started = time.perf_counter()
        trace = FieldTrace()
        trace.add(
            "RUN START",
            "build=%s | fetch=%s | resolution_input=%s | exportFieldData=%s"
            % (COMPONENT_BUILD, bool(fetch), resolution, bool(exportFieldData)))
        try:
            site, anchor_lat, anchor_lon, frames = parse_site_data(SiteDataJson)
            selected_variable = variable_name(variable)
            if selected_variable not in VARIABLES:
                raise ValueError(
                    "unknown variable '%s'. Choose: %s" % (
                        selected_variable, ", ".join(sorted(VARIABLES))))
            specification = VARIABLES[selected_variable]
            export_enabled = bool(exportFieldData)
            index, time_selection_mode = resolve_frame(timeIndex, frames)
            requested_timestamp = str(frames[index].get("time_utc") or "")
            if not requested_timestamp:
                raise ValueError("selected SiteData frame has no timestamp.")
            if specification["daily"]:
                requested_timestamp = daily_time(requested_timestamp)
            if selected_variable == "satellite_chlorophyll":
                layer = ocean_colour_layer(site["layers"], requested_timestamp)
            else:
                layer = site["layers"].get(specification["layer"])
            if not layer:
                raise ValueError(
                    "siteData regional catalogue has no '%s' layer."
                    % specification["layer"])
            sample_depth = abs(finite_number(depth) or 0.0)
            requested_resolution = 12 if resolution is None else int(resolution)
            count = max(4, min(MAX_PREVIEW_RESOLUTION, requested_resolution))
            if count != requested_resolution:
                trace.add(
                    "RESOLUTION CLAMP",
                    "requested=%d | applied=%d | allowed=4..%d"
                    % (requested_resolution, count, MAX_PREVIEW_RESOLUTION))
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
            (samples, nx, ny, du, dv, width_m, height_m, domain_anchor,
             field_bbox, domain_request) = grid_samples(
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
            trace.add(
                "INPUT CONTRACT",
                "variable=%s | frame=%d/%d | selection=%s | requested_time=%s | "
                "depth=%s | preview_resolution=%d"
                % (selected_variable, index, len(frames), time_selection_mode,
                   requested_timestamp,
                   "surface" if specification["surface"]
                   else "%.6g m" % sample_depth,
                   count))
            trace.add(
                "LAYER CONTRACT",
                "catalogue_layer=%s | kind=%s | lookback_days=%d"
                % (layer, specification["kind"],
                   int(specification.get("lookback_days") or 0)))
            trace.add(
                "BOUNDARY SAMPLE",
                "source=closed Rhino curve | points=%d | signature=%s | "
                "model_size=%.3f x %.3f m | metres_per_unit=%.12g"
                % (domain_request["boundary_sample_count"],
                   domain_request["boundary_signature"], width_m, height_m,
                   metres_per_unit))
            trace.add(
                "PREVIEW GRID",
                "columns=%d | rows=%d | candidate_cells=%d | inside_boundary=%d"
                % (nx, ny, nx*ny, len(samples)))
            trace.add(
                "GEOREFERENCE",
                "anchor_lat=%.8f | anchor_lon=%.8f | bbox=%.8f,%.8f,%.8f,%.8f"
                % (anchor_lat, anchor_lon,
                   field_bbox["minimum_longitude"],
                   field_bbox["minimum_latitude"],
                   field_bbox["maximum_longitude"],
                   field_bbox["maximum_latitude"]))
        except Exception as exception:
            trace.add_exception("INPUT/BOUNDARY", exception)
            return empty_outputs("INVALID_INPUT", str(exception), trace)

        query_depth = None if specification["surface"] else sample_depth
        lookback_days = int(specification.get("lookback_days") or 0)
        timestamp = None
        fallback_days = None
        missing_before = 0
        cached_before = 0
        sdk_call_count = 0
        backend = "COPERNICUS TOOLBOX NATIVE PATCH"

        trace.add(
            "BACKEND IMPORT START",
            "lazy import after GH_ScriptInstance construction")
        toolbox, toolbox_error, runtime_info = load_copernicus_toolbox()
        trace.add(
            "EMBEDDED RUNTIME",
            "host_executable=%s | worker_executable=%s | dask_scheduler=%s"
            % (runtime_info.get("host_executable"),
               runtime_info.get("worker_executable") or "unavailable",
               runtime_info.get("dask_scheduler")))
        if runtime_info.get("multiprocessing_warning"):
            trace.add(
                "RUNTIME WARNING",
                runtime_info["multiprocessing_warning"])
        if toolbox is None:
            trace.add("BACKEND IMPORT FAILED", toolbox_error)
            return empty_outputs(
                "BULK_BACKEND_UNAVAILABLE",
                "The optional Copernicus Toolbox is not available (%s). No "
                "package installation and no per-cell WMTS fallback were attempted."
                % toolbox_error, trace)
        trace.add(
            "BACKEND READY",
            "copernicusmarine=%s | read_dataframe=callable | module=%s"
            % (getattr(toolbox, "__version__", "unknown"),
               getattr(toolbox, "__file__", "unknown")))

        field_started = time.perf_counter()
        try:
            patch = fetch_native_patch(
                toolbox, layer, specification["kind"], field_bbox,
                requested_timestamp, lookback_days, query_depth, bool(fetch),
                trace=trace)
            mapping_stats = apply_native_patch(
                samples, patch["records"], anchor_lat)
            timestamp = patch.get("selected_time") or requested_timestamp
            sdk_call_count = int(patch.get("sdk_call_count") or 0)
            if patch.get("cache_hit"):
                cached_before = 1
            else:
                missing_before = 1
            trace.add(
                "LOCAL RESAMPLE",
                "preview_samples=%d | patch_records=%d | unique_source_cells=%d | "
                "mean_nearest=%.3f m | max_nearest=%.3f m | per_cell_http_calls=%d"
                % (mapping_stats["preview_samples"],
                   mapping_stats["patch_records"],
                   mapping_stats["unique_source_cells"],
                   mapping_stats["mean_nearest_distance_m"],
                   mapping_stats["maximum_nearest_distance_m"],
                   mapping_stats["per_cell_http_calls"]))
            requested_number = timestamp_number(requested_timestamp)
            selected_number = timestamp_number(timestamp)
            if requested_number is not None and selected_number is not None:
                fallback_days = max(
                    0, int(round((requested_number-selected_number)/86400.0)))
        except Exception as patch_exception:
            trace.add_exception("BULK PATCH", patch_exception)
            if isinstance(
                    patch_exception, CopernicusAuthenticationRequired):
                return empty_outputs(
                    "AUTH_REQUIRED", str(patch_exception), trace)
            if not fetch:
                return empty_outputs(
                    "WAITING",
                    "Press fetch once for the native numerical bounding-box patch. "
                    +str(patch_exception), trace)
            return empty_outputs(
                "BULK_FETCH_ERROR",
                "%s: %s. No per-cell WMTS fallback was attempted."
                % (type(patch_exception).__name__, str(patch_exception)), trace)
        field_seconds = time.perf_counter()-field_started

        try:
            valid_samples, low, high = normalize_samples(samples)
        except Exception as exception:
            trace.add_exception("NORMALIZE", exception)
            return empty_outputs(
                "NORMALIZE_ERROR", "%s: %s" % (
                    type(exception).__name__, str(exception)), trace)
        if not valid_samples:
            trace.add(
                "NO VALID DATA",
                "patch_records=%d | preview_samples=%d"
                % (len(patch.get("records") or []), len(samples)))
            return empty_outputs(
                "NO_DATA",
                "The bulk patch returned no valid %s values for this "
                "domain/time/depth." % selected_variable, trace)
        trace.add(
            "VALUES NORMALIZED",
            "valid=%d/%d | raw_min=%.8g | raw_max=%.8g"
            % (len(valid_samples), len(samples), low, high))

        geometry_started = time.perf_counter()
        try:
            stops = color_stops(colors)
            field_mesh, hotspot_mesh, circles, points, values, normalized = (
                make_geometry(
                    valid_samples, display_plane, du, dv, power, stops))
        except Exception as exception:
            trace.add_exception("GEOMETRY", exception)
            return empty_outputs(
                "GEOMETRY_ERROR", "%s: %s" % (
                    type(exception).__name__, str(exception)), trace)
        geometry_seconds = time.perf_counter()-geometry_started
        trace.add(
            "GEOMETRY READY",
            "field_cells=%d | hotspot_circles=%d | points=%d | %.3fs"
            % (len(valid_samples), len(circles), len(points), geometry_seconds))
        unique_cells = source_cell_count(valid_samples)
        region = site.get("regional_product") or {}
        json_started = time.perf_counter()
        if export_enabled:
            field_document = {
                "schema": "musselflow.site_field.1.0.0",
                "build": COMPONENT_BUILD,
                "variable": selected_variable,
                "units": specification["units"],
                "data_source": specification.get("source", "Copernicus regional product"),
                "nominal_resolution_m": specification.get("nominal_resolution_m"),
                "time_utc": timestamp,
                "requested_time_utc": requested_timestamp,
                "time_fallback_days": fallback_days,
                "source_frame_index": index,
                "time_selection_mode": time_selection_mode,
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
                "source_patch": {
                    "backend": backend,
                    "dataset_id": patch.get("dataset_id"),
                    "wmts_version_hint": patch.get("wmts_version_hint"),
                    "wgs84_bbox": field_bbox,
                    "component_sdk_call_count": sdk_call_count,
                    "per_cell_http_call_count": 0,
                    "cache_hit": bool(patch.get("cache_hit")),
                    "cache_id": patch.get("cache_id"),
                    "native_record_count": len(patch.get("records") or []),
                    "mapping": mapping_stats,
                },
                "domain_request": domain_request,
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
                        "source_time": sample.get("source_time", timestamp),
                        "source_depth": sample.get("source_depth"),
                        "source_cell_id": sample.get("source_cell_id"),
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
    
            try:
                field_data_json = canonical_json(field_document)
            except Exception as exception:
                trace.add_exception("FIELD JSON", exception)
                return empty_outputs(
                    "FIELD_JSON_ERROR", "%s: %s" % (
                        type(exception).__name__, str(exception)), trace)
        else:
            field_data_json = "{}"
        json_seconds = time.perf_counter()-json_started
        trace.add(
            "FIELD JSON",
            "enabled=%s | characters=%d | %.3fs"
            % (export_enabled, len(field_data_json), json_seconds))
        legend = [
            "VARIABLE | %s" % selected_variable,
            "SOURCE | %s" % specification.get(
                "source", "Copernicus regional product"),
            "NOMINAL RESOLUTION | %s" % (
                "not declared" if specification.get("nominal_resolution_m") is None
                else "%.6g m" % specification["nominal_resolution_m"]),
            "UNITS | %s" % specification["units"],
            "RAW RANGE | %.6g to %.6g" % (low, high),
            "DISPLAY | linear colour; hotspot radius = normalized^%.3g" % power,
            "TIME | %s" % timestamp,
            "REQUESTED TIME | %s" % requested_timestamp,
            "DEPTH | %s" % (
                "surface product" if specification["surface"]
                else "%.6g m" % sample_depth),
        ]
        report = [
            "MUSSELFLOW SITE FIELD | build %s | %s | %d/%d valid | %.3fs"
            % (COMPONENT_BUILD, selected_variable, len(valid_samples),
               len(samples), time.perf_counter()-started),
            "MODE | %s"
            % ("VISUAL + FIELD JSON" if export_enabled else "VISUAL ONLY"),
            "DATA CONTRACT | SiteDataJson input active | FieldDataJson output %s"
            % ("enabled" if export_enabled else "disabled"),
            "TIMING | field %.3fs | geometry %.3fs | json %.3fs"
            % (field_seconds, geometry_seconds, json_seconds),
            "DATA BACKEND | %s" % backend,
            "DATASET | %s | WMTS version hint %s"
            % (patch.get("dataset_id"),
               patch.get("wmts_version_hint") or "none"),
            "REQUEST PLAN | one boundary patch | read_dataframe calls this run %d | "
            "per-cell GetFeatureInfo calls 0"
            % sdk_call_count,
            "REQUEST REDUCTION | %d bulk SDK call(s) instead of up to %d legacy "
            "point requests for this preview"
            % (sdk_call_count, len(samples)),
            "WGS84 BBOX | %.8f, %.8f to %.8f, %.8f"
            % (field_bbox["minimum_longitude"], field_bbox["minimum_latitude"],
               field_bbox["maximum_longitude"], field_bbox["maximum_latitude"]),
            "DOMAIN BOUNDARY | input closed curve | %d sampled boundary points | signature %s"
            % (domain_request["boundary_sample_count"],
               domain_request["boundary_signature"]),
            "GEOREFERENCE | SiteData %.8f, %.8f is fixed at the input domain bounding-box centre"
            % (anchor_lat, anchor_lon),
            "REGION | %s" % (region.get("name") or "from SiteData catalogue"),
            "DOMAIN | %.3f x %.3f m | grid %d x %d | %d inside"
            % (width_m, height_m, nx, ny, len(samples)),
            "PREVIEW RESOLUTION | requested %d | applied %d | maximum %d"
            % (requested_resolution, count, MAX_PREVIEW_RESOLUTION),
            "PLACEMENT | Rhino %.3f, %.3f, %.3f | SiteData WGS84 unchanged"
            % (display_anchor.X, display_anchor.Y, display_anchor.Z),
            "CACHE | %d uncached patch request(s) | %d cached"
            % (missing_before, cached_before),
            "PATCH CACHE KEY | %s | independent of preview resolution"
            % patch.get("cache_id"),
            "SOURCE CELLS | %d distinct Copernicus cells represented"
            % unique_cells,
            "LOCAL MAPPING | %d preview samples from %d native records | "
            "nearest mean %.3f m | max %.3f m"
            % (mapping_stats["preview_samples"], mapping_stats["patch_records"],
               mapping_stats["mean_nearest_distance_m"],
               mapping_stats["maximum_nearest_distance_m"]),
            "TIME SELECTION | %s | frame %d | %s"
            % (time_selection_mode, index, requested_timestamp),
        ]
        if fallback_days:
            report.append(
                "TIME FALLBACK | requested %s | used %s | %d day(s) earlier"
                % (requested_timestamp, timestamp, fallback_days))
        elif lookback_days:
            report.append(
                "TIME FALLBACK | exact requested daily observation available")
        if query_depth is not None:
            native_depths = patch.get("source_depths") or []
            report.append(
                "DEPTH SELECTION | requested %.6g m | nearest native layer %s"
                % (query_depth,
                   ", ".join("%.6g m" % item for item in native_depths)
                   if native_depths else "not reported"))
        if selected_variable == "chlorophyll":
            report.append(
                "SOURCE MODE | depth-aware regional BGC model | nominal 7 km; "
                "use satellite_chlorophyll for the 300 m Sentinel-3 OLCI surface field.")
        elif selected_variable == "satellite_chlorophyll":
            report.append(
                "SOURCE MODE | Sentinel-3A/B OLCI open-ocean colour | nominal 300 m | "
                "surface only; depth input ignored.")
        if unit_warning:
            report.append("UNIT WARNING | "+unit_warning)
        if unique_cells <= 1 and len(valid_samples) > 1:
            report.append(
                "RESOLUTION WARNING | the Rhino domain is smaller than this "
                "product's spatial support; all samples resolve to one source cell.")
        if count >= HIGH_RESOLUTION_WARNING:
            report.append(
                "PERFORMANCE NOTE | resolution %d produced %d represented cells; "
                "Rhino mesh, circle, and JSON work grows approximately with "
                "resolution squared. Copernicus requests remain one-or-zero."
                % (count, len(valid_samples)))
        if abs(high-low) <= 1e-15:
            report.append(
                "DISPLAY NOTE | all valid samples contain one physical source value; "
                "Normalized is 0.5 by definition, not a missing-data placeholder.")
        if max(width_m, height_m) > 100000.0:
            report.append(
                "GEOREFERENCE WARNING | domain exceeds 100 km; the local tangent "
                "coordinate approximation should be replaced by a map projection.")
        report.extend([
            "RESOLUTION NOTE | higher preview resolution adds local Rhino samples; "
            "it does not increase the native Copernicus product resolution.",
            "DATA LIMIT | processed model/ocean-colour product, not raw multispectral bands.",
            "VALIDATION LIMIT | regional fields require product-QC and local observations.",
        ])
        invariant_failures = []
        if sdk_call_count not in (0, 1):
            invariant_failures.append("read_dataframe call count is not 0 or 1")
        if mapping_stats["per_cell_http_calls"] != 0:
            invariant_failures.append("per-cell HTTP calls were detected")
        if mapping_stats["preview_samples"] != len(samples):
            invariant_failures.append("not every preview sample was mapped")
        if patch.get("requested_bbox") != field_bbox:
            invariant_failures.append("returned patch key does not match boundary bbox")
        dataset_tail = str(patch.get("dataset_id") or "").rpartition("_")[2]
        if (len(dataset_tail) == 6 and dataset_tail.isdigit() and
                dataset_tail.startswith("20")):
            invariant_failures.append(
                "Toolbox dataset_id still contains a WMTS version suffix")
        if invariant_failures:
            verdict = "FAIL | "+"; ".join(invariant_failures)
        else:
            verdict = (
                "PASS | one-or-zero bulk SDK calls, zero per-cell calls, "
                "boundary bbox preserved, all preview samples mapped locally")
        report.append("BACKTEST VERDICT | "+verdict)
        trace.add("BACKTEST VERDICT", verdict)
        trace.add(
            "REQUEST REDUCTION",
            "bulk_sdk_calls=%d | legacy_point_call_upper_bound=%d | "
            "actual_per_cell_calls=0"
            % (sdk_call_count, len(samples)))
        trace.add(
            "RESOLUTION SEMANTICS",
            "preview_resolution=%d changes local sample density only | "
            "boundary cache_key=%s"
            % (count, patch.get("cache_id")))
        trace.add(
            "RUN COMPLETE",
            "valid=%d/%d | source_cells=%d | total=%.3fs"
            % (len(valid_samples), len(samples), unique_cells,
               time.perf_counter()-started))
        report.extend(trace.lines())
        return (
            field_mesh, hotspot_mesh, circles, points, values, normalized,
            field_data_json, legend, report)

    def AfterRunScript(self):
        pass
