"""
MusselFlow Site Data — Copernicus Regional Router
==================================================

Self-contained Rhino 8 Grasshopper Python 3 SDK-mode component. It samples
Copernicus Marine WMTS JSON at one WGS84 site and compiles the result into the
existing MusselFlow ecological case. It uses only Python's standard library,
RhinoCommon, and Grasshopper: no external executable, sidecar process, desktop
application, account credentials, NumPy, or operating-system-specific path.

This implementation routes Baltic coordinates to Baltic products and North Sea
coordinates to Northwest Shelf products. The HTTP backend is isolated in small
functions so a compiled Rhino plugin can reuse the same data contract on both
Windows and macOS.

Inputs
------
fetch : item / bool
    Connect a Button. True performs the request. Identical requests are cached
    in memory so a held Toggle cannot repeatedly contact Copernicus.
latitude : item / float
    Site latitude in WGS84 decimal degrees.
longitude : item / float
    Site longitude in WGS84 decimal degrees.
startTime : item / str
    First UTC time, for example ``2026-04-20T00:00:00Z``.
endTime : item / str
    Last UTC time. Leave empty to sample 24 hours after startTime.
depths : list / float
    Positive water depths in metres, such as ``[1, 5, 10]``. Values are
    sampled at Copernicus' nearest available model layers and averaged equally.
frameCount : item / int
    Number of chronological states, clamped to 1–10. Six is recommended.
BaseModelJson : item / str
    Complete existing ``musselflow_grammar.json`` text. The component replaces
    only its forcing timeline and initial oxygen; biological constants remain.
northVector : item / Rhino.Geometry.Vector3d
    Optional Rhino-plan vector for geographic north. Defaults to World +Y.
FetchRequestJson : item / str
    Optional request from Copernicus Data Browser. When connected, its site,
    time and depths replace the corresponding manual inputs. This makes the
    browser-to-fetcher connection one JSON wire.

Outputs — create seven ports once in this exact order
-----------------------------------------------------
FlowVectors : list / Rhino.Geometry.Vector3d
    Depth-averaged ambient currents in m/s, ready for MusselFlow.flowVectors.
Times : list / str
    UTC timestamps corresponding one-to-one with FlowVectors.
SiteDataJson : item / str
    Canonical JSON containing profiles, converted values, provenance, and gaps.
SimulationCaseJson : item / str
    BaseModelJson with a site-derived forcing timeline, ready for MusselFlow.
Values : list / str
    Human-readable per-frame temperature, salinity, oxygen, chlorophyll, TSM.
Report : list / str
    Coverage, requested/actual coordinates, missing values, and limitations.
FetchLog : list / str
    Chronological trace of routing, cache decisions, HTTP URLs, request timing,
    response sizes, data assembly, and case compilation for the current run.

SDK setup
---------
Create a Rhino 8 Python 3 component, convert the default component with
``Convert To GH_ScriptInstance``, and only then replace its generated text with
this complete file. The typed RunScript signature creates and type-hints all
ten inputs, including List access for depths. Grasshopper cannot derive seven
output names from a returned tuple: add seven outputs once in the order above.
BeforeRunScript then applies every human-readable name and hover tooltip.

Scientific scope
----------------
Copernicus values are ambient model/satellite boundary conditions, not
farm-scale CFD and not site measurements. The regional model cannot resolve
one-metre socks. Satellite TSM is a surface, weather-dependent observation and
may be unavailable for the requested day. Every fallback is reported.
"""

import concurrent.futures
import copy
import datetime
import hashlib
import json
import math
import threading
import time
import urllib.parse
import urllib.request

import Grasshopper
import Rhino


COMPONENT_BUILD = "2026-08-08c"
SITE_REQUEST_SCHEMA = "copernicus.site_request.1.0"
WMTS_ENDPOINT = "https://wmts.marine.copernicus.eu/teroWmts/"
TILE_LEVEL = 10
HTTP_TIMEOUT_S = 25
MAX_WORKERS = 8

# Current layer versions are discovered from the public WMTS capabilities.
# Keep these in one table so a future plugin can update catalogue metadata
# without changing the component's data contract.
BALTIC_LAYERS = {
    "current": (
        "BALTICSEA_ANALYSISFORECAST_PHY_003_006/"
        "cmems_mod_bal_phy_anfc_PT1H-i_202411/sea_water_velocity"),
    "temperature": (
        "BALTICSEA_ANALYSISFORECAST_PHY_003_006/"
        "cmems_mod_bal_phy_anfc_PT1H-i_202411/thetao"),
    "salinity": (
        "BALTICSEA_ANALYSISFORECAST_PHY_003_006/"
        "cmems_mod_bal_phy_anfc_PT1H-i_202411/so"),
    "oxygen": (
        "BALTICSEA_ANALYSISFORECAST_BGC_003_007/"
        "cmems_mod_bal_bgc_anfc_P1D-m_202411/o2"),
    "chlorophyll": (
        "BALTICSEA_ANALYSISFORECAST_BGC_003_007/"
        "cmems_mod_bal_bgc_anfc_P1D-m_202411/chl"),
    "phytoplankton_carbon": (
        "BALTICSEA_ANALYSISFORECAST_BGC_003_007/"
        "cmems_mod_bal_bgc_anfc_P1D-m_202411/phyc"),
    "nitrate": (
        "BALTICSEA_ANALYSISFORECAST_BGC_003_007/"
        "cmems_mod_bal_bgc_anfc_P1D-m_202411/no3"),
    "phosphate": (
        "BALTICSEA_ANALYSISFORECAST_BGC_003_007/"
        "cmems_mod_bal_bgc_anfc_P1D-m_202411/po4"),
    "tsm": (
        "OCEANCOLOUR_BAL_BGC_HR_L4_NRT_009_208/"
        "cmems_obs_oc_bal_bgc_tur-spm-chl_nrt_l4-hr-mosaic_"
        "P1D-m_202107/SPM"),
    "satellite_chlorophyll": (
        "OCEANCOLOUR_BAL_BGC_L3_MY_009_133/"
        "cmems_obs-oc_bal_bgc-plankton_my_l3-olci-300m_P1D-m/CHL"),
    "satellite_chlorophyll_nrt": (
        "OCEANCOLOUR_BAL_BGC_L3_NRT_009_131/"
        "cmems_obs-oc_bal_bgc-plankton_nrt_l3-olci-300m_P1D/CHL"),
    "turbidity": (
        "OCEANCOLOUR_BAL_BGC_HR_L4_NRT_009_208/"
        "cmems_obs_oc_bal_bgc_tur-spm-chl_nrt_l4-hr-mosaic_"
        "P1D-m_202107/TUR"),
}

NORTHWEST_SHELF_LAYERS = {
    "current": (
        "NWSHELF_ANALYSISFORECAST_PHY_004_013/"
        "cmems_mod_nws_phy-cur_anfc_1.5km-3D_PT1H-i_202511/"
        "sea_water_velocity"),
    "temperature": (
        "NWSHELF_ANALYSISFORECAST_PHY_004_013/"
        "cmems_mod_nws_phy-tem_anfc_1.5km-3D_PT1H-i_202511/thetao"),
    "salinity": (
        "NWSHELF_ANALYSISFORECAST_PHY_004_013/"
        "cmems_mod_nws_phy-sal_anfc_1.5km-3D_PT1H-i_202511/so"),
    "oxygen": (
        "NWSHELF_ANALYSISFORECAST_BGC_004_002/"
        "cmems_mod_nws_bgc-o2_anfc_7km-3D_P1D-m_202511/o2"),
    "chlorophyll": (
        "NWSHELF_ANALYSISFORECAST_BGC_004_002/"
        "cmems_mod_nws_bgc-chl_anfc_7km-3D_P1D-m_202511/chl"),
    "phytoplankton_carbon": (
        "NWSHELF_ANALYSISFORECAST_BGC_004_002/"
        "cmems_mod_nws_bgc-phyc_anfc_7km-3D_P1D-m_202511/phyc"),
    "nitrate": (
        "NWSHELF_ANALYSISFORECAST_BGC_004_002/"
        "cmems_mod_nws_bgc-no3_anfc_7km-3D_P1D-m_202511/no3"),
    "phosphate": (
        "NWSHELF_ANALYSISFORECAST_BGC_004_002/"
        "cmems_mod_nws_bgc-po4_anfc_7km-3D_P1D-m_202511/po4"),
    "tsm": (
        "OCEANCOLOUR_NWS_BGC_HR_L3_NRT_009_203/"
        "cmems_obs_oc_nws_bgc_tur-spm-chl_nrt_l3-hr-mosaic_"
        "P1D-m_202107/SPM"),
    "satellite_chlorophyll": (
        "OCEANCOLOUR_ATL_BGC_L3_MY_009_113/"
        "cmems_obs-oc_atl_bgc-plankton_my_l3-olci-300m_P1D/CHL"),
    "satellite_chlorophyll_nrt": (
        "OCEANCOLOUR_ATL_BGC_L3_NRT_009_111/"
        "cmems_obs-oc_atl_bgc-plankton_nrt_l3-olci-300m_P1D/CHL"),
    "turbidity": (
        "OCEANCOLOUR_NWS_BGC_HR_L3_NRT_009_203/"
        "cmems_obs_oc_nws_bgc_tur-spm-chl_nrt_l3-hr-mosaic_"
        "P1D-m_202107/TUR"),
}

PRODUCTS = (
    {
        "id": "baltic",
        "name": "Baltic Sea",
        "bounds": (53.0, 66.0, 9.0, 31.0),
        "layers": BALTIC_LAYERS,
    },
    {
        "id": "northwest_shelf",
        "name": "Northwest European Shelf",
        "bounds": (46.0, 63.0, -17.0, 13.0),
        "layers": NORTHWEST_SHELF_LAYERS,
    },
)

COMPONENT_METADATA = {
    "name": "MusselFlow Site Data",
    "nickname": "SiteData",
    "description": (
        "Fetch regional Copernicus current and ecological boundary conditions "
        "for one WGS84 site; cache them and compile a MusselFlow timeline."),
}

INPUT_METADATA = (
    ("fetch", "fetch", "Button: fetch once; identical requests use memory cache."),
    ("latitude", "latitude", "WGS84 latitude in decimal degrees."),
    ("longitude", "longitude", "WGS84 longitude in decimal degrees."),
    ("startTime", "startTime", "First UTC timestamp, e.g. 2026-04-20T00:00:00Z."),
    ("endTime", "endTime", "Last UTC timestamp; empty means start + 24 hours."),
    ("depths", "depths", "Positive depth list in metres; nearest model layers are used."),
    ("frameCount", "frameCount", "Chronological states, 1–10; six recommended."),
    ("BaseModelJson", "BaseModel", "Ecological model, coefficients, objectives and constraints before site forcing."),
    ("northVector", "northVector", "Optional Rhino-plan vector for geographic north."),
    ("FetchRequestJson", "FetchRequest", "Not-yet-fetched request from Data Browser; supplies site, time, depths and selected product."),
)

OUTPUT_METADATA = (
    ("FlowVectors", "FlowVectors", "Depth-averaged current vectors in m/s."),
    ("Times", "Times", "UTC timestamps matching FlowVectors."),
    ("SiteDataJson", "SiteData", "Actual sampled environmental frames, units, layers, gaps and provenance."),
    ("SimulationCaseJson", "SimulationCase", "Base model patched with site forcing; connect to Optimizer.SimulationCaseJson."),
    ("Values", "Values", "Readable environmental values for each frame."),
    ("Report", "Report", "Coverage, fallbacks, warnings, and scientific limits."),
    ("FetchLog", "FetchLog", "Chronological cache, routing, HTTP request, response, and assembly trace."),
)

_MEMORY_CACHE = {}


class FetchTrace(object):
    """Thread-safe chronological trace shared by concurrent HTTP workers."""

    def __init__(self):
        self.started = time.perf_counter()
        self.events = []
        self.sequence = 0
        self.request_sequence = 0
        self.lock = threading.Lock()

    def add(self, event, detail=""):
        with self.lock:
            self.sequence += 1
            elapsed = time.perf_counter()-self.started
            line = "%04d | +%.3fs | %s" % (
                self.sequence, elapsed, str(event))
            if detail:
                line += " | "+str(detail)
            self.events.append(line)

    def request_id(self):
        with self.lock:
            self.request_sequence += 1
            return self.request_sequence

    def lines(self):
        with self.lock:
            return list(self.events)


def apply_component_metadata(component):
    """Apply component, port, and hover descriptions without changing wires."""
    if component is None:
        return
    component.Name = COMPONENT_METADATA["name"]
    component.NickName = COMPONENT_METADATA["nickname"]
    component.Description = COMPONENT_METADATA["description"]
    component.Message = "Regional WMTS"
    for index, (name, nickname, description) in enumerate(INPUT_METADATA):
        if index >= component.Params.Input.Count:
            break
        parameter = component.Params.Input[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description
        parameter.Optional = index != 0
        if index == 5:
            # SDK signature synchronization is not reliable after every
            # copy/paste into an existing component. Enforce depths as one
            # Grasshopper list so RunScript executes once for the full profile.
            parameter.Access = Grasshopper.Kernel.GH_ParamAccess.list
    for index, (name, nickname, description) in enumerate(OUTPUT_METADATA):
        if index >= component.Params.Output.Count:
            break
        parameter = component.Params.Output[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description


def empty_outputs(status, message, base_case="", fetch_log=None):
    return ([], [], "{}", base_case or "{}", [], [
        "MUSSELFLOW SITE DATA | build "+COMPONENT_BUILD+" | "+status,
        message,
    ], list(fetch_log or []))


def finite_number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def parse_site_request(value):
    """Parse the browser bridge without confusing a request with site data."""
    text = str(value or "").strip()
    if not text:
        return {}
    try:
        document = json.loads(text)
    except Exception as exception:
        raise ValueError("siteRequestJson is invalid JSON: %s" % exception)
    if not isinstance(document, dict):
        raise ValueError("siteRequestJson root must be a JSON object.")
    if document.get("schema") != SITE_REQUEST_SCHEMA:
        raise ValueError(
            "siteRequestJson schema must be %s." % SITE_REQUEST_SCHEMA)
    context = document.get("search_context")
    selection = document.get("selection")
    if not isinstance(context, dict):
        raise ValueError("siteRequestJson is missing search_context.")
    if not isinstance(selection, dict) or not selection.get("id"):
        raise ValueError("siteRequestJson is missing an accepted selection.")
    return document


def request_override(request, latitude, longitude, start_time, end_time, depths):
    """Resolve browser values first and retain manual inputs as fallbacks."""
    if not request:
        return latitude, longitude, start_time, end_time, depths
    context = request.get("search_context") or {}
    site = context.get("site") or {}
    period = context.get("time") or {}
    request_depths = context.get("depths_m") or []
    return (
        site.get("latitude", latitude),
        site.get("longitude", longitude),
        period.get("start") or start_time,
        period.get("end") or end_time,
        request_depths or depths,
    )


def select_product(latitude, longitude):
    """Return the regional product profile covering one WGS84 coordinate."""
    for product in PRODUCTS:
        min_lat, max_lat, min_lon, max_lon = product["bounds"]
        if (min_lat <= latitude <= max_lat and
                min_lon <= longitude <= max_lon):
            return product
    return None


def depth_list(values):
    candidates = []

    def visit(value, depth=0):
        if value is None or depth > 12:
            return
        if isinstance(value, (str, bytes, int, float)):
            candidates.append(value)
            return
        branches = getattr(value, "Branches", None)
        if branches is not None:
            for branch in branches:
                visit(branch, depth+1)
            return
        script_variable = getattr(value, "ScriptVariable", None)
        if callable(script_variable):
            try:
                unwrapped = script_variable()
                if unwrapped is not value:
                    visit(unwrapped, depth+1)
                    return
            except Exception:
                pass
        if hasattr(value, "Value"):
            try:
                unwrapped = value.Value
                if unwrapped is not value:
                    visit(unwrapped, depth+1)
                    return
            except Exception:
                pass
        try:
            for item in value:
                visit(item, depth+1)
        except TypeError:
            candidates.append(value)

    visit(values)
    result = []
    for value in candidates:
        number = finite_number(value)
        if number is None:
            continue
        number = abs(number)
        if number <= 1000.0 and number not in result:
            result.append(number)
    return sorted(result) if result else [1.0]


def parse_utc(value, label):
    text = "" if value is None else str(value).strip()
    if not text:
        raise ValueError(label+" is empty.")
    normalized = text[:-1]+"+00:00" if text.endswith("Z") else text
    try:
        result = datetime.datetime.fromisoformat(normalized)
    except ValueError:
        try:
            result = datetime.datetime.strptime(text, "%Y-%m-%d")
        except ValueError:
            raise ValueError(
                label+" must be ISO text such as 2026-04-20T00:00:00Z.")
    if result.tzinfo is None:
        result = result.replace(tzinfo=datetime.timezone.utc)
    return result.astimezone(datetime.timezone.utc)


def utc_text(value, daily=False):
    if daily:
        value = value.replace(hour=0, minute=0, second=0, microsecond=0)
    else:
        value = value.replace(minute=0, second=0, microsecond=0)
    return value.strftime("%Y-%m-%dT%H:%M:%SZ")


def frame_times(start, end, count):
    if end < start:
        raise ValueError("endTime precedes startTime.")
    if count == 1:
        return [start]
    span = (end-start).total_seconds()
    return [start+datetime.timedelta(seconds=span*i/float(count-1))
            for i in range(count)]


def wmts_position(latitude, longitude, level=TILE_LEVEL):
    """Return EPSG:4326 WMTS tile row/column and in-tile pixel."""
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


def ocean_colour_layer(layers, timestamp):
    """Choose OLCI MY for archives and NRT only for the recent rolling window."""
    archive = layers["satellite_chlorophyll"]
    recent = layers.get("satellite_chlorophyll_nrt")
    if not recent:
        return archive
    try:
        requested = datetime.datetime.strptime(
            str(timestamp)[:10], "%Y-%m-%d").date()
        age_days = (datetime.datetime.utcnow().date()-requested).days
    except Exception:
        return archive
    return recent if -2 <= age_days <= 35 else archive


def url_summary(url):
    """Readable WMTS request identity without hiding the exact URL."""
    try:
        query = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
        one = lambda key, default="-": query.get(key, [default])[0]
        return (
            "layer=%s | time=%s | depth=%s | tile=%s/%s/%s | pixel=%s,%s"
            % (one("layer"), one("time"), one("elevation", "surface"),
               one("tilematrix"), one("tilerow"), one("tilecol"),
               one("i"), one("j")))
    except Exception:
        return "unparsed WMTS URL"


def http_json(url, trace=None, request_id=None, phase="data"):
    request_label = "id=%s | phase=%s" % (request_id, phase)
    if trace is not None:
        trace.add("HTTP START", request_label+" | "+url_summary(url))
        trace.add("HTTP URL", "id=%s | %s" % (request_id, url))
    started = time.perf_counter()
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "MusselFlow-Rhino/"+COMPONENT_BUILD})
    try:
        with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_S) as response:
            raw_payload = response.read()
            status = getattr(response, "status", "unknown")
        payload = raw_payload.decode("utf-8")
        document = json.loads(payload)
    except Exception as exception:
        if trace is not None:
            trace.add(
                "HTTP ERROR",
                "%s | %.3fs | %s: %s" % (
                    request_label, time.perf_counter()-started,
                    type(exception).__name__, str(exception)))
        raise
    if trace is not None:
        trace.add(
            "HTTP OK",
            "%s | status=%s | bytes=%d | %.3fs" % (
                request_label, status, len(raw_payload),
                time.perf_counter()-started))
    return document


def feature_properties(document):
    try:
        features = document.get("features") or []
        properties = features[0].get("properties") or {}
        return properties
    except (AttributeError, IndexError):
        return {}


def scalar_value(properties):
    return finite_number(properties.get("value"))


def vector_value(properties):
    east = finite_number(properties.get("component1Value"))
    north = finite_number(properties.get("component2Value"))
    if east is None or north is None:
        return None
    return east, north


def fetch_urls(urls, trace=None, phase="data"):
    """Fetch each unique fixed-endpoint URL concurrently and retain errors."""
    unique = sorted(set(urls))
    documents = {}
    errors = {}
    if not unique:
        if trace is not None:
            trace.add("HTTP BATCH SKIP", "phase=%s | no URLs" % phase)
        return documents, errors
    workers = min(MAX_WORKERS, len(unique))
    if trace is not None:
        trace.add(
            "HTTP BATCH START",
            "phase=%s | submitted=%d | unique=%d | workers=%d"
            % (phase, len(urls), len(unique), workers))
    batch_started = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        future_urls = {}
        for url in unique:
            request_id = trace.request_id() if trace is not None else None
            future = pool.submit(http_json, url, trace, request_id, phase)
            future_urls[future] = url
        for future in concurrent.futures.as_completed(future_urls):
            url = future_urls[future]
            try:
                documents[url] = future.result()
            except Exception as exception:
                errors[url] = "%s: %s" % (
                    type(exception).__name__, str(exception))
    if trace is not None:
        trace.add(
            "HTTP BATCH END",
            "phase=%s | successful=%d | failed=%d | %.3fs"
            % (phase, len(documents), len(errors),
               time.perf_counter()-batch_started))
    return documents, errors


def haversine_km(lat_a, lon_a, lat_b, lon_b):
    radius = 6371.0088
    p1, p2 = math.radians(lat_a), math.radians(lat_b)
    dp = math.radians(lat_b-lat_a)
    dl = math.radians(lon_b-lon_a)
    term = (math.sin(dp/2.0)**2+
            math.cos(p1)*math.cos(p2)*math.sin(dl/2.0)**2)
    return 2.0*radius*math.asin(min(1.0, math.sqrt(term)))


def ring_points(latitude, longitude, radius_km):
    lat_scale = 1.0/111.32
    lon_scale = 1.0/max(1e-6, 111.32*math.cos(math.radians(latitude)))
    result = []
    for angle in range(0, 360, 45):
        radians = math.radians(angle)
        result.append((
            latitude+radius_km*math.sin(radians)*lat_scale,
            longitude+radius_km*math.cos(radians)*lon_scale))
    return result


def locate_wet_cell(latitude, longitude, timestamp, depth, layers, trace=None):
    """Find the requested or nearest sampled cell containing current data."""
    candidates_by_ring = [[(latitude, longitude)]]
    candidates_by_ring.extend(
        ring_points(latitude, longitude, radius) for radius in (1.0, 2.0, 4.0, 8.0))
    request_errors = []
    if trace is not None:
        trace.add(
            "WET CELL SEARCH START",
            "requested=%.6f,%.6f | time=%s | depth=%.3g m"
            % (latitude, longitude, timestamp, depth))
    for ring_index, candidates in enumerate(candidates_by_ring):
        radius = (0.0, 1.0, 2.0, 4.0, 8.0)[ring_index]
        if trace is not None:
            trace.add(
                "WET CELL RING",
                "radius=%.3g km | candidates=%d" % (radius, len(candidates)))
        urls = [feature_url(
            layers["current"], lat, lon, timestamp, depth)
            for lat, lon in candidates]
        documents, errors = fetch_urls(
            urls, trace, "wet_cell_%.3gkm" % radius)
        request_errors.extend(errors.values())
        valid = []
        for (lat, lon), url in zip(candidates, urls):
            properties = feature_properties(documents.get(url, {}))
            if vector_value(properties) is None:
                continue
            sampled_lat = finite_number(properties.get("lat"))
            sampled_lon = finite_number(properties.get("lon"))
            sampled_lat = lat if sampled_lat is None else sampled_lat
            sampled_lon = lon if sampled_lon is None else sampled_lon
            valid.append((
                haversine_km(latitude, longitude, sampled_lat, sampled_lon),
                sampled_lat, sampled_lon, properties))
        if valid:
            valid.sort(key=lambda item: item[0])
            if trace is not None:
                trace.add(
                    "WET CELL FOUND",
                    "sampled=%.6f,%.6f | offset=%.3f km | radius=%.3g km"
                    % (valid[0][1], valid[0][2], valid[0][0], radius))
            return valid[0], request_errors
        if trace is not None:
            trace.add(
                "WET CELL RING EMPTY",
                "radius=%.3g km | HTTP errors=%d"
                % (radius, len(errors)))
    if trace is not None:
        trace.add("WET CELL SEARCH FAILED", "no valid current cell within 8 km")
    return None, request_errors


def mean_or_none(values):
    accepted = [float(value) for value in values if value is not None]
    return sum(accepted)/len(accepted) if accepted else None


def round_numbers(value, digits=6):
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, float):
        return round(value, digits)
    if isinstance(value, list):
        return [round_numbers(item, digits) for item in value]
    if isinstance(value, dict):
        return {key: round_numbers(item, digits)
                for key, item in value.items()}
    return value


def canonical_json(value):
    return json.dumps(
        round_numbers(value), sort_keys=True, separators=(",", ":"),
        ensure_ascii=False, allow_nan=False)


def north_axes(vector):
    try:
        north = Rhino.Geometry.Vector3d(vector.X, vector.Y, 0.0)
        valid = north.IsValid and north.Length > 1e-9
    except Exception:
        valid = False
    if not valid:
        north = Rhino.Geometry.Vector3d.YAxis
    north.Unitize()
    east = Rhino.Geometry.Vector3d(north.Y, -north.X, 0.0)
    return east, north


def rhino_vector(east_speed, north_speed, north_vector):
    east_axis, north_axis = north_axes(north_vector)
    return east_axis*east_speed+north_axis*north_speed


def profile_record(documents, urls, depth):
    """Parse one depth profile and convert units used by the grammar."""
    properties = {
        name: feature_properties(documents.get(url, {}))
        for name, url in urls.items()}
    current = vector_value(properties["current"])
    result = {
        "depth_m": depth,
        "eastward_current_m_s": None if current is None else current[0],
        "northward_current_m_s": None if current is None else current[1],
        "temperature_c": scalar_value(properties["temperature"]),
        "salinity_psu": scalar_value(properties["salinity"]),
        # Copernicus regional BGC oxygen is mmol/m3; x 0.032 = mg/L.
        "dissolved_oxygen_mg_l": None,
        # mg/m3 is numerically equal to micrograms/L.
        "chlorophyll_a_ug_l": scalar_value(properties["chlorophyll"]),
        "phytoplankton_carbon_mmol_m3": scalar_value(
            properties["phytoplankton_carbon"]),
        "nitrate_mmol_m3": scalar_value(properties["nitrate"]),
        "phosphate_mmol_m3": scalar_value(properties["phosphate"]),
    }
    oxygen = scalar_value(properties["oxygen"])
    if oxygen is not None:
        result["dissolved_oxygen_mg_l"] = oxygen*0.032
    return result


def aggregate_profiles(
        profiles, tsm=None, turbidity=None, satellite_chlorophyll=None):
    east = mean_or_none([p["eastward_current_m_s"] for p in profiles])
    north = mean_or_none([p["northward_current_m_s"] for p in profiles])
    return {
        "eastward_current_m_s": east,
        "northward_current_m_s": north,
        "current_speed_m_s": (
            None if east is None or north is None
            else math.hypot(east, north)),
        "temperature_c": mean_or_none([p["temperature_c"] for p in profiles]),
        "salinity_psu": mean_or_none([p["salinity_psu"] for p in profiles]),
        "dissolved_oxygen_mg_l": mean_or_none([
            p["dissolved_oxygen_mg_l"] for p in profiles]),
        "chlorophyll_a_ug_l": mean_or_none([
            p["chlorophyll_a_ug_l"] for p in profiles]),
        # Copernicus HR SPM uses g/m3, numerically equal to mg/L.
        "tsm_mg_l": tsm,
        "turbidity_fnu": turbidity,
        "satellite_chlorophyll_a_ug_l": satellite_chlorophyll,
        "phytoplankton_carbon_mmol_m3": mean_or_none([
            p["phytoplankton_carbon_mmol_m3"] for p in profiles]),
        "nitrate_mmol_m3": mean_or_none([
            p["nitrate_mmol_m3"] for p in profiles]),
        "phosphate_mmol_m3": mean_or_none([
            p["phosphate_mmol_m3"] for p in profiles]),
    }


def existing_boundary(case_document):
    try:
        steps = case_document["forcing"]["steps"]
        return copy.deepcopy(steps[0]["boundary"])
    except (KeyError, IndexError, TypeError):
        return {
            "temperature_c": 12.0,
            "salinity_psu": 20.0,
            "dissolved_oxygen_mg_l": 9.0,
            "chlorophyll_a_ug_l": 5.0,
            "tsm_mg_l": 3.0,
        }


def patch_case(base_text, frames, duration_h):
    if not base_text or not str(base_text).strip():
        return None, "BaseModelJson is empty; SiteDataJson was fetched but no SimulationCaseJson was compiled."
    try:
        case = json.loads(str(base_text))
    except Exception as exception:
        return None, "BaseModelJson is invalid JSON: %s" % exception
    if not isinstance(case, dict):
        return None, "BaseModelJson root must be a JSON object."
    fallback = existing_boundary(case)
    # The current fitness timeline evaluates several flow vectors but supports
    # one shared environmental boundary. Preserve every frame-specific value
    # in SiteData, and use period means here so CaseJson remains executable.
    # A future per-step evaluator can consume the stored profiles directly.
    mapping = (
        ("temperature_c", "temperature_c"),
        ("salinity_psu", "salinity_psu"),
        ("dissolved_oxygen_mg_l", "dissolved_oxygen_mg_l"),
        ("chlorophyll_a_ug_l", "chlorophyll_a_ug_l"),
        ("tsm_mg_l", "tsm_mg_l"),
    )
    shared_boundary = copy.deepcopy(fallback)
    for target, source in mapping:
        period_mean = mean_or_none([
            frame["aggregate"].get(source) for frame in frames])
        if period_mean is not None:
            shared_boundary[target] = period_mean

    steps = []
    for index, frame in enumerate(frames):
        steps.append({
            "id": "site_%02d_%s" % (
                index, frame["time_utc"].replace("-", "").replace(":", "")),
            "order": index,
            "flow_vector_index": index,
            "duration_h": duration_h,
            "boundary": copy.deepcopy(shared_boundary),
        })
    forcing = case.setdefault("forcing", {})
    forcing["mode"] = "timeline"
    # Keep the grammar's strict public source label. The physical values still
    # come directly from this component's FlowVectors output, connected to the
    # optimizer input with the same name.
    forcing["vector_source"] = "Grasshopper.flowVectors"
    forcing["steps"] = steps
    forcing["repeat_count"] = 1
    initial_oxygen = shared_boundary.get("dissolved_oxygen_mg_l")
    if initial_oxygen is not None:
        case.setdefault("site", {})[
            "initial_dissolved_oxygen_mg_l"] = initial_oxygen
    return case, None


def format_value(frame):
    aggregate = frame["aggregate"]
    def show(key, unit):
        value = aggregate.get(key)
        return "missing" if value is None else "%.6g %s" % (value, unit)
    return (
        "%s | speed %s | T %s | salinity %s | DO %s | chl-a %s | TSM %s"
        % (frame["time_utc"], show("current_speed_m_s", "m/s"),
           show("temperature_c", "degC"), show("salinity_psu", "PSU"),
           show("dissolved_oxygen_mg_l", "mg/L"),
           show("chlorophyll_a_ug_l", "ug/L"), show("tsm_mg_l", "mg/L")))


def request_key(latitude, longitude, start, end, depths, count, base_text,
                north_vector, product, site_request=None):
    try:
        north = (north_vector.X, north_vector.Y, north_vector.Z)
    except Exception:
        north = (0.0, 1.0, 0.0)
    payload = canonical_json({
        "latitude": latitude,
        "longitude": longitude,
        "start": start.isoformat(),
        "end": end.isoformat(),
        "depths": depths,
        "count": count,
        "base_hash": hashlib.sha256(
            (base_text or "").encode("utf-8")).hexdigest(),
        "north": north,
        "product": product["id"],
        "catalogue_selection": (
            (site_request or {}).get("selection") or {}).get("id", ""),
    })
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def fetch_site(latitude, longitude, start, end, depths, count,
               base_case_text, north_vector, product, trace=None,
               site_request=None):
    started = time.perf_counter()
    trace = trace or FetchTrace()
    layers = product["layers"]
    trace.add(
        "FETCH SITE START",
        "region=%s | coordinate=%.6f,%.6f | frames=%d | depths=%s"
        % (product["id"], latitude, longitude, count,
           ",".join("%.3g" % depth for depth in depths)))
    trace.add(
        "LAYER ROUTE",
        " | ".join("%s=%s" % (name, layer)
                   for name, layer in sorted(layers.items())))
    requested_times = frame_times(start, end, count)
    # The physics product is hourly. Closely spaced requested frames can round
    # to the same model timestamp; collapse them so case IDs and vector indices
    # remain strictly one-to-one and unique.
    time_pairs = []
    seen_physics_times = set()
    for value in requested_times:
        physics_time = utc_text(value)
        if physics_time in seen_physics_times:
            continue
        seen_physics_times.add(physics_time)
        time_pairs.append((physics_time, utc_text(value, daily=True)))
    physics_times = [pair[0] for pair in time_pairs]
    daily_times = [pair[1] for pair in time_pairs]
    trace.add(
        "TIME PLAN",
        "requested=%d | unique hourly=%d | start=%s | end=%s"
        % (len(requested_times), len(physics_times),
           utc_text(start), utc_text(end)))

    wet, location_errors = locate_wet_cell(
        latitude, longitude, physics_times[0], depths[0], layers, trace)
    if wet is None:
        message = (
            "No wet current cell was found within 8 km. Move the coordinate "
            "offshore or choose a time covered by the selected regional product.")
        if location_errors:
            message += " First HTTP error: "+location_errors[0]
        trace.add("FETCH SITE STOP", "NO_WET_CELL")
        return empty_outputs(
            "NO_WET_CELL", message, base_case_text, trace.lines())

    distance_km, sampled_latitude, sampled_longitude, _ = wet
    requests = []
    frame_requests = []
    for frame_index, (physics_time, daily_time) in enumerate(
            zip(physics_times, daily_times)):
        profiles = []
        for depth in depths:
            urls = {
                "current": feature_url(
                    layers["current"], sampled_latitude, sampled_longitude,
                    physics_time, depth),
                "temperature": feature_url(
                    layers["temperature"], sampled_latitude, sampled_longitude,
                    physics_time, depth),
                "salinity": feature_url(
                    layers["salinity"], sampled_latitude, sampled_longitude,
                    physics_time, depth),
                "oxygen": feature_url(
                    layers["oxygen"], sampled_latitude, sampled_longitude,
                    daily_time, depth),
                "chlorophyll": feature_url(
                    layers["chlorophyll"], sampled_latitude, sampled_longitude,
                    daily_time, depth),
                "phytoplankton_carbon": feature_url(
                    layers["phytoplankton_carbon"], sampled_latitude,
                    sampled_longitude, daily_time, depth),
                "nitrate": feature_url(
                    layers["nitrate"], sampled_latitude, sampled_longitude,
                    daily_time, depth),
                "phosphate": feature_url(
                    layers["phosphate"], sampled_latitude, sampled_longitude,
                    daily_time, depth),
            }
            requests.extend(urls.values())
            profiles.append((depth, urls))
        surface_urls = {
            "tsm": feature_url(
                layers["tsm"], sampled_latitude, sampled_longitude, daily_time),
            "turbidity": feature_url(
                layers["turbidity"], sampled_latitude, sampled_longitude,
                daily_time),
            "satellite_chlorophyll": feature_url(
                ocean_colour_layer(layers, daily_time), sampled_latitude,
                sampled_longitude, daily_time),
        }
        requests.extend(surface_urls.values())
        frame_requests.append((frame_index, physics_time, profiles, surface_urls))

    trace.add(
        "DATA REQUEST PLAN",
        "frames=%d | depths=%d | submitted URLs=%d | unique URLs=%d"
        % (len(frame_requests), len(depths), len(requests), len(set(requests))))
    documents, request_errors = fetch_urls(requests, trace, "site_data")
    trace.add(
        "FRAME ASSEMBLY START",
        "documents=%d | HTTP errors=%d" % (len(documents), len(request_errors)))
    frames = []
    missing_counts = {}
    for _, timestamp, profile_specs, surface_urls in frame_requests:
        profiles = [profile_record(documents, urls, depth)
                    for depth, urls in profile_specs]
        tsm = scalar_value(feature_properties(
            documents.get(surface_urls["tsm"], {})))
        turbidity = scalar_value(feature_properties(
            documents.get(surface_urls["turbidity"], {})))
        satellite_chlorophyll = scalar_value(feature_properties(
            documents.get(surface_urls["satellite_chlorophyll"], {})))
        aggregate = aggregate_profiles(
            profiles, tsm, turbidity, satellite_chlorophyll)
        for key, value in aggregate.items():
            if value is None:
                missing_counts[key] = missing_counts.get(key, 0)+1
        if (aggregate["eastward_current_m_s"] is None or
                aggregate["northward_current_m_s"] is None):
            continue
        frames.append({
            "time_utc": timestamp,
            "depths_m": list(depths),
            "profiles": profiles,
            "aggregate": aggregate,
        })
        trace.add(
            "FRAME ASSEMBLED",
            "time=%s | profiles=%d | speed=%s m/s"
            % (timestamp, len(profiles),
               "missing" if aggregate["current_speed_m_s"] is None
               else "%.6f" % aggregate["current_speed_m_s"]))

    if not frames:
        detail = ""
        if request_errors:
            detail = " First HTTP error: "+next(iter(request_errors.values()))
        trace.add("FETCH SITE STOP", "NO_CURRENT_DATA")
        return empty_outputs(
            "NO_CURRENT_DATA",
            "The sampled wet cell returned no valid currents for the requested times."
            +detail,
            base_case_text,
            trace.lines())

    total_hours = max(1.0, (end-start).total_seconds()/3600.0)
    duration_h = total_hours/len(frames)
    patched_case, case_error = patch_case(
        base_case_text, frames, duration_h)
    trace.add(
        "CASE COMPILE",
        "frames=%d | duration_per_step=%.6f h | status=%s"
        % (len(frames), duration_h, "warning" if case_error else "ok"))
    case_text = "{}" if patched_case is None else canonical_json(patched_case)
    vectors = [rhino_vector(
        frame["aggregate"]["eastward_current_m_s"],
        frame["aggregate"]["northward_current_m_s"], north_vector)
        for frame in frames]
    output_times = [frame["time_utc"] for frame in frames]

    site_document = {
        "schema": "musselflow.site_data.1.0.0",
        "build": COMPONENT_BUILD,
        "data_state": "sampled",
        "source_type": {
            "physics": "regional numerical analysis/forecast model",
            "biogeochemistry": "regional numerical analysis/forecast model",
            "tsm_turbidity": "satellite-derived surface observation",
        },
        "regional_product": {
            "id": product["id"],
            "name": product["name"],
        },
        "requested": {
            "latitude": latitude,
            "longitude": longitude,
            "start_utc": utc_text(start),
            "end_utc": utc_text(end),
            "depths_m": list(depths),
            "frame_count": count,
        },
        "resolved_frame_count": len(frames),
        "sampled": {
            "latitude": sampled_latitude,
            "longitude": sampled_longitude,
            "distance_from_request_km": distance_km,
        },
        "layers": copy.deepcopy(layers),
        "frames": frames,
        "missing_by_aggregate_field": missing_counts,
        "http_error_count": len(request_errors),
        "scientific_status": "UNVALIDATED_BOUNDARY_DATA",
        "limitations": [
            "Ambient grid data do not resolve farm-scale wakes or socks.",
            "Depth profiles are equally averaged; biomass weighting is not yet active.",
            "SimulationCaseJson uses period-mean chemistry because the current fitness timeline requires one shared boundary; frame-specific chemistry remains in SiteDataJson.",
            "Satellite TSM/turbidity are surface products and may be missing.",
            "WMTS values are point samples, not local field validation.",
        ],
    }
    if site_request:
        site_document["catalogue_request"] = {
            "selection": copy.deepcopy(site_request.get("selection") or {}),
            "search_context": copy.deepcopy(
                site_request.get("search_context") or {}),
        }
    report = [
        "MUSSELFLOW SITE DATA | build %s | %d valid frames | %.3fs"
        % (COMPONENT_BUILD, len(frames), time.perf_counter()-started),
        "SITE | requested %.6f, %.6f | sampled %.6f, %.6f | offset %.3f km"
        % (latitude, longitude, sampled_latitude, sampled_longitude, distance_km),
        "REGION | %s | automatic product routing" % product["name"],
        "DEPTHS | %s m | equally averaged | nearest model layers"
        % ", ".join("%.3g" % depth for depth in depths),
        "SOURCE | Copernicus regional physics hourly + BGC daily + ocean colour",
        "CASE | flow varies by frame; chemistry uses period means in the current fitness bridge",
        "CONTRACT | %d FlowVectors | %d SimulationCaseJson timeline steps"
        % (len(vectors), len(frames)),
        "CACHE | in-memory request cache; no network call during Galapagos",
    ]
    if site_request:
        report.insert(
            3, "CATALOGUE SELECTION | %s" % (
                (site_request.get("selection") or {}).get("id") or "unknown"))
    if distance_km > 0.5:
        report.append(
            "LOCATION WARNING | requested point was dry/masked; nearest wet cell used.")
    if len(physics_times) < count:
        report.append(
            "TIME NOTE | %d requested frames collapsed to %d unique hourly states."
            % (count, len(physics_times)))
    if missing_counts:
        report.append(
            "MISSING | "+", ".join(
                "%s %d/%d" % (key, value, len(frames))
                for key, value in sorted(missing_counts.items())))
    if request_errors:
        report.append(
            "HTTP WARNING | %d requests failed; successful fields were retained."
            % len(request_errors))
    if case_error:
        report.append("CASE WARNING | "+case_error)
    report.extend([
        "MODEL LIMIT | ambient regional data, not farm-scale CFD or measurements.",
        "VALIDATION LIMIT | inspect product quality and local monitoring before use.",
    ])
    trace.add(
        "FETCH SITE COMPLETE",
        "valid frames=%d | missing fields=%d | HTTP errors=%d | %.3fs"
        % (len(frames), len(missing_counts), len(request_errors),
           time.perf_counter()-started))
    return (
        vectors,
        output_times,
        canonical_json(site_document),
        case_text,
        [format_value(frame) for frame in frames],
        report,
        trace.lines(),
    )


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def BeforeRunScript(self):
        apply_component_metadata(getattr(self, "Component", None))

    def RunScript(
            self,
            fetch: bool,
            latitude: float,
            longitude: float,
            startTime: str,
            endTime: str,
            depths: list[float],
            frameCount: int,
            baseCaseJson: str,
            northVector: Rhino.Geometry.Vector3d,
            siteRequestJson: str):
        """
        Inputs:
            fetch: Button to fetch/cache this site request {item,bool}
            latitude: WGS84 latitude decimal degrees {item,float}
            longitude: WGS84 longitude decimal degrees {item,float}
            startTime: First UTC ISO timestamp {item,str}
            endTime: Last UTC ISO timestamp; empty means +24h {item,str}
            depths: Positive depth samples in metres {list,float}
            frameCount: Number of chronological states, 1-10 {item,int}
            BaseModelJson: Existing complete MusselFlow grammar JSON {item,str}
            northVector: Optional Rhino-plan geographic north {item,Vector3d}
            FetchRequestJson: Browser request overriding site/time/depths {item,str}

        Outputs, in this exact order:
            FlowVectors: Ambient current vectors in m/s {list,Vector3d}
            Times: UTC timestamps matching FlowVectors {list,str}
            SiteDataJson: Canonical environmental/provenance JSON {item,str}
            SimulationCaseJson: Site-patched model JSON {item,str}
            Values: Human-readable frame values {list,str}
            Report: Coverage, warnings, and limitations {list,str}
            FetchLog: Chronological cache, HTTP, and assembly trace {list,str}
        """
        trace = FetchTrace()
        trace.add("RUN START", "fetch=%s" % bool(fetch))
        try:
            site_request = parse_site_request(siteRequestJson)
            latitude, longitude, startTime, endTime, depths = request_override(
                site_request, latitude, longitude, startTime, endTime, depths)
            if site_request:
                trace.add(
                    "SITE REQUEST APPLIED",
                    "selection=%s | browser context overrides site/time/depths"
                    % ((site_request.get("selection") or {}).get("id") or "unknown"))
        except Exception as exception:
            trace.add("SITE REQUEST ERROR", str(exception))
            return empty_outputs(
                "INVALID_SITE_REQUEST", str(exception), baseCaseJson,
                trace.lines())
        lat = finite_number(latitude)
        lon = finite_number(longitude)
        if lat is None or lon is None:
            trace.add("VALIDATION ERROR", "latitude or longitude is not finite")
            return empty_outputs(
                "INVALID_COORDINATE", "Connect finite latitude and longitude.",
                baseCaseJson, trace.lines())
        if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
            trace.add("VALIDATION ERROR", "coordinate is outside WGS84 bounds")
            return empty_outputs(
                "INVALID_COORDINATE",
                "Latitude must be -90..90 and longitude -180..180.",
                baseCaseJson, trace.lines())
        product = select_product(lat, lon)
        if product is None:
            trace.add("ROUTING ERROR", "coordinate is outside configured products")
            return empty_outputs(
                "UNSUPPORTED_REGION",
                "No regional product profile covers this coordinate yet. "
                "Current profiles: Baltic Sea and Northwest European Shelf.",
                baseCaseJson, trace.lines())
        trace.add(
            "REGION SELECTED",
            "%s (%s)" % (product["name"], product["id"]))
        try:
            start = parse_utc(startTime, "startTime")
            if endTime is None or not str(endTime).strip():
                end = start+datetime.timedelta(hours=24)
            else:
                end = parse_utc(endTime, "endTime")
            count = 6 if frameCount is None else int(frameCount)
            count = max(1, min(10, count))
            requested_depths = depth_list(depths)
            if end < start:
                raise ValueError("endTime precedes startTime.")
        except Exception as exception:
            trace.add("VALIDATION ERROR", str(exception))
            return empty_outputs(
                "INVALID_REQUEST", str(exception), baseCaseJson, trace.lines())

        trace.add(
            "REQUEST NORMALIZED",
            "coordinate=%.6f,%.6f | start=%s | end=%s | frames=%d | depths=%s"
            % (lat, lon, utc_text(start), utc_text(end), count,
               ",".join("%.3g" % depth for depth in requested_depths)))

        key = request_key(
            lat, lon, start, end, requested_depths, count,
            baseCaseJson or "", northVector, product, site_request)
        trace.add("CACHE LOOKUP", "key=%s" % key[:16])
        if key in _MEMORY_CACHE:
            cached = _MEMORY_CACHE[key]
            report = list(cached[5])
            report.insert(1, "CACHE HIT | no Copernicus request was made.")
            trace.add("CACHE HIT", "no HTTP request made | key=%s" % key[:16])
            trace.add("RUN COMPLETE", "served from memory cache")
            return (
                cached[0], cached[1], cached[2], cached[3], cached[4], report,
                trace.lines())
        trace.add("CACHE MISS", "network fetch required | key=%s" % key[:16])
        if not fetch:
            trace.add("RUN WAITING", "fetch is False; no HTTP request made")
            return empty_outputs(
                "WAITING", "Press the fetch Button once for this request.",
                baseCaseJson, trace.lines())

        try:
            result = fetch_site(
                lat, lon, start, end, requested_depths, count,
                baseCaseJson or "", northVector, product, trace, site_request)
        except Exception as exception:
            trace.add(
                "FETCH EXCEPTION",
                "%s: %s" % (type(exception).__name__, str(exception)))
            return empty_outputs(
                "FETCH_ERROR", "%s: %s" % (
                    type(exception).__name__, str(exception)), baseCaseJson,
                trace.lines())
        if result[0]:
            trace.add("CACHE STORE", "key=%s | successful result stored" % key[:16])
            trace.add("RUN COMPLETE", "network result returned and cached")
            result = result[:6]+(trace.lines(),)
            _MEMORY_CACHE[key] = result
        else:
            trace.add("RUN COMPLETE", "fetch ended without valid FlowVectors")
            result = result[:6]+(trace.lines(),)
        return result

    def AfterRunScript(self):
        pass
