"""
Copernicus Search Context
=========================

Small Rhino 8 / Grasshopper Python 3 SDK-mode component that packages a site
and sampling intention into one canonical JSON wire. Connect
``SearchContextJson`` to the Copernicus Data Browser input with the same name.
The browser uses it to narrow and
rank catalogue products; it does not claim that measurements have been fetched.

Inputs
------
latitude : item / float
    WGS84 latitude in decimal degrees.
longitude : item / float
    WGS84 longitude in decimal degrees.
startTime : item / str
    Optional first UTC timestamp, e.g. 2026-04-20T00:00:00Z.
endTime : item / str
    Optional final UTC timestamp.
depths : list / float
    Positive sampling depths in metres.
variables : list / str
    Desired fields, e.g. current, temperature, salinity, oxygen, chlorophyll,
    phytoplankton_carbon, nitrate, phosphate, TSM or turbidity.
modelLayer : item / str
    Optional preferred product, dataset or complete WMTS layer identifier.
pixelSize : item / float
    Optional desired spatial pixel/grid size in metres. This is a catalogue
    preference; the source product still determines its actual resolution.
keywords : list / str
    Optional extra catalogue terms.

Outputs -- create three ports once in this exact order
------------------------------------------------------
SearchContextJson : item / str
    Canonical search-context JSON; connect to Copernicus Data
    Browser.searchContext.
Query : item / str
    Primary live-catalogue query chosen from the context.
Report : list / str
    Region inference, normalized parameters and wiring guidance.
"""

import datetime
import json
import math
import re

import Grasshopper


COMPONENT_BUILD = "2026-08-08c"
SCHEMA = "copernicus.search_context.1.0"

COMPONENT_METADATA = (
    "Copernicus Search Context",
    "CopernicusContext",
    "Package site, time, depth, variables and layer preferences into one "
    "catalogue-search context for the Copernicus Data Browser.",
)

INPUT_METADATA = (
    ("latitude", "lat", "WGS84 latitude in decimal degrees."),
    ("longitude", "lon", "WGS84 longitude in decimal degrees."),
    ("startTime", "start", "Optional first UTC ISO timestamp."),
    ("endTime", "end", "Optional final UTC ISO timestamp."),
    ("depths", "depths", "Positive sampling-depth list in metres."),
    ("variables", "variables", "Desired fields: current, oxygen, chlorophyll, nutrients, etc."),
    ("modelLayer", "layer", "Optional product, dataset or complete WMTS layer identifier."),
    ("pixelSize", "pixel", "Optional desired spatial pixel/grid size in metres."),
    ("keywords", "keywords", "Optional additional catalogue search terms."),
)

OUTPUT_METADATA = (
    ("SearchContextJson", "SearchContext", "Catalogue intent only; connect to Data Browser.SearchContextJson."),
    ("Query", "Query", "Primary query used when searching the live catalogue."),
    ("Report", "Report", "Normalized context, inferred region and wiring guidance."),
)

REGIONS = (
    ("Baltic Sea", (53.0, 66.0, 9.0, 31.0), ["Baltic"]),
    ("Northwest European Shelf", (46.0, 63.0, -17.0, 13.0),
     ["Northwest Shelf", "North Sea"]),
)

VARIABLE_ALIASES = {
    "flow": "current",
    "flow vector": "current",
    "flow vectors": "current",
    "current speed": "current",
    "sea water velocity": "current",
    "velocity": "current",
    "temp": "temperature",
    "o2": "oxygen",
    "dissolved oxygen": "oxygen",
    "chl": "chlorophyll",
    "chlorophyll a": "chlorophyll",
    "phyc": "phytoplankton_carbon",
    "phytoplankton": "phytoplankton_carbon",
    "no3": "nitrate",
    "po4": "phosphate",
    "spm": "tsm",
    "suspended particulate matter": "tsm",
    "tur": "turbidity",
}


def apply_metadata(component):
    if component is None:
        return
    component.Name, component.NickName, component.Description = COMPONENT_METADATA
    component.Message = "Search context"
    for index, (name, nickname, description) in enumerate(INPUT_METADATA):
        if index >= component.Params.Input.Count:
            break
        parameter = component.Params.Input[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description
        parameter.Optional = index in (2, 3, 4, 5, 6, 7, 8)
        if index in (4, 5, 8):
            parameter.Access = Grasshopper.Kernel.GH_ParamAccess.list
    for index, (name, nickname, description) in enumerate(OUTPUT_METADATA):
        if index >= component.Params.Output.Count:
            break
        parameter = component.Params.Output[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description
        if index == 2:
            parameter.Access = Grasshopper.Kernel.GH_ParamAccess.list


def finite(value, name, required=False):
    if value is None or (isinstance(value, str) and not value.strip()):
        if required:
            raise ValueError(name+" is required.")
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(name+" must be a finite number.")
    if not math.isfinite(number):
        raise ValueError(name+" must be a finite number.")
    return number


def list_values(value):
    """Flatten GH lists, DataTrees, branches and Goo wrappers safely."""
    result = []

    def visit(item, depth=0):
        if item is None or depth > 12:
            return
        if isinstance(item, (str, bytes, int, float)):
            result.append(item)
            return
        branches = getattr(item, "Branches", None)
        if branches is not None:
            for branch in branches:
                visit(branch, depth+1)
            return
        script_variable = getattr(item, "ScriptVariable", None)
        if callable(script_variable):
            try:
                unwrapped = script_variable()
                if unwrapped is not item:
                    visit(unwrapped, depth+1)
                    return
            except Exception:
                pass
        if hasattr(item, "Value"):
            try:
                unwrapped = item.Value
                if unwrapped is not item:
                    visit(unwrapped, depth+1)
                    return
            except Exception:
                pass
        try:
            iterator = iter(item)
        except TypeError:
            result.append(item)
            return
        for child in iterator:
            visit(child, depth+1)

    visit(value)
    return result


def normalized_texts(values, aliases=None):
    result = []
    for value in list_values(values):
        text = " ".join(str(value or "").strip().lower().split())
        if not text:
            continue
        text = (aliases or {}).get(text, text)
        if text not in result:
            result.append(text)
    return result


def normalized_depths(values):
    result = []
    for raw_value in list_values(values):
        candidates = [raw_value]
        if isinstance(raw_value, str):
            candidates = [value for value in re.split(
                r"[,;\s]+", raw_value.strip()) if value]
        for value in candidates:
            try:
                number = finite(value, "depth")
            except ValueError:
                continue
            if number is None:
                continue
            number = abs(number)
            if number <= 1000.0 and number not in result:
                result.append(number)
    return sorted(result)


def parse_time(value, name):
    text = str(value or "").strip()
    if not text:
        return "", None
    normalized = text[:-1]+"+00:00" if text.endswith("Z") else text
    try:
        moment = datetime.datetime.fromisoformat(normalized)
    except ValueError:
        raise ValueError(name+" must be ISO text, e.g. 2026-04-20T00:00:00Z.")
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=datetime.timezone.utc)
    moment = moment.astimezone(datetime.timezone.utc)
    return moment.strftime("%Y-%m-%dT%H:%M:%SZ"), moment


def infer_region(latitude, longitude):
    for name, bounds, terms in REGIONS:
        min_lat, max_lat, min_lon, max_lon = bounds
        if min_lat <= latitude <= max_lat and min_lon <= longitude <= max_lon:
            return {"name": name, "catalogue_terms": list(terms)}
    return {"name": "unconfigured/global", "catalogue_terms": []}


def product_hint(layer):
    text = str(layer or "").strip()
    if not text:
        return ""
    return text.split("/", 1)[0].strip()


def build_context(latitude, longitude, start_time, end_time, depths,
                  variables, model_layer, pixel_size, keywords):
    lat = finite(latitude, "latitude", True)
    lon = finite(longitude, "longitude", True)
    if not -90.0 <= lat <= 90.0 or not -180.0 <= lon <= 180.0:
        raise ValueError("latitude/longitude is outside WGS84 bounds.")
    start_text, start = parse_time(start_time, "startTime")
    end_text, end = parse_time(end_time, "endTime")
    if start is not None and end is not None and end < start:
        raise ValueError("endTime precedes startTime.")

    depth_values = normalized_depths(depths)
    variable_values = normalized_texts(variables, VARIABLE_ALIASES)
    keyword_values = normalized_texts(keywords)
    layer = str(model_layer or "").strip()
    pixel = finite(pixel_size, "pixelSize")
    if pixel is not None and pixel <= 0.0:
        raise ValueError("pixelSize must be greater than zero.")
    region = infer_region(lat, lon)
    product = product_hint(layer)

    primary_query = (
        variable_values[0] if variable_values else
        product if product else
        region["catalogue_terms"][0] if region["catalogue_terms"] else
        keyword_values[0] if keyword_values else "")
    context = {
        "schema": SCHEMA,
        "build": COMPONENT_BUILD,
        "role": "catalogue_filter_not_measurements",
        "site": {
            "latitude": lat,
            "longitude": lon,
            "region_hint": region,
        },
        "time": {"start": start_text, "end": end_text},
        "depths_m": depth_values,
        "variables": variable_values,
        "model_layer": layer,
        "product_hint": product,
        "desired_pixel_size_m": pixel,
        "keywords": keyword_values,
        "suggested_live_query": primary_query,
    }
    return context, primary_query


def canonical_json(value):
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def BeforeRunScript(self):
        apply_metadata(getattr(self, "Component", None))

    def RunScript(
            self,
            latitude: float,
            longitude: float,
            startTime: str,
            endTime: str,
            depths: list[float],
            variables: list[str],
            modelLayer: str,
            pixelSize: float,
            keywords: list[str]):
        """Build one reusable Copernicus catalogue-search context."""
        try:
            context, query = build_context(
                latitude, longitude, startTime, endTime, depths, variables,
                modelLayer, pixelSize, keywords)
        except Exception as exception:
            return "{}", "", [
                "COPERNICUS SEARCH CONTEXT | build %s | INVALID" % COMPONENT_BUILD,
                "%s: %s" % (type(exception).__name__, str(exception)),
            ]

        region = context["site"]["region_hint"]["name"]
        report = [
            "COPERNICUS SEARCH CONTEXT | build %s | READY" % COMPONENT_BUILD,
            "SITE | %.6f, %.6f | %s" % (
                context["site"]["latitude"], context["site"]["longitude"], region),
            "TIME | %s -> %s" % (
                context["time"]["start"] or "unspecified",
                context["time"]["end"] or "unspecified"),
            "DEPTHS | %s m" % (
                ", ".join("%g" % value for value in context["depths_m"])
                or "unspecified"),
            "VARIABLES | %s" % (", ".join(context["variables"]) or "unspecified"),
            "LAYER | %s" % (context["model_layer"] or "automatic/unspecified"),
            "PIXEL PREFERENCE | %s m" % (
                "%g" % context["desired_pixel_size_m"]
                if context["desired_pixel_size_m"] is not None else "unspecified"),
            "WIRING | SearchContextJson -> Copernicus Data Browser.SearchContextJson",
            "SCOPE | filters catalogue products only; Site Data performs the actual sampling.",
        ]
        return canonical_json(context), query, report

    def AfterRunScript(self):
        pass
