"""
Copernicus Data Browser
=======================

Clean Rhino 8 / Grasshopper Python 3 SDK-mode browser for Copernicus data.
It is the native Eto.Forms companion to ``Copernicus Universal Data Broker``.

The browser is useful before any web request: six verified Baltic and Northwest
Shelf products are bundled as a starter register. Text typed into the search
box filters locally. ``Search Copernicus Live`` searches the official CDSE STAC
and Copernicus Marine CSW catalogues on a background thread, leaving Rhino
interactive. A failed live request never destroys the local or cached results.

Inputs
------
openBrowser : item / bool
    Button or Toggle. True opens the modeless browser or brings it forward.
data : item / str
    Optional Data JSON from Copernicus Universal Data Broker.
collections : list / str
    Optional Collections output from the broker.
items : list / str
    Optional Items output from the broker.
assets : list / str
    Optional Assets output from the broker.
fetchLog : list / str
    Optional FetchLog output from the broker.
windowTitle : item / str
    Optional title. Empty uses ``Copernicus Data Browser``.
SearchContextJson : item / str
    SearchContextJson from Copernicus Search Context. This explicit final input was
    added without shifting the seven original browser inputs.

Outputs -- create six ports once in this exact order
---------------------------------------------------
SelectedType : item / str
    ``collection``, ``item`` or ``asset``.
SelectedId : item / str
    Accepted collection ID, item ID or asset key.
SelectedUrl : item / str
    Best available URL for the accepted selection.
FetchRequestJson : item / str
    Search context and accepted catalogue product merged into one explicit
    ``copernicus.site_request.1.0`` request. Connect this to MusselFlow Site
    Data ``FetchRequestJson``. The sampled values are returned by that component.
IsOpen : item / bool
    True while the modeless browser is visible.
Report : list / str
    Window, payload, search and wiring state.

SDK setup
---------
Paste this file into a Rhino 8 Python 3 component and choose
``Convert To GH_ScriptInstance``. The typed RunScript signature creates the
eight inputs. Add the six outputs once in the order above. BeforeRunScript
assigns names, list access and human-readable hover tooltips.

Architecture
------------
``CatalogueModel`` owns data and filtering. ``CatalogueClient`` owns bounded,
allow-listed HTTP access. ``CopernicusBrowserForm`` owns Eto presentation.
``BrowserSession`` owns one Grasshopper component's sticky lifetime and output
state. The layers communicate through narrow method calls, not shared edits.
"""

import concurrent.futures
import hashlib
import json
import threading
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import Eto.Drawing as drawing
import Eto.Forms as forms
import Grasshopper
import Rhino
import System
import scriptcontext as sc


# -----------------------------------------------------------------------------
# Contract and configuration
# -----------------------------------------------------------------------------

COMPONENT_BUILD = "2026-08-08g"
DEFAULT_TITLE = "Copernicus Data Browser"
SEARCH_CONTEXT_SCHEMA = "copernicus.search_context.1.0"
SITE_REQUEST_SCHEMA = "copernicus.site_request.1.0"

CDSE_URL = "https://stac.dataspace.copernicus.eu/v1/collections"
CDSE_HOST = "stac.dataspace.copernicus.eu"
MARINE_URL = (
    "https://csw.marine.copernicus.eu/geonetwork/"
    "csw-MYOCEAN-CORE-PRODUCTS/eng/csw")
MARINE_HOST = "csw.marine.copernicus.eu"
ALLOWED_HOSTS = frozenset((CDSE_HOST, MARINE_HOST))

HTTP_TIMEOUT_S = 30
MAX_RESPONSE_BYTES = 20*1024*1024
MAX_STAC_PAGES = 50
MAX_COLLECTIONS = 5000

COMPONENT_METADATA = (
    "Copernicus Data Browser",
    "CopernicusBrowser",
    "Browse, inspect and select Copernicus products in a native modeless Eto "
    "window. Includes an offline starter register and non-blocking live search.",
)

INPUT_METADATA = (
    ("openBrowser", "open", "Button: open or bring forward the Copernicus browser."),
    ("data", "data", "Optional Data JSON from Copernicus Universal Data Broker."),
    ("collections", "collections", "Optional Collections list from the data broker."),
    ("items", "items", "Optional Items list from the data broker."),
    ("assets", "assets", "Optional Assets list from the data broker."),
    ("fetchLog", "log", "Optional FetchLog list from the data broker."),
    ("windowTitle", "title", "Optional browser title; empty uses the default."),
    ("SearchContextJson", "SearchContext", "Catalogue intent from Copernicus Search Context; contains no sampled values."),
)

OUTPUT_METADATA = (
    ("SelectedType", "Type", "Accepted selection type: collection, item or asset."),
    ("SelectedId", "Id", "Accepted collection ID, item ID or asset key."),
    ("SelectedUrl", "Url", "Best available URL for the accepted selection."),
    ("FetchRequestJson", "FetchRequest", "Accepted product plus search context; connect to Site Data.FetchRequestJson."),
    ("IsOpen", "Open", "True while this component's browser is visible."),
    ("Report", "Report", "Browser state, search state and downstream wiring notes."),
)

STARTER_PRODUCTS = (
    (
        "BALTICSEA_ANALYSISFORECAST_PHY_003_006",
        "[Starter | Baltic] Physics forecast",
        "Currents, sea-water velocity, temperature and salinity.",
        "Baltic current currents velocity flow physics temperature salinity "
        "hydrodynamics forecast",
    ),
    (
        "BALTICSEA_ANALYSISFORECAST_BGC_003_007",
        "[Starter | Baltic] Biogeochemistry forecast",
        "Oxygen, chlorophyll, phytoplankton carbon, nitrate and phosphate.",
        "Baltic biogeochemistry oxygen chlorophyll phytoplankton nitrate "
        "phosphate nutrients eutrophication food",
    ),
    (
        "OCEANCOLOUR_BAL_BGC_HR_L4_NRT_009_208",
        "[Starter | Baltic] High-resolution ocean colour",
        "Satellite chlorophyll, suspended particulate matter and turbidity.",
        "Baltic ocean colour satellite chlorophyll phytoplankton turbidity "
        "TSM SPM water quality",
    ),
    (
        "NWSHELF_ANALYSISFORECAST_PHY_004_013",
        "[Starter | Northwest Shelf] Physics forecast",
        "Currents, sea-water velocity, temperature and salinity.",
        "Northwest Shelf North Sea current currents velocity flow physics "
        "temperature salinity hydrodynamics",
    ),
    (
        "NWSHELF_ANALYSISFORECAST_BGC_004_002",
        "[Starter | Northwest Shelf] Biogeochemistry forecast",
        "Oxygen, chlorophyll, phytoplankton carbon, nitrate and phosphate.",
        "Northwest Shelf North Sea biogeochemistry oxygen chlorophyll "
        "phytoplankton nitrate phosphate nutrients eutrophication food",
    ),
    (
        "OCEANCOLOUR_NWS_BGC_HR_L3_NRT_009_203",
        "[Starter | Northwest Shelf] High-resolution ocean colour",
        "Satellite chlorophyll, suspended particulate matter and turbidity.",
        "Northwest Shelf North Sea ocean colour satellite chlorophyll "
        "phytoplankton turbidity TSM SPM water quality",
    ),
)

SESSION_KEY = "copernicus.data.browser.sessions"
CACHE_KEY = "copernicus.data.browser.catalogue.cache"
SESSIONS = sc.sticky.setdefault(SESSION_KEY, {})
CATALOGUE_CACHE = sc.sticky.setdefault(CACHE_KEY, {})


# -----------------------------------------------------------------------------
# Pure data helpers
# -----------------------------------------------------------------------------

def text_list(value):
    """Convert Grasshopper item/list/tree-like values without splitting text."""
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        return [str(value)]
    try:
        return [str(item) for item in value if item is not None]
    except TypeError:
        return [str(value)]


def parse_json(source):
    if source is None or not str(source).strip():
        return {}
    try:
        value = json.loads(str(source))
        return value if isinstance(value, dict) else {"response": value}
    except Exception:
        return {"raw_text": str(source)}


def pretty(value):
    try:
        return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)
    except Exception:
        return str(value)


def payload_hash(data, collections, items, assets, fetch_log, context):
    digest = hashlib.sha256()
    for value in [str(data or ""), str(context or "")] + collections + items + assets + fetch_log:
        digest.update(str(value).encode("utf-8", "replace"))
        digest.update(b"\x00")
    return digest.hexdigest()


def response_object(document):
    response = document.get("response") if isinstance(document, dict) else None
    return response if isinstance(response, dict) else {}


def search_context(document):
    """Return a validated-enough context object without treating it as data."""
    if not isinstance(document, dict):
        return {}
    candidate = document.get("search_context")
    if not isinstance(candidate, dict):
        candidate = document
    return candidate if candidate.get("schema") == SEARCH_CONTEXT_SCHEMA else {}


def context_summary(context):
    if not context:
        return "none"
    site = context.get("site") or {}
    region = site.get("region_hint") or {}
    variables = context.get("variables") or []
    return "%s | %.5f, %.5f | %s" % (
        region.get("name") or "global",
        float(site.get("latitude", 0.0)), float(site.get("longitude", 0.0)),
        ", ".join(str(value) for value in variables) or "all variables")


def make_site_request(context, kind, identifier, url, detail):
    """Merge a catalogue selection with its spatial/temporal search context."""
    metadata = parse_json(detail)
    if not isinstance(metadata, dict):
        metadata = {"raw": str(detail or "")}
    return {
        "schema": SITE_REQUEST_SCHEMA,
        "build": COMPONENT_BUILD,
        "role": "site_data_fetch_request_not_measurements",
        "data_state": "not_fetched",
        "search_context": context if isinstance(context, dict) else {},
        "selection": {
            "type": str(kind or ""),
            "id": str(identifier or ""),
            "url": str(url or ""),
            "metadata": metadata,
        },
        "fetch_plan": {
            "mode": "mussel_site_bundle",
            "primary_product_id": str(identifier or ""),
            "requested_variables": list((context or {}).get("variables") or []),
            "complementary_sources": [
                "regional physics", "regional biogeochemistry", "ocean colour"],
            "note": (
                "MusselFlow Site Data adds the complementary regional products "
                "required for flow and ecological boundary conditions."),
        },
    }


def context_matches(context, value, row=""):
    """Apply region as a boundary and variables/keywords as inclusive families."""
    if not context:
        return True
    text = (str(row)+" "+pretty(value)).lower().replace("_", " ")
    site = context.get("site") or {}
    region = site.get("region_hint") or {}
    region_terms = [str(term).lower() for term in region.get("catalogue_terms") or []]
    variables = [str(term).lower().replace("_", " ")
                 for term in context.get("variables") or []]
    keywords = [str(term).lower() for term in context.get("keywords") or []]
    product = str(context.get("product_hint") or "").lower().replace("_", " ")

    if product and product not in text:
        return False
    if region_terms and not any(term in text for term in region_terms):
        return False
    if variables and not any(
            term in text or term.split(" ", 1)[0] in text for term in variables):
        return False
    if keywords and not any(term in text for term in keywords):
        return False
    return True


def collection_objects(document):
    """Read STAC collections, Marine products or project-neutral OGC layers."""
    values = response_object(document).get("collections") or []
    if not values and isinstance(document, dict):
        values = document.get("layers") or []
    return [value for value in values if isinstance(value, dict)]


def feature_objects(document):
    values = response_object(document).get("features") or []
    return [value for value in values if isinstance(value, dict)]


def starter_document():
    products = []
    for product_id, title, description, keywords in STARTER_PRODUCTS:
        products.append({
            "id": product_id,
            "title": title,
            "description": description,
            "keywords": keywords,
            "source": "MusselFlow verified starter register",
            "links": [{
                "rel": "self",
                "type": "text/html",
                "href": (
                    "https://data.marine.copernicus.eu/product/%s/description"
                    % urllib.parse.quote(product_id, safe="")),
            }],
            "register_status": (
                "Bundled starter metadata; refresh the live catalogue before "
                "production use."),
        })
    return {
        "mode": "starter_register",
        "build": COMPONENT_BUILD,
        "response": {"collections": products},
    }


def with_collections(base, *sources):
    """Attach deduplicated collections while preserving base features/metadata."""
    result = dict(base) if isinstance(base, dict) else {}
    response = dict(response_object(result))
    merged = []
    seen = set()
    for source in sources:
        for value in collection_objects(source):
            identifier = str(value.get("id") or "").strip()
            if not identifier or identifier in seen:
                continue
            seen.add(identifier)
            merged.append(value)
    response["collections"] = merged
    result["response"] = response
    return result


def collection_rows(document):
    rows = []
    for value in collection_objects(document):
        identifier = str(value.get("id") or "").strip()
        if not identifier:
            continue
        title = " ".join(str(
            value.get("title") or value.get("description") or "").split())
        if len(title) > 180:
            title = title[:177]+"..."
        rows.append(identifier+(" | "+title if title else ""))
    return rows


def row_identifier(row):
    return str(row or "").split(" | ", 1)[0].strip()


def find_collection(document, identifier):
    for value in collection_objects(document):
        if str(value.get("id") or "") == identifier:
            return value
    return None


def find_item(document, identifier):
    for value in feature_objects(document):
        if str(value.get("id") or "") == identifier:
            return value
    return None


def best_link(value):
    if not isinstance(value, dict):
        return ""
    preferred = ("self", "items", "data")
    links = [link for link in value.get("links") or [] if isinstance(link, dict)]
    for relation in preferred:
        for link in links:
            if str(link.get("rel") or "").lower() == relation:
                href = str(link.get("href") or "").strip()
                if href:
                    return href
    return ""


def describe_selection(document, kind, row):
    """Return accepted ID, URL and JSON for one stable browser row."""
    parts = [part.strip() for part in str(row or "").split(" | ")]
    identifier = parts[0] if parts else ""
    selected = None
    url = ""

    if kind == "collection":
        selected = find_collection(document, identifier)
        url = best_link(selected)
    elif kind == "item":
        selected = find_item(document, identifier)
        url = best_link(selected)
    elif kind == "asset":
        item_id = parts[0] if parts else ""
        asset_key = parts[1] if len(parts) > 1 else ""
        identifier = asset_key
        item = find_item(document, item_id)
        asset = (item.get("assets") or {}).get(asset_key) if item else None
        if isinstance(asset, dict):
            selected = {"item_id": item_id, "asset_key": asset_key, "asset": asset}
            url = str(asset.get("href") or "")
        elif len(parts) > 3:
            url = parts[3]

    return identifier, url, pretty(selected if selected is not None else str(row))


# -----------------------------------------------------------------------------
# Catalogue model
# -----------------------------------------------------------------------------

class CatalogueModel:
    """Single owner of browseable data; contains no Rhino or Eto operations."""

    KINDS = ("collection", "item", "asset")

    def __init__(self):
        self.document = {}
        self.context = {}
        self.rows = {kind: [] for kind in self.KINDS}
        self.log = []
        self.origin = "empty"
        self.load_starter()

    def load_starter(self, keep_non_collections=False):
        base = self.document if keep_non_collections else {}
        self.document = with_collections(base, starter_document())
        self.rows["collection"] = collection_rows(self.document)
        if not keep_non_collections:
            self.rows["item"] = []
            self.rows["asset"] = []
        self.log = [
            "LOCAL STARTER REGISTER",
            "No network request was made.",
            "Type to filter locally or press Search Copernicus Live.",
        ]
        self.origin = "starter"

    def load_inputs(self, data, collections, items, assets, fetch_log,
                    explicit_context=""):
        parsed = parse_json(data)
        self.context = (
            search_context(parse_json(explicit_context)) or search_context(parsed))
        base = dict(parsed)
        if self.context:
            base["search_context"] = self.context
        self.rows["item"] = list(items)
        self.rows["asset"] = list(assets)
        self.log = list(fetch_log)

        if collections:
            self.document = base
            self.rows["collection"] = list(collections)
            self.origin = "broker"
            return

        cached = CATALOGUE_CACHE.get("document")
        source = cached if isinstance(cached, dict) else starter_document()
        self.document = with_collections(base, source)
        self.rows["collection"] = collection_rows(self.document)
        if cached:
            self.log = list(CATALOGUE_CACHE.get("trace") or self.log)
            self.origin = "cache"
        else:
            self.log = [
                "LOCAL STARTER REGISTER",
                "No network request was made.",
                "Type to filter locally or press Search Copernicus Live.",
            ]
            self.origin = "starter"

    def load_live(self, live_document, trace, query):
        self.document = with_collections(
            self.document, live_document, starter_document())
        self.document["catalogue_search"] = live_document
        self.rows["collection"] = collection_rows(self.document)
        self.log = list(trace)
        self.origin = "live"
        CATALOGUE_CACHE.clear()
        CATALOGUE_CACHE.update({
            "document": live_document,
            "trace": list(trace),
            "query": str(query or ""),
            "timestamp": time.time(),
        })

    def filtered(self, query):
        needle = str(query or "").strip().lower()
        result = {}
        for kind in self.KINDS:
            visible = []
            for row in self.rows[kind]:
                haystack = str(row).lower()
                value = None
                if kind == "collection":
                    value = find_collection(self.document, row_identifier(row))
                    if not context_matches(self.context, value, row):
                        continue
                    if needle and value is not None:
                        haystack += " "+pretty(value).lower()
                if not needle or needle in haystack:
                    visible.append(row)
            result[kind] = visible
        return result

    def live_query(self, manual_query):
        manual = str(manual_query or "").strip()
        if manual or not self.context:
            return manual
        return str(self.context.get("suggested_live_query") or "").strip()

    def context_note(self):
        return context_summary(self.context)

    def selection(self, kind, row):
        return describe_selection(self.document, kind, row)


# -----------------------------------------------------------------------------
# Bounded official catalogue client
# -----------------------------------------------------------------------------

class CatalogueSearchError(RuntimeError):
    def __init__(self, errors, trace):
        self.errors = list(errors)
        self.trace = list(trace)
        super().__init__("; ".join(self.errors) or "No catalogue source responded.")


class CatalogueClient:
    """Search approved Copernicus endpoints; never mutates UI or GH state."""

    @staticmethod
    def _read(url):
        parsed = urllib.parse.urlparse(url)
        if parsed.scheme != "https" or parsed.hostname not in ALLOWED_HOSTS:
            raise ValueError("Catalogue URL is outside the approved official hosts.")
        request = urllib.request.Request(url, headers={
            "Accept": (
                "application/geo+json, application/json, application/xml, "
                "text/xml"),
            "User-Agent": "MusselFlow-CopernicusBrowser/%s" % COMPONENT_BUILD,
        })
        with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_S) as response:
            length = response.headers.get("Content-Length")
            if length and int(length) > MAX_RESPONSE_BYTES:
                raise ValueError("Catalogue response exceeds 20 MB.")
            payload = response.read(MAX_RESPONSE_BYTES+1)
        if len(payload) > MAX_RESPONSE_BYTES:
            raise ValueError("Catalogue response exceeds 20 MB.")
        return payload

    @staticmethod
    def _next_stac_url(current_url, document):
        for link in document.get("links", []) if isinstance(document, dict) else []:
            if not isinstance(link, dict):
                continue
            if str(link.get("rel") or "").lower() != "next":
                continue
            if str(link.get("method") or "GET").upper() != "GET":
                return ""
            candidate = urllib.parse.urljoin(
                current_url, str(link.get("href") or ""))
            parsed = urllib.parse.urlparse(candidate)
            return candidate if (
                parsed.scheme == "https" and parsed.hostname == CDSE_HOST) else ""
        return ""

    @classmethod
    def search_cdse(cls, query):
        clean = str(query or "").strip()
        url = CDSE_URL
        if clean:
            url += "?"+urllib.parse.urlencode({"q": clean})
        seen_urls = set()
        seen_ids = set()
        products = []
        requests = []
        byte_count = 0

        while url and len(requests) < MAX_STAC_PAGES:
            if url in seen_urls:
                break
            seen_urls.add(url)
            requests.append(url)
            payload = cls._read(url)
            byte_count += len(payload)
            page = json.loads(payload.decode("utf-8"))
            values = page.get("collections") if isinstance(page, dict) else None
            if not isinstance(values, list):
                raise ValueError("STAC response has no collections list.")
            for value in values:
                if not isinstance(value, dict):
                    continue
                identifier = str(value.get("id") or "").strip()
                if not identifier or identifier in seen_ids:
                    continue
                seen_ids.add(identifier)
                products.append(value)
                if len(products) >= MAX_COLLECTIONS:
                    url = ""
                    break
            else:
                url = cls._next_stac_url(url, page)

        document = {
            "mode": "stac_collections",
            "endpoint": CDSE_URL,
            "query": clean,
            "request_urls": requests,
            "response": {"collections": products},
        }
        trace = [
            "SOURCE | official CDSE STAC",
            "QUERY | %s" % (clean or "all collections"),
            "HTTP | %d page(s) | %.3f MB" % (
                len(requests), byte_count/(1024.0*1024.0)),
            "RESULT | %d collection(s)" % len(products),
        ]
        return document, trace

    @staticmethod
    def _element_text(record, path, namespaces):
        node = record.find(path, namespaces)
        return str(node.text or "").strip() if node is not None else ""

    @staticmethod
    def _marine_id(datasets, fallback):
        for value in datasets:
            parts = [part for part in urllib.parse.urlparse(
                str(value.get("href") or "")).path.split("/") if part]
            for marker in ("metadata", "teroWmts"):
                if marker in parts and parts.index(marker)+1 < len(parts):
                    return parts[parts.index(marker)+1]
        return fallback

    @classmethod
    def parse_marine(cls, payload):
        root = ET.fromstring(payload)
        ns = {
            "csw": "http://www.opengis.net/cat/csw/2.0.2",
            "dc": "http://purl.org/dc/elements/1.1/",
            "dct": "http://purl.org/dc/terms/",
        }
        records = list(root.findall(".//csw:Record", ns))
        records += list(root.findall(".//csw:SummaryRecord", ns))
        products = []

        for record in records:
            csw_id = cls._element_text(record, "dc:identifier", ns)
            title = cls._element_text(record, "dc:title", ns)
            description = cls._element_text(record, "dct:abstract", ns)
            if not description:
                description = cls._element_text(record, "dc:description", ns)
            keywords = [
                str(node.text or "").strip()
                for node in record.findall("dc:subject", ns)
                if str(node.text or "").strip()]
            datasets = []
            for node in record.findall("dc:URI", ns):
                href = str(node.text or "").strip()
                if href:
                    datasets.append({
                        "href": href,
                        "protocol": str(node.get("protocol") or ""),
                        "dataset_id": str(node.get("name") or ""),
                        "description": str(node.get("description") or ""),
                    })
            product_id = cls._marine_id(datasets, csw_id)
            products.append({
                "id": product_id,
                "title": "[Copernicus Marine] "+(title or product_id),
                "description": description,
                "keywords": keywords,
                "source": "Copernicus Marine CSW",
                "csw_identifier": csw_id,
                "datasets": datasets,
                "links": [{
                    "href": value["href"],
                    "rel": "data",
                    "type": value["protocol"],
                    "title": value["dataset_id"] or value["description"],
                } for value in datasets],
            })
        return products

    @classmethod
    def search_marine(cls, query):
        clean = str(query or "").strip()
        parameters = {
            "service": "CSW",
            "request": "GetRecords",
            "version": "2.0.2",
            "resultType": "results",
            "ElementSetName": "full" if clean else "summary",
            "maxRecords": 500,
            "startPosition": 1,
            "typeNames": "csw:Record",
        }
        if clean:
            escaped = clean.replace("'", "''")
            parameters.update({
                "constraintLanguage": "CQL_TEXT",
                "constraint_language_version": "1.1.0",
                "constraint": "AnyText LIKE '%%%s%%'" % escaped,
            })
        url = MARINE_URL+"?"+urllib.parse.urlencode(parameters)
        payload = cls._read(url)
        products = cls.parse_marine(payload)
        document = {
            "mode": "marine_csw_collections",
            "endpoint": MARINE_URL,
            "query": clean,
            "request_urls": [url],
            "response": {"collections": products},
        }
        trace = [
            "SOURCE | official Copernicus Marine CSW",
            "QUERY | %s" % (clean or "all products"),
            "HTTP | 1 request | %.3f MB" % (len(payload)/(1024.0*1024.0)),
            "RESULT | %d marine product(s)" % len(products),
        ]
        return document, trace

    @classmethod
    def search(cls, query):
        sources = (
            ("CDSE STAC", cls.search_cdse),
            ("Copernicus Marine CSW", cls.search_marine),
        )
        results = {}
        errors = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            jobs = {pool.submit(function, query): name for name, function in sources}
            for job in concurrent.futures.as_completed(jobs):
                name = jobs[job]
                try:
                    results[name] = job.result()
                except Exception as exception:
                    errors[name] = "%s: %s" % (
                        type(exception).__name__, str(exception))

        documents = []
        trace = []
        for name, _function in sources:
            if name in results:
                document, source_trace = results[name]
                documents.append(document)
                trace.extend(source_trace)
            else:
                trace.append("SOURCE ERROR | %s | %s" % (name, errors[name]))
        if not documents:
            messages = ["%s | %s" % (name, errors[name]) for name, _ in sources]
            raise CatalogueSearchError(messages, trace)

        combined = with_collections({}, *documents)
        products = collection_objects(combined)
        products.sort(key=lambda value: (
            str(value.get("title") or "").lower(),
            str(value.get("id") or "").lower()))
        envelope = {
            "mode": "copernicus_catalogue_search",
            "query": str(query or "").strip(),
            "sources": documents,
            "source_errors": [
                "%s | %s" % (name, errors[name]) for name in errors],
            "response": {"collections": products},
        }
        trace.append("COMBINED RESULT | %d collection/product(s)" % len(products))
        return envelope, trace


# -----------------------------------------------------------------------------
# Rhino / Grasshopper / Eto infrastructure
# -----------------------------------------------------------------------------

def apply_metadata(component):
    if component is None:
        return
    component.Name, component.NickName, component.Description = COMPONENT_METADATA
    component.Message = "Eto browser"
    for index, (name, nickname, description) in enumerate(INPUT_METADATA):
        if index >= component.Params.Input.Count:
            break
        parameter = component.Params.Input[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description
        if index in (2, 3, 4, 5):
            parameter.Access = Grasshopper.Kernel.GH_ParamAccess.list
        if index > 0:
            parameter.Optional = True
    for index, (name, nickname, description) in enumerate(OUTPUT_METADATA):
        if index >= component.Params.Output.Count:
            break
        parameter = component.Params.Output[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description
        if index == 5:
            parameter.Access = Grasshopper.Kernel.GH_ParamAccess.list


def component_id(component):
    try:
        return str(component.InstanceGuid)
    except Exception:
        return str(id(component))


def schedule_component(component):
    """Expire this component and ask its GH document for one future solution."""
    if component is None:
        return False, "component reference is unavailable"
    try:
        document = component.OnPingDocument()
        if document is not None:
            # Expire first, then use the simple ScheduleSolution(Int32)
            # overload. This avoids a short-lived Python callback/delegate
            # being collected before Grasshopper invokes it.
            component.ExpireSolution(False)
            document.ScheduleSolution(1)
            return True, "scheduled through GH_Document"
    except Exception as exception:
        schedule_error = "%s: %s" % (
            type(exception).__name__, str(exception))
    else:
        schedule_error = "component is not attached to a GH_Document"
    try:
        component.ExpireSolution(True)
        return True, "expired directly"
    except Exception as exception:
        return False, "%s | direct expiry failed: %s: %s" % (
            schedule_error, type(exception).__name__, str(exception))


def on_ui(action):
    callback = System.Action(action)
    application = forms.Application.Instance
    if application is not None:
        application.AsyncInvoke(callback)
    else:
        Rhino.RhinoApp.InvokeOnUiThread(callback)


def label(text):
    control = forms.Label()
    control.Text = str(text)
    return control


def button(text):
    control = forms.Button()
    control.Text = str(text)
    return control


def list_item(text):
    item = forms.ListItem()
    item.Text = str(text)
    return item


def form_visible(form):
    if form is None:
        return False
    try:
        return bool(form.Visible)
    except Exception:
        return True


def show_modeless(form):
    result = {"method": "RhinoEtoApp owner + Form.Show", "error": ""}
    try:
        form.Owner = Rhino.UI.RhinoEtoApp.MainWindow
        try:
            Grasshopper.GUI.GH_EtoUtil.CenterFormOnEditor(form, True)
        except Exception:
            pass
        form.Show()
    except Exception as exception:
        result["error"] = "%s: %s" % (
            type(exception).__name__, str(exception))
        return result
    try:
        form.BringToFront()
    except Exception:
        pass
    return result


# -----------------------------------------------------------------------------
# Eto presentation
# -----------------------------------------------------------------------------

class CopernicusBrowserForm(forms.Form):
    """Modeless view/controller. All data remains owned by CatalogueModel."""

    def __init__(self, title, session):
        super().__init__()
        self.Title = title
        self.ClientSize = drawing.Size(1040, 720)
        self.MinimumSize = drawing.Size(760, 500)
        self.Resizable = True
        self.ShowInTaskbar = False
        self.Padding = drawing.Padding(12)

        self.session = session
        self.model = CatalogueModel()
        self.visible_rows = {kind: [] for kind in CatalogueModel.KINDS}
        self.selected_kind = ""
        self.selected_row = ""
        self.search_token = 0
        self.searching = False
        self.closed = False

        self.search_box = forms.TextBox()
        self.search_box.PlaceholderText = (
            "Filter locally: phytoplankton, oxygen, currents, chlorophyll...")
        self.search_box.TextChanged += self._filter_changed
        self.live_button = button("Search Copernicus Live")
        self.live_button.Click += self._live_search
        self.starter_button = button("Starter Products")
        self.starter_button.Click += self._show_starter
        self.status = label("Loading starter products...")
        self.help = label(
            "Choose a product, click Use Selection, then connect FetchRequestJson to MusselFlow Site Data.")

        self.lists = {}
        self.tabs = forms.TabControl()
        for kind, title_text in (
                ("collection", "Collections"),
                ("item", "Items"),
                ("asset", "Assets")):
            control = forms.ListBox()
            control.Size = drawing.Size(640, 300)
            control.SelectedIndexChanged += self._selection_handler(kind)
            self.lists[kind] = control
            page = forms.TabPage()
            page.Text = title_text
            page.Content = control
            self.tabs.Pages.Add(page)

        self.raw = forms.TextArea()
        self.raw.ReadOnly = True
        self.raw.Wrap = False
        raw_page = forms.TabPage()
        raw_page.Text = "Raw Data"
        raw_page.Content = self.raw
        self.tabs.Pages.Add(raw_page)

        self.log = forms.TextArea()
        self.log.ReadOnly = True
        self.log.Wrap = False
        log_page = forms.TabPage()
        log_page.Text = "Fetch Log"
        log_page.Content = self.log
        self.tabs.Pages.Add(log_page)

        self.details = forms.TextArea()
        self.details.ReadOnly = True
        self.details.Wrap = False
        self.details.Height = 150

        self.use_button = button("Use Selection")
        self.use_button.Enabled = False
        self.use_button.Click += self._accept
        self.copy_button = button("Copy Details")
        self.copy_button.Click += self._copy
        self.close_button = button("Close")
        self.close_button.Click += self._close

        footer = forms.DynamicLayout()
        footer.Spacing = drawing.Size(8, 8)
        footer.AddRow(self.use_button, self.copy_button, None, self.close_button)

        layout = forms.DynamicLayout()
        layout.Padding = drawing.Padding(5)
        layout.Spacing = drawing.Size(8, 8)
        layout.AddRow(
            label("Search"), self.search_box, self.live_button,
            self.starter_button)
        layout.AddRow(self.status)
        layout.AddRow(self.help)
        layout.Add(self.tabs, True, True)
        layout.AddRow(label("Selected catalogue metadata; Use Selection builds the site request"))
        layout.Add(self.details, True, False)
        layout.AddRow(footer)
        self.Content = layout
        self.Closed += self._closed

    def load_inputs(self, data, collections, items, assets, fetch_log,
                    search_context=""):
        self.model.load_inputs(
            data, collections, items, assets, fetch_log, search_context)
        self._clear_selection()
        self._render()
        counts = [len(self.model.rows[kind]) for kind in CatalogueModel.KINDS]
        if self.model.context:
            self.status.Text = (
                "CONTEXT FILTER | %d matching products | %s"
                % (len(self.visible_rows["collection"]), self.model.context_note()))
        elif self.model.origin == "starter":
            self.status.Text = (
                "STARTER REGISTER | 6 verified products | no network request")
        elif self.model.origin == "cache":
            self.status.Text = "SESSION CACHE | %d products" % counts[0]
        else:
            self.status.Text = (
                "BROKER PAYLOAD | %d collections | %d items | %d assets"
                % tuple(counts))
        self.session.catalogue_view(
            self.model.origin, counts[0],
            len(self.visible_rows["collection"]))

    def _render(self):
        self.visible_rows = self.model.filtered(self.search_box.Text)
        for kind in CatalogueModel.KINDS:
            control = self.lists[kind]
            control.Items.Clear()
            for row in self.visible_rows[kind]:
                control.Items.Add(list_item(row))
        self.raw.Text = pretty(self.model.document)
        self.log.Text = "\n".join(self.model.log)

    def _clear_selection(self):
        self.selected_kind = ""
        self.selected_row = ""
        self.details.Text = ""
        self.use_button.Enabled = False

    def _selection_handler(self, kind):
        def selected(_sender, _event):
            index = int(self.lists[kind].SelectedIndex)
            rows = self.visible_rows[kind]
            if index < 0 or index >= len(rows):
                return
            self.selected_kind = kind
            self.selected_row = rows[index]
            identifier, url, detail = self.model.selection(kind, self.selected_row)
            self.details.Text = detail
            self.use_button.Enabled = True
            self.status.Text = "%s | %s%s" % (
                kind.upper(), identifier, " | "+url if url else "")
        return selected

    def _filter_changed(self, _sender, _event):
        self._clear_selection()
        self._render()
        query = str(self.search_box.Text or "").strip()
        total = len(self.model.rows["collection"])
        shown = len(self.visible_rows["collection"])
        if query:
            self.status.Text = (
                "LOCAL FILTER | %d of %d products | Search Live for more"
                % (shown, total))
        elif self.model.context:
            self.status.Text = (
                "CONTEXT FILTER | %d of %d products | %s"
                % (shown, total, self.model.context_note()))
        self.session.catalogue_view(self.model.origin, total, shown)

    def _show_starter(self, _sender, _event):
        self.search_token += 1
        self.model.load_starter(keep_non_collections=True)
        self.search_box.Text = ""
        self._clear_selection()
        self._render()
        self.status.Text = "STARTER REGISTER | 6 verified products"
        self.session.catalogue_view(
            self.model.origin, len(self.model.rows["collection"]),
            len(self.visible_rows["collection"]))

    def _live_search(self, _sender, _event):
        if self.searching:
            return
        query = self.model.live_query(self.search_box.Text)
        self.searching = True
        self.search_token += 1
        token = self.search_token
        self.live_button.Enabled = False
        self.live_button.Text = "Searching..."
        self.status.Text = "CONTACTING COPERNICUS | Rhino remains interactive"

        def worker():
            started = time.perf_counter()
            try:
                document, trace = CatalogueClient.search(query)
                error = None
            except Exception as exception:
                document = None
                trace = list(getattr(exception, "trace", []))
                error = "%s: %s" % (type(exception).__name__, str(exception))
            elapsed = time.perf_counter()-started

            def finish():
                if self.closed or token != self.search_token:
                    return
                self.searching = False
                self.live_button.Enabled = True
                self.live_button.Text = "Search Copernicus Live"
                if error is None:
                    self.model.load_live(document, trace, query)
                    self._clear_selection()
                    self._render()
                    shown = len(self.visible_rows["collection"])
                    total = len(self.model.rows["collection"])
                    self.status.Text = (
                        "LIVE CATALOGUE | %d products | %d matching | %.2fs"
                        % (total, shown, elapsed))
                    self.session.catalogue_view("live", total, shown)
                    self.session.catalogue_updated(query, total, elapsed, "")
                else:
                    self.model.log = trace + [
                        "LIVE CATALOGUE ERROR | "+error,
                        "Existing starter/cached results were preserved.",
                        "Retry when Rhino has DNS/network access.",
                        "CDSE | "+CDSE_URL,
                        "MARINE | "+MARINE_URL,
                    ]
                    self.log.Text = "\n".join(self.model.log)
                    total = len(self.model.rows["collection"])
                    shown = len(self.visible_rows["collection"])
                    self.status.Text = (
                        "LIVE SEARCH UNAVAILABLE | %d local products | "
                        "%d matching" % (total, shown))
                    self.session.catalogue_updated(query, 0, elapsed, error)

            on_ui(finish)

        thread = threading.Thread(target=worker)
        thread.daemon = True
        thread.start()

    def _accept(self, _sender, _event):
        if not self.selected_kind or not self.selected_row:
            self.status.Text = "Select one collection, item or asset first."
            return
        identifier, url, detail = self.model.selection(
            self.selected_kind, self.selected_row)
        self.session.accept(self.selected_kind, identifier, url, detail)
        self.status.Text = "SENT TO GRASSHOPPER | %s | %s | site request ready" % (
            self.selected_kind.upper(), identifier)

    def _copy(self, _sender, _event):
        value = str(self.details.Text or self.raw.Text or "")
        try:
            forms.Clipboard.Instance.Text = value
            self.status.Text = "Copied details to the clipboard."
        except Exception as exception:
            self.status.Text = "Clipboard error: %s" % exception

    def _close(self, _sender, _event):
        self.Close()

    def _closed(self, _sender, _event):
        self.closed = True
        self.search_token += 1
        try:
            self.session.form_closed(self)
        finally:
            try:
                self.Dispose()
            except Exception:
                pass


# -----------------------------------------------------------------------------
# One sticky session per Grasshopper component
# -----------------------------------------------------------------------------

class BrowserSession:
    def __init__(self, component):
        self.build = COMPONENT_BUILD
        self.component = component
        self.form = None
        self.last_open = False
        self.signature = ""
        self.payload = ("", [], [], [], [], "")
        self.selection_type = ""
        self.selection_id = ""
        self.selection_url = ""
        self.selection_json = ""
        self.current_context = {}
        self.ui_method = "not requested"
        self.ui_error = ""
        self.refresh_state = "not requested"
        self.catalogue_origin = "starter"
        self.catalogue_count = len(STARTER_PRODUCTS)
        self.visible_count = len(STARTER_PRODUCTS)
        self.query = ""
        self.result_count = 0
        self.search_seconds = 0.0
        self.search_error = ""
        self.context_note = "none"

    def dispose(self):
        form = self.form
        self.form = None
        if form is not None:
            try:
                form.Close()
            except Exception:
                try:
                    form.Dispose()
                except Exception:
                    pass

    def sync(self, data, collections, items, assets, fetch_log, context=""):
        signature = payload_hash(
            data, collections, items, assets, fetch_log, context)
        changed = signature != self.signature
        self.payload = (data, collections, items, assets, fetch_log, context)
        self.context_note = context_summary(
            search_context(parse_json(context)) or search_context(parse_json(data)))
        self.current_context = (
            search_context(parse_json(context)) or search_context(parse_json(data)))
        if changed:
            self.selection_type = ""
            self.selection_id = ""
            self.selection_url = ""
            self.selection_json = ""
            self.refresh_state = "waiting for selection"
        if self.form is not None and changed:
            self.form.load_inputs(*self.payload)
        self.signature = signature

    def open(self, title):
        if self.form is None:
            try:
                self.form = CopernicusBrowserForm(title, self)
                self.form.load_inputs(*self.payload)
                result = show_modeless(self.form)
                self.ui_method = result["method"]
                self.ui_error = result["error"]
                if self.ui_error:
                    failed_form = self.form
                    self.form = None
                    try:
                        failed_form.Dispose()
                    except Exception:
                        pass
            except Exception as exception:
                self.form = None
                self.ui_method = "form construction"
                self.ui_error = "%s: %s" % (
                    type(exception).__name__, str(exception))
        elif not form_visible(self.form):
            result = show_modeless(self.form)
            self.ui_method = result["method"]
            self.ui_error = result["error"]
        elif not self.last_open:
            try:
                self.form.BringToFront()
                self.form.Focus()
            except Exception:
                pass

    def accept(self, kind, identifier, url, detail):
        self.selection_type = str(kind or "")
        self.selection_id = str(identifier or "")
        self.selection_url = str(url or "")
        self.selection_json = pretty(make_site_request(
            self.current_context, kind, identifier, url, detail))
        refreshed, note = schedule_component(self.component)
        self.refresh_state = ("OK | " if refreshed else "ERROR | ")+note

    def catalogue_view(self, origin, total, visible):
        self.catalogue_origin = str(origin or "unknown")
        self.catalogue_count = int(total)
        self.visible_count = int(visible)

    def catalogue_updated(self, query, count, elapsed, error):
        self.query = str(query or "")
        self.result_count = int(count)
        self.search_seconds = float(elapsed)
        self.search_error = str(error or "")
        refreshed, note = schedule_component(self.component)
        self.refresh_state = ("OK | " if refreshed else "ERROR | ")+note

    def form_closed(self, form):
        if self.form is form:
            self.form = None
        refreshed, note = schedule_component(self.component)
        self.refresh_state = ("OK | " if refreshed else "ERROR | ")+note

    def report(self, counts):
        return [
            "COPERNICUS DATA BROWSER | build %s | %s" % (
                COMPONENT_BUILD, "OPEN" if form_visible(self.form) else "CLOSED"),
            "INPUT PAYLOAD | %d external collections | %d items | %d assets | %d log events"
            % tuple(counts),
            "CATALOGUE VIEW | %s | %d available | %d matching context/filter" % (
                self.catalogue_origin, self.catalogue_count, self.visible_count),
            "SELECTION | %s | %s" % (
                self.selection_type or "none", self.selection_id or "none"),
            "OUTPUT REFRESH | %s" % self.refresh_state,
            "WINDOW | %s%s" % (
                self.ui_method or "unknown",
                " | "+self.ui_error if self.ui_error else ""),
            "CATALOGUE | %d live result(s) | query %s | %.3fs%s" % (
                self.result_count, self.query or "not requested",
                self.search_seconds,
                " | "+self.search_error if self.search_error else ""),
            "SEARCH CONTEXT | %s" % self.context_note,
            "OFFLINE | six verified Baltic/Northwest Shelf starter products.",
            "NETWORK | live search uses official CDSE STAC and Marine CSW only.",
            "DATA LEVEL | FetchRequestJson is intent + product metadata; data_state=not_fetched.",
            "NEXT | FetchRequestJson -> MusselFlow Site Data.FetchRequestJson; press fetch for actual values.",
            "MUSSELFLOW WIRING | Site Data.SiteDataJson -> Site Field.SiteDataJson; Site Data.FlowVectors -> Optimizer.flowVectors.",
            "GENERAL WIRING | SelectedId -> Universal Data Broker.collection for a compatible STAC/WMTS request.",
        ]


def retire_session(value):
    """Close sessions/forms left by an older pasted build."""
    if value is None:
        return
    try:
        value.dispose()
        return
    except Exception:
        pass
    if isinstance(value, dict):
        form = value.get("form")
        if form is not None:
            try:
                form.Close()
            except Exception:
                try:
                    form.Dispose()
                except Exception:
                    pass


def get_session(component):
    key = component_id(component)
    session = SESSIONS.get(key)
    if getattr(session, "build", None) != COMPONENT_BUILD:
        retire_session(session)
        session = BrowserSession(component)
        SESSIONS[key] = session
    session.component = component
    return session


# -----------------------------------------------------------------------------
# Grasshopper SDK entry point
# -----------------------------------------------------------------------------

class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def BeforeRunScript(self):
        apply_metadata(getattr(self, "Component", None))

    def RunScript(
            self,
            openBrowser: bool,
            data: str,
            collections: list[str],
            items: list[str],
            assets: list[str],
            fetchLog: list[str],
            windowTitle: str,
            searchContext: str):
        """Open/update the browser and expose its accepted site-data request."""
        component = getattr(self, "Component", None)
        session = get_session(component)

        collection_rows_input = text_list(collections)
        item_rows_input = text_list(items)
        asset_rows_input = text_list(assets)
        log_rows_input = text_list(fetchLog)
        session.sync(
            data, collection_rows_input, item_rows_input,
            asset_rows_input, log_rows_input, searchContext)

        requested = bool(openBrowser)
        if requested:
            session.open(str(windowTitle or "").strip() or DEFAULT_TITLE)
        session.last_open = requested

        counts = (
            len(collection_rows_input), len(item_rows_input),
            len(asset_rows_input), len(log_rows_input))
        return (
            session.selection_type,
            session.selection_id,
            session.selection_url,
            session.selection_json,
            form_visible(session.form),
            session.report(counts),
        )

    def AfterRunScript(self):
        pass
