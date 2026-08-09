"""Regression tests for the standalone Rhino/Grasshopper Copernicus browser."""

import ast
import importlib.util
import json
import sys
import types
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "copernicus_data_browser_gh_sdk.py"
CONTEXT_SOURCE = HERE / "copernicus_search_context_gh_sdk.py"
SITE_DATA_SOURCE = HERE / "musselflow_site_data_gh_sdk.py"


def _module_stub(name, **values):
    module = types.ModuleType(name)
    for key, value in values.items():
        setattr(module, key, value)
    sys.modules[name] = module
    return module


def load_browser_module():
    """Load pure browser logic without requiring Rhino on the test machine."""
    form_base = type("Form", (), {})
    drawing = _module_stub("Eto.Drawing")
    forms = _module_stub("Eto.Forms", Form=form_base)
    eto = _module_stub("Eto", Drawing=drawing, Forms=forms)
    eto.Drawing, eto.Forms = drawing, forms

    script_base = type("GH_ScriptInstance", (), {})
    kernel = types.SimpleNamespace(GH_ScriptInstance=script_base)
    _module_stub("Grasshopper", Kernel=kernel)
    _module_stub("Rhino")
    _module_stub("System")
    _module_stub("scriptcontext", sticky={})

    spec = importlib.util.spec_from_file_location("browser_under_test", SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_context_module():
    if "Grasshopper" not in sys.modules:
        script_base = type("GH_ScriptInstance", (), {})
        kernel = types.SimpleNamespace(GH_ScriptInstance=script_base)
        _module_stub("Grasshopper", Kernel=kernel)
    spec = importlib.util.spec_from_file_location("context_under_test", CONTEXT_SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_site_data_module():
    script_base = type("GH_ScriptInstance", (), {})
    kernel = types.SimpleNamespace(GH_ScriptInstance=script_base)
    _module_stub("Grasshopper", Kernel=kernel)
    vector = type("Vector3d", (), {})
    _module_stub("Rhino", Geometry=types.SimpleNamespace(Vector3d=vector))
    spec = importlib.util.spec_from_file_location(
        "site_data_under_test", SITE_DATA_SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class BrowserContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SOURCE.read_text(encoding="utf-8")
        cls.tree = ast.parse(cls.source)
        cls.browser = load_browser_module()
        cls.context_builder = load_context_module()
        cls.site_data = load_site_data_module()

    def test_sdk_signature_and_return_count_are_stable(self):
        script_class = next(
            node for node in self.tree.body
            if isinstance(node, ast.ClassDef) and node.name == "Script_Instance")
        run = next(
            node for node in script_class.body
            if isinstance(node, ast.FunctionDef) and node.name == "RunScript")
        self.assertEqual(
            [arg.arg for arg in run.args.args],
            ["self", "openBrowser", "data", "collections", "items", "assets",
             "fetchLog", "windowTitle", "searchContext"],
        )
        returns = [node for node in ast.walk(run) if isinstance(node, ast.Return)]
        public_return = returns[-1].value
        self.assertIsInstance(public_return, ast.Tuple)
        self.assertEqual(len(public_return.elts), 6)

    def test_starter_register_and_local_search(self):
        model = self.browser.CatalogueModel()
        self.assertEqual(len(model.rows["collection"]), 6)
        for query in ("phytoplankton", "oxygen", "currents", "chlorophyll",
                      "Baltic", "North Sea"):
            self.assertTrue(model.filtered(query)["collection"], query)

    def test_report_distinguishes_catalogue_metadata_from_site_values(self):
        session = self.browser.BrowserSession(None)
        report = "\n".join(session.report((0, 0, 0, 0)))
        self.assertIn("data_state=not_fetched", report)
        self.assertIn("Site Data.SiteDataJson -> Site Field.SiteDataJson", report)
        self.assertIn("Site Data.FlowVectors -> Optimizer.flowVectors", report)

    def test_accepted_selection_builds_one_complete_site_request(self):
        context, _query = self.context_builder.build_context(
            54.95117, 7.60938, "2026-04-20T00:00:00Z", "",
            [1, 5, 10], ["phosphate"], "", None, [])
        request = self.browser.make_site_request(
            context, "collection", "NWSHELF_ANALYSISFORECAST_BGC_004_002",
            "https://example.invalid/product", json.dumps({"title": "BGC"}))
        self.assertEqual(request["schema"], "copernicus.site_request.1.0")
        self.assertEqual(request["data_state"], "not_fetched")
        self.assertEqual(request["search_context"]["depths_m"], [1.0, 5.0, 10.0])
        self.assertEqual(
            request["selection"]["id"],
            "NWSHELF_ANALYSISFORECAST_BGC_004_002")
        self.assertEqual(request["selection"]["metadata"]["title"], "BGC")
        self.assertEqual(
            request["fetch_plan"]["requested_variables"], ["phosphate"])

    def test_site_data_accepts_browser_request_as_final_input(self):
        tree = ast.parse(SITE_DATA_SOURCE.read_text(encoding="utf-8"))
        script_class = next(
            node for node in tree.body
            if isinstance(node, ast.ClassDef) and node.name == "Script_Instance")
        run = next(
            node for node in script_class.body
            if isinstance(node, ast.FunctionDef) and node.name == "RunScript")
        self.assertEqual(run.args.args[-1].arg, "siteRequestJson")

    def test_site_request_overrides_manual_site_time_and_depths(self):
        context, _query = self.context_builder.build_context(
            54.95117, 7.60938, "2026-04-20T00:00:00Z", "",
            [1, 5, 10], ["phosphate"], "", None, [])
        request = self.browser.make_site_request(
            context, "collection", "NWSHELF_ANALYSISFORECAST_BGC_004_002",
            "https://example.invalid/product", "{}")
        parsed = self.site_data.parse_site_request(json.dumps(request))
        values = self.site_data.request_override(
            parsed, 0.0, 0.0, "manual-start", "manual-end", [99])
        self.assertEqual(values[:2], (54.95117, 7.60938))
        self.assertEqual(values[2], "2026-04-20T00:00:00Z")
        self.assertEqual(values[3], "manual-end")
        self.assertEqual(values[4], [1.0, 5.0, 10.0])

    def test_search_context_filters_region_and_variable_families(self):
        context, query = self.context_builder.build_context(
            54.95, 12.0, "2026-04-20T00:00:00Z",
            "2026-04-21T00:00:00Z", [1, 5],
            ["flow vectors", "oxygen"], "", 1000, [])
        self.assertEqual(query, "current")
        self.assertEqual(context["site"]["region_hint"]["name"], "Baltic Sea")
        self.assertEqual(context["variables"], ["current", "oxygen"])

        model = self.browser.CatalogueModel()
        model.load_inputs("", [], [], [], [], json.dumps(context))
        rows = model.filtered("")["collection"]
        self.assertTrue(any("BALTICSEA_ANALYSISFORECAST_PHY" in row for row in rows))
        self.assertTrue(any("BALTICSEA_ANALYSISFORECAST_BGC" in row for row in rows))
        self.assertFalse(any("NWSHELF" in row for row in rows))

    def test_nested_grasshopper_depth_tree_is_flattened(self):
        class Goo:
            def __init__(self, value):
                self.Value = value

        class Tree:
            Branches = [[Goo(1.0), Goo("5")], [Goo("10; 20")]]

        self.assertEqual(
            self.context_builder.normalized_depths(Tree()),
            [1.0, 5.0, 10.0, 20.0],
        )

    def test_preferred_product_hint_can_select_one_product(self):
        context, _query = self.context_builder.build_context(
            54.95, 12.0, "", "", [1], ["current"],
            "BALTICSEA_ANALYSISFORECAST_PHY_003_006/dataset/layer", None, [])
        model = self.browser.CatalogueModel()
        model.load_inputs("", [], [], [], [], json.dumps(context))
        rows = model.filtered("")["collection"]
        self.assertEqual(len(rows), 1)
        self.assertTrue(rows[0].startswith("BALTICSEA_ANALYSISFORECAST_PHY_003_006"))

    def test_merge_deduplicates_and_preserves_features(self):
        base = {"response": {"features": [{"id": "frame-1"}]}}
        source_a = {"response": {"collections": [{"id": "A"}]}}
        source_b = {"response": {"collections": [{"id": "A"}, {"id": "B"}]}}
        result = self.browser.with_collections(base, source_a, source_b)
        self.assertEqual([value["id"] for value in result["response"]["collections"]],
                         ["A", "B"])
        self.assertEqual(result["response"]["features"], [{"id": "frame-1"}])

    def test_marine_csw_parser_extracts_product(self):
        xml = b'''<csw:GetRecordsResponse
          xmlns:csw="http://www.opengis.net/cat/csw/2.0.2"
          xmlns:dc="http://purl.org/dc/elements/1.1/"
          xmlns:dct="http://purl.org/dc/terms/">
          <csw:SearchResults><csw:Record>
            <dc:identifier>fallback-id</dc:identifier>
            <dc:title>Phytoplankton product</dc:title>
            <dct:abstract>Example</dct:abstract>
            <dc:subject>phytoplankton</dc:subject>
            <dc:URI protocol="WWW:LINK-1.0-http--link">https://x/metadata/PRODUCT_ID</dc:URI>
          </csw:Record></csw:SearchResults>
        </csw:GetRecordsResponse>'''
        products = self.browser.CatalogueClient.parse_marine(xml)
        self.assertEqual(len(products), 1)
        self.assertEqual(products[0]["id"], "PRODUCT_ID")
        self.assertIn("phytoplankton", products[0]["keywords"])


if __name__ == "__main__":
    unittest.main()
