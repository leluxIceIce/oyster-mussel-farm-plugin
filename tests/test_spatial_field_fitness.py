"""Regression tests for spatial ecology in the MusselFlow fitness core.

These tests establish code-level causality: changing a raw environmental field
must change the biological quantity governed by that field.  They do not claim
that the reduced-order equations are calibrated predictions of a real farm.
"""

import ast
from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULES = ROOT / "Miesmuschel"
if str(MODULES) not in sys.path:
    sys.path.insert(0, str(MODULES))

import musselflow_bio_optimizer_core as optimizer  # noqa: E402
import musselflow_case_core as case_core  # noqa: E402


def compiled_case(obstacle_count=3):
    source = (MODULES / "musselflow_grammar.json").read_text(
        encoding="utf-8")
    case, _warnings = case_core.parse_case(
        source, flow_count=2, obstacle_count=obstacle_count)
    config, flow_indices, _warnings = case_core.compile_timeline(
        case, obstacle_count, 2)
    return config, flow_indices


def layout():
    return np.asarray([
        [0.0, 0.0, 1.2, 0.5, 0.0, 3.0],
        [4.0, 1.0, 1.2, 0.5, 0.2, 3.0],
        [8.0, -1.0, 1.2, 0.5, -0.2, 3.0],
    ], dtype=float)


DOMAIN = np.asarray([
    [-10.0, -10.0],
    [20.0, -10.0],
    [20.0, 10.0],
    [-10.0, 10.0],
], dtype=float)
PROBES = np.empty((0, 2), dtype=float)
FLOWS = np.asarray([[0.30, 0.0], [-0.20, 0.0]], dtype=float)


def field(name, values):
    return {
        name: {
            "values": np.asarray(values, dtype=float),
            "metadata": {"test_field": name},
        }
    }


class SpatialBiologyTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.obstacles = layout()
        cls.config, indices = compiled_case(len(cls.obstacles))
        cls.flows = FLOWS[np.asarray(indices, dtype=int)]

    def evaluate(self, environment_fields=None):
        return optimizer.evaluate_layout(
            self.obstacles, DOMAIN, PROBES, self.flows, self.config,
            environment_fields=environment_fields)

    def test_no_field_matches_explicit_case_baselines(self):
        baseline = self.evaluate()
        count = len(self.obstacles)
        explicit = {
            "chlorophyll_ug_l": np.full(
                count, self.config["site.chlorophyll_ug_l"]),
            "tsm_mg_l": np.full(count, self.config["site.tsm_mg_l"]),
            "temperature_c": np.full(
                count, self.config["site.temperature_c"]),
            "salinity_psu": np.full(
                count, self.config["site.salinity_psu"]),
            "boundary_do_mg_l": np.full(
                count, self.config["site.boundary_do_mg_l"]),
        }
        conditioned = self.evaluate(explicit)
        for key in (
                "fitness", "chlorophyll_capture_g_day",
                "particulate_capture_kg_day",
                "mussel_respiration_kg_o2_day"):
            self.assertAlmostEqual(baseline[key], conditioned[key], places=12)
        self.assertEqual(
            baseline["environment"]["summary"]["source"],
            "SimulationCaseJson")
        self.assertEqual(
            conditioned["environment"]["summary"]["source"],
            "FieldDataJson")

    def test_chlorophyll_changes_capture(self):
        low = self.evaluate(field("chlorophyll_ug_l", [0.8, 0.8, 0.8]))
        high = self.evaluate(field("chlorophyll_ug_l", [8.0, 8.0, 8.0]))
        self.assertGreater(
            high["chlorophyll_capture_g_day"],
            low["chlorophyll_capture_g_day"])
        self.assertIn(
            "chlorophyll_ug_l",
            high["biology_trace"]["environment_fields"])

    def test_tsm_changes_pseudofaeces_and_biodeposition(self):
        low = self.evaluate(field("tsm_mg_l", [1.0, 1.0, 1.0]))
        high = self.evaluate(field("tsm_mg_l", [20.0, 20.0, 20.0]))
        self.assertGreater(
            np.mean(high["pseudofaeces_fraction_by_obstacle"]),
            np.mean(low["pseudofaeces_fraction_by_obstacle"]))
        self.assertGreater(
            high["biodeposit_organic_kg_day"],
            low["biodeposit_organic_kg_day"])

    def test_temperature_changes_clearance_and_respiration(self):
        cold = self.evaluate(field("temperature_c", [5.0, 5.0, 5.0]))
        warm = self.evaluate(field("temperature_c", [20.0, 20.0, 20.0]))
        self.assertGreater(
            warm["biology_trace"]["maximum_clearance_l_h_total"],
            cold["biology_trace"]["maximum_clearance_l_h_total"])
        self.assertGreater(
            warm["mussel_respiration_kg_o2_day"],
            cold["mussel_respiration_kg_o2_day"])

    def test_salinity_changes_clearance_activity(self):
        marginal = self.evaluate(field("salinity_psu", [5.0, 5.0, 5.0]))
        optimal = self.evaluate(field("salinity_psu", [25.0, 25.0, 25.0]))
        self.assertGreater(
            np.mean(optimal["salinity_activity_by_obstacle"]),
            np.mean(marginal["salinity_activity_by_obstacle"]))
        self.assertGreater(
            optimal["biology_trace"]["maximum_clearance_l_h_total"],
            marginal["biology_trace"]["maximum_clearance_l_h_total"])

    def test_oxygen_changes_activity_and_boundary_balance(self):
        low = self.evaluate(field("boundary_do_mg_l", [1.0, 1.0, 1.0]))
        high = self.evaluate(field("boundary_do_mg_l", [9.0, 9.0, 9.0]))
        self.assertGreater(
            np.mean(high["oxygen_activity_by_obstacle"]),
            np.mean(low["oxygen_activity_by_obstacle"]))
        self.assertAlmostEqual(low["oxygen"]["boundary_mg_l"], 1.0)
        self.assertAlmostEqual(high["oxygen"]["boundary_mg_l"], 9.0)

    def test_raw_environment_changes_final_galapagos_fitness(self):
        low_oxygen = self.evaluate(field(
            "boundary_do_mg_l", [1.0, 1.0, 1.0]))
        high_oxygen = self.evaluate(field(
            "boundary_do_mg_l", [9.0, 9.0, 9.0]))
        self.assertGreater(high_oxygen["fitness"], low_oxygen["fitness"])

        low_food = self.evaluate(field(
            "chlorophyll_ug_l", [0.8, 0.8, 0.8]))
        high_food = self.evaluate(field(
            "chlorophyll_ug_l", [8.0, 8.0, 8.0]))
        self.assertNotAlmostEqual(
            low_food["fitness"], high_food["fitness"], places=8)

    def test_biology_trace_is_numerically_tied_to_result(self):
        result = self.evaluate({
            "chlorophyll_ug_l": np.asarray([1.0, 4.0, 8.0]),
            "tsm_mg_l": np.asarray([2.0, 5.0, 12.0]),
        })
        trace = result["biology_trace"]
        self.assertTrue(trace["executed"])
        self.assertGreaterEqual(len(trace["equations"]), 9)
        self.assertAlmostEqual(
            trace["chlorophyll_capture_g_day"],
            result["chlorophyll_capture_g_day"])
        self.assertAlmostEqual(
            trace["mussel_respiration_kg_o2_day"],
            result["mussel_respiration_kg_o2_day"])
        self.assertAlmostEqual(trace["fitness"], result["fitness"])

    def test_biology_trace_separates_evidence_from_proxies(self):
        evidence = self.evaluate()["biology_trace"]["evidence"]
        self.assertIn(
            "clearance_allometry_mohlenberg_riisgard_1979",
            evidence["source_backed"])
        self.assertIn(
            "one_box_oxygen_mass_balance",
            evidence["mechanistic_reduced_order"])
        self.assertIn(
            "current_activity_curve",
            evidence["calibration_proxies"])
        self.assertIn(
            "normalized_weighted_objective",
            evidence["decision_logic"])
        self.assertNotIn(
            "normalized_weighted_objective",
            evidence["source_backed"])


class GrasshopperSpatialFieldContractTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.path = MODULES / "musselflow_component_gh_sdk.py"
        cls.source = cls.path.read_text(encoding="utf-8")
        cls.tree = ast.parse(cls.source, filename=str(cls.path))
        cls.run_script = next(
            node for node in ast.walk(cls.tree)
            if isinstance(node, ast.FunctionDef) and node.name == "RunScript")

    def test_field_input_is_appended_without_reordering_existing_ports(self):
        names = [argument.arg for argument in self.run_script.args.args]
        self.assertEqual(names, [
            "self", "run", "obstacles", "domain", "probes", "flowVectors",
            "SimulationCaseJson", "qualityMode", "speedMode",
            "FieldDataJson",
        ])
        annotation = ast.unparse(self.run_script.args.args[-1].annotation)
        self.assertIn("Grasshopper.DataTree", annotation)

    def test_site_field_record_accepts_current_xyz_object_contract(self):
        helper = next(
            node for node in ast.walk(self.tree)
            if isinstance(node, ast.FunctionDef)
            and node.name == "_field_point_xy")
        helper_source = ast.unparse(helper)
        self.assertIn("isinstance(coordinates, dict)", helper_source)
        self.assertIn("coordinates.get('x')", helper_source)

    def test_solver_reads_raw_field_value_not_preview_normalization(self):
        resolver = next(
            node for node in ast.walk(self.tree)
            if isinstance(node, ast.FunctionDef)
            and node.name == "resolve_environment_fields")
        resolver_source = ast.unparse(resolver)
        self.assertIn("record.get('value')", resolver_source)
        self.assertNotIn("record.get('normalized')", resolver_source)

    def test_field_values_enter_solver_without_replacing_case_json(self):
        run_source = ast.unparse(self.run_script)
        self.assertIn(
            "resolve_environment_fields(FieldDataJson", run_source)
        self.assertIn(
            "compile_case_cached(case_core, SimulationCaseJson", run_source)
        self.assertIn("environment_fields=environment_fields", run_source)


if __name__ == "__main__":
    unittest.main()
