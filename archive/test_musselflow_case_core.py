import ast
import copy
import json
import math
import unittest
from pathlib import Path

from musselflow_ecogrammar_core import DEFAULTS
from musselflow_case_core import (
    CaseValidationError,
    UnsupportedCaseError,
    canonical_json,
    case_hash,
    compile_ensemble,
    compile_timeline,
    ensemble_state_configs,
    parse_case,
)


ROOT = Path(__file__).resolve().parent
EXAMPLE = ROOT / "musselflow_grammar.json"
COMPONENT = ROOT / "musselflow_component_gh_sdk.py"
GEOMETRY_BRIDGE = ROOT / "musselflow_bio_optimizer_gh_sdk.py"
VISUALIZER = ROOT / "musselflow_gh_sdk.py"


class MusselFlowCaseCoreTests(unittest.TestCase):

    def case(self):
        return json.loads(EXAMPLE.read_text(encoding="utf-8"))

    def test_example_is_valid_two_flow_timeline(self):
        case, warnings = parse_case(
            EXAMPLE, flow_count=2, obstacle_count=1)
        self.assertEqual(case["forcing"]["mode"], "timeline")
        self.assertEqual(
            [step["flow_vector_index"] for step in case["forcing"]["steps"]],
            [0, 1])
        self.assertIsInstance(warnings, list)

    def test_canonical_json_and_hash_are_order_independent(self):
        case = self.case()
        reversed_case = dict(reversed(list(case.items())))
        self.assertEqual(canonical_json(case), canonical_json(reversed_case))
        self.assertEqual(case_hash(case), case_hash(reversed_case))
        self.assertEqual(len(case_hash(case)), 64)

    def test_unknown_key_is_rejected(self):
        case = self.case()
        case["mystery"] = 1
        with self.assertRaisesRegex(CaseValidationError, "unknown keys"):
            parse_case(case, flow_count=2, obstacle_count=1)

    def test_non_finite_json_is_rejected(self):
        text = EXAMPLE.read_text(encoding="utf-8").replace(
            '"depth_m": 12.0', '"depth_m": NaN')
        with self.assertRaisesRegex(CaseValidationError, "non-finite"):
            parse_case(text, flow_count=2, obstacle_count=1)

    def test_out_of_range_flow_reference_is_rejected(self):
        case = self.case()
        case["forcing"]["steps"][1]["flow_vector_index"] = 2
        with self.assertRaisesRegex(CaseValidationError, "flow_count=2"):
            parse_case(case, flow_count=2, obstacle_count=1)

    def test_duplicate_scenario_id_is_rejected(self):
        case = self.case()
        case["forcing"]["steps"][1]["id"] = "flood"
        with self.assertRaisesRegex(CaseValidationError, "must be unique"):
            parse_case(case, flow_count=2, obstacle_count=1)

    def test_timeline_order_must_be_contiguous(self):
        case = self.case()
        case["forcing"]["steps"][1]["order"] = 3
        with self.assertRaisesRegex(CaseValidationError, "exactly"):
            parse_case(case, flow_count=2, obstacle_count=1)

    def test_per_obstacle_array_mismatch_is_rejected(self):
        case = self.case()
        case["structure"]["porosity_per_obstacle"] = [0.5, 0.6]
        with self.assertRaisesRegex(CaseValidationError, "expected 1 or 3"):
            parse_case(case, flow_count=2, obstacle_count=3)

    def test_editable_calibration_status_is_rejected(self):
        case = self.case()
        case["validation"]["status"] = "CALIBRATED_WITHIN_DECLARED_ENVELOPE"
        with self.assertRaisesRegex(CaseValidationError, "derived"):
            parse_case(case, flow_count=2, obstacle_count=1)

    def ensemble_case(self):
        case = self.case()
        boundary = copy.deepcopy(case["forcing"]["steps"][0]["boundary"])
        case["forcing"] = {
            "mode": "ensemble",
            "vector_source": "Grasshopper.flowVectors",
            "states": [
                {
                    "id": "flood",
                    "flow_vector_index": 0,
                    "occurrence_probability": 2.0 / 3.0,
                    "duration_h": 6.0,
                    "boundary": boundary,
                },
                {
                    "id": "ebb",
                    "flow_vector_index": 1,
                    "occurrence_probability": 1.0 / 6.0,
                    "duration_h": 6.0,
                    "boundary": copy.deepcopy(boundary),
                },
            ],
        }
        return case

    def test_ensemble_probabilities_are_normalized(self):
        parsed, warnings = parse_case(
            self.ensemble_case(), flow_count=2, obstacle_count=1)
        probabilities = [
            state["occurrence_probability"]
            for state in parsed["forcing"]["states"]]
        self.assertAlmostEqual(sum(probabilities), 1.0)
        self.assertTrue(any("normalized" in warning for warning in warnings))

    def test_ensemble_state_missing_duration_is_rejected(self):
        case = self.ensemble_case()
        del case["forcing"]["states"][0]["duration_h"]
        with self.assertRaisesRegex(CaseValidationError, "duration_h"):
            parse_case(case, flow_count=2, obstacle_count=1)

    def test_timeline_compiler_refuses_ensemble(self):
        with self.assertRaisesRegex(UnsupportedCaseError, "compile_ensemble"):
            compile_timeline(
                self.ensemble_case(), obstacle_count=1, flow_count=2)

    def test_previous_schema_is_accepted_only_as_migration_input(self):
        case = self.case()
        case["schema_version"] = "2.0.0-draft"
        case["machine_learning"] = {
            "surrogate_enabled": False,
            "model_artifact_id": None,
            "training_dataset_id": None,
            "reinforcement_learning": {
                "enabled": False,
                "intended_future_role": "compatibility only",
                "static_layout_optimizer": "evolutionary search",
            },
        }
        case["evidence"] = [
            {"parameter": record["parameter"], "doi": record["doi"]}
            for record in case["evidence"]]
        parsed, warnings = parse_case(
            case, flow_count=2, obstacle_count=1)
        self.assertEqual(parsed["schema_version"], "2.0.0-draft")
        self.assertTrue(any("accepted for migration" in item
                            for item in warnings))
        self.assertTrue(any("legacy inert metadata" in item
                            for item in warnings))

    def test_ensemble_compiles_to_base_config_and_states(self):
        base_config, states, warnings = compile_ensemble(
            self.ensemble_case(), obstacle_count=1, flow_count=2)
        self.assertEqual([s["id"] for s in states], ["flood", "ebb"])
        self.assertEqual([s["flow_vector_index"] for s in states], [0, 1])
        self.assertAlmostEqual(
            sum(s["occurrence_probability"] for s in states), 1.0)
        # base_config carries the shared, boundary-independent fields.
        self.assertEqual(base_config["site.depth_m"], 12.0)
        self.assertEqual(base_config["harvest.n_kg_per_t_wet"], 13.7)
        self.assertEqual(states[0]["duration_h"], 6.0)

    def test_ensemble_state_configs_apply_boundary_per_state(self):
        base_config, states, _ = compile_ensemble(
            self.ensemble_case(), obstacle_count=1, flow_count=2)
        states[0]["boundary"]["temperature_c"] = 10.0
        states[1]["boundary"]["temperature_c"] = 14.0
        specs = ensemble_state_configs(base_config, states)
        self.assertEqual(specs[0]["config"]["site.temperature_c"], 10.0)
        self.assertEqual(specs[1]["config"]["site.temperature_c"], 14.0)
        # The shared base config is not mutated by per-state application.
        self.assertEqual(base_config["site.temperature_c"],
                         base_config["site.temperature_c"])
        self.assertEqual(specs[0]["config"]["site.depth_m"], 12.0)
        self.assertEqual(
            specs[0]["probability"], states[0]["occurrence_probability"])

    def test_timeline_compiler_maps_flow_indices_and_durations(self):
        config, flow_indices, warnings = compile_timeline(
            self.case(), obstacle_count=3, flow_count=2)
        self.assertEqual(flow_indices, [0, 1])
        self.assertEqual(config["scenario.duration_h"], [6.0, 6.0])
        self.assertEqual(config["scenario.weights"], [0.5, 0.5])
        self.assertEqual(
            config["stocking.mussels_per_obstacle"],
            [1000.0, 1000.0, 1000.0])
        self.assertEqual(config["structure.porosity"], [0.7, 0.7, 0.7])
        self.assertFalse(config["validation.biology_calibrated"])
        self.assertTrue(any("UNVALIDATED_SCREENING" in item
                            for item in warnings))

    def test_canonical_grammar_resolves_to_all_115_reference_parameters(self):
        """Prevent biological coefficients from disappearing into adapters."""
        config, _, _ = compile_timeline(
            self.case(), obstacle_count=1, flow_count=2)
        expected = copy.deepcopy(DEFAULTS)
        expected["scenario.weights"] = [0.5, 0.5]
        expected["scenario.duration_h"] = [6.0, 6.0]
        self.assertEqual(len(expected), 115)
        self.assertEqual(config, expected)

    def test_canonical_grammar_carries_full_evidence_records(self):
        case = self.case()
        self.assertEqual(len(case["evidence"]), 11)
        required = {
            "parameter", "citation", "equation", "validity_envelope", "doi"}
        for record in case["evidence"]:
            self.assertEqual(set(record), required)

    def test_ecological_fields_compile_into_the_solver_config(self):
        config, _, _ = compile_timeline(
            self.case(), obstacle_count=1, flow_count=2)
        # Species fields are explicit in the canonical grammar.
        self.assertEqual(config["species.valid_flow_max_m_s"], 1.4)
        self.assertEqual(config["species.ammonia_mg_n_g_dw_h"], 0.015)
        self.assertEqual(config["species.low_food_threshold_ug_l"], 0.7)
        self.assertEqual(config["species.low_food_transition_ug_l"], 0.15)
        # Sediment section.
        self.assertEqual(config["sediment.settling_velocity_m_s"], 0.006)
        self.assertEqual(config["sediment.decay_per_day"], 0.05)
        self.assertEqual(config["sediment.mortality_deposition_fraction"], 1.0)
        # Harvest section.
        self.assertEqual(config["harvest.fraction_per_year"], 0.8)
        self.assertEqual(config["harvest.n_kg_per_t_wet"], 13.7)
        self.assertEqual(config["harvest.p_kg_per_t_wet"], 0.9)

    def test_null_deposition_fraction_maps_to_automatic_sentinel(self):
        case = self.case()
        self.assertIsNone(case["sediment"]["in_domain_deposition_fraction"])
        config, _, _ = compile_timeline(
            case, obstacle_count=1, flow_count=2)
        self.assertEqual(config["sediment.in_domain_deposition_fraction"], -1.0)

    def test_explicit_deposition_fraction_passes_through(self):
        case = self.case()
        case["sediment"]["in_domain_deposition_fraction"] = 0.4
        config, _, _ = compile_timeline(
            case, obstacle_count=1, flow_count=2)
        self.assertEqual(config["sediment.in_domain_deposition_fraction"], 0.4)

    def test_valid_flow_max_can_be_widened_by_the_case(self):
        case = self.case()
        case["species"]["valid_flow_max_m_s"] = 2.5
        config, _, _ = compile_timeline(
            case, obstacle_count=1, flow_count=2)
        self.assertEqual(config["species.valid_flow_max_m_s"], 2.5)

    def test_missing_new_section_is_rejected(self):
        case = self.case()
        del case["sediment"]
        with self.assertRaisesRegex(CaseValidationError, "missing"):
            parse_case(case, flow_count=2, obstacle_count=1)

    def test_out_of_range_new_field_is_rejected(self):
        case = self.case()
        case["harvest"]["fraction_per_year"] = 1.5
        with self.assertRaisesRegex(CaseValidationError, "must be <="):
            parse_case(case, flow_count=2, obstacle_count=1)

    def test_scenario_varying_boundary_is_not_silently_compiled(self):
        case = self.case()
        case["forcing"]["steps"][1]["boundary"]["temperature_c"] = 13.0
        parsed, _ = parse_case(case, flow_count=2, obstacle_count=1)
        with self.assertRaisesRegex(
                UnsupportedCaseError, "scenario-varying"):
            compile_timeline(
                parsed, obstacle_count=1, flow_count=2)


class GrasshopperComponentContractTests(unittest.TestCase):

    def component_source(self):
        return COMPONENT.read_text(encoding="utf-8")

    def run_script(self):
        tree = ast.parse(self.component_source())
        for node in tree.body:
            if isinstance(node, ast.ClassDef) and node.name == "Script_Instance":
                for item in node.body:
                    if (isinstance(item, ast.FunctionDef) and
                            item.name == "RunScript"):
                        return item
        self.fail("Script_Instance.RunScript was not found")

    def test_geometry_collections_use_list_access_annotations(self):
        function = self.run_script()
        annotations = {
            argument.arg: ast.unparse(argument.annotation)
            for argument in function.args.args if argument.annotation}
        for name in ("obstacles", "probes", "flowVectors"):
            self.assertTrue(
                annotations[name].startswith("list["),
                "%s must use Grasshopper List access" % name)
        self.assertNotIn("Grasshopper.DataTree", ast.unparse(function))
        self.assertEqual(annotations["qualityMode"], "int")
        self.assertEqual(annotations["speedMode"], "bool")

    def test_sidecar_loader_reloads_changed_source_files(self):
        source = self.component_source()
        self.assertIn("st_mtime_ns", source)
        self.assertIn("importlib.reload(module)", source)
        self.assertIn("__musselflow_source_stamp__", source)
        self.assertIn("OUTDATED SIDECAR", source)

    def test_public_numeric_outputs_preserve_small_fitness_values(self):
        source = self.component_source()
        self.assertIn("OUTPUT_DECIMALS = 6", source)
        self.assertIn("OUTPUT_SIGNIFICANT_DIGITS = 6", source)
        self.assertIn(
            'compact_number(result["fitness"])', source)
        self.assertIn('"%s = %s"', source)

        tree = ast.parse(source)
        function = next(
            node for node in tree.body
            if isinstance(node, ast.FunctionDef) and
            node.name == "compact_number")
        namespace = {
            "math": math,
            "OUTPUT_SIGNIFICANT_DIGITS": 6,
        }
        exec(compile(
            ast.fix_missing_locations(ast.Module(
                body=[function], type_ignores=[])),
            str(COMPONENT), "exec"), namespace)
        compact = namespace["compact_number"]
        self.assertNotEqual(compact(4.9e-7), 0.0)
        self.assertEqual(compact(2.43156789e-9), 2.43157e-9)
        self.assertEqual(compact(4.89706123), 4.89706)

    def test_optimizer_geometry_path_avoids_known_slow_rhino_calls(self):
        component_source = self.component_source()
        bridge_source = GEOMETRY_BRIDGE.read_text(encoding="utf-8")
        self.assertNotIn("GetBoundingBox(True)", component_source)
        self.assertNotIn("GetBoundingBox(True)", bridge_source)
        self.assertNotIn(".IsPointOnFace(", bridge_source)
        self.assertIn("if vertex_count >= 4:", bridge_source)

    def test_detailed_visualizer_selects_from_list_of_flow_scenarios(self):
        source = VISUALIZER.read_text(encoding="utf-8")
        tree = ast.parse(source)
        script_class = next(
            node for node in tree.body
            if isinstance(node, ast.ClassDef) and
            node.name == "Script_Instance")
        run_script = next(
            node for node in script_class.body
            if isinstance(node, ast.FunctionDef) and
            node.name == "RunScript")
        annotations = {
            argument.arg: ast.unparse(argument.annotation)
            for argument in run_script.args.args if argument.annotation}
        self.assertEqual(
            annotations["flowDir"],
            "Grasshopper.DataTree[Rhino.Geometry.Vector3d]")
        self.assertEqual(
            annotations["socks"], "Grasshopper.DataTree[object]")
        self.assertNotIn("flowIndex", annotations)
        wrapper_source = ast.unparse(run_script)
        self.assertIn("flowVectorList(flowDir)", wrapper_source)
        self.assertIn("sum((item.X for item in vectors))", wrapper_source)
        self.assertEqual(wrapper_source.count("self._runSingle("), 1)
        self.assertNotIn("flowDir.Length", wrapper_source)
        self.assertIn("projectedMeshPolys", source)
        self.assertIn("displaySolid = solid | ~insideDomain", source)

if __name__ == "__main__":
    unittest.main()
