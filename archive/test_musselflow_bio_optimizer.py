"""Tests for the ecological grammar and MusselFlow biological screening core."""

import ast
import copy
from pathlib import Path
import time
import tempfile
import unittest

import numpy as np

from musselflow_bio_optimizer_core import (
    current_activity,
    evaluate_ensemble,
    evaluate_layout,
    feeding_state,
    oxygen_saturation_mg_l,
    particle_retention,
)
from musselflow_ecogrammar_core import (
    grammar_template,
    parse_grammar,
    resolved_lists,
)
from musselflow_surrogate_core import (
    LAYOUT_FEATURE_NAMES,
    extract_layout_features,
    fit_residual_ensemble,
    load_model,
    predict_residual,
    save_model,
)


def rectangle(width=40.0, height=30.0):
    return np.array([
        [0.0, 0.0],
        [width, 0.0],
        [width, height],
        [0.0, height],
    ])


def layout(count=12):
    x = np.linspace(8.0, 32.0, count)
    y = 15.0+4.0*np.sin(np.linspace(0.0, 3.0*np.pi, count))
    return np.column_stack((
        x,
        y,
        np.full(count, 2.0),
        np.full(count, 0.5),
        np.linspace(0.0, np.pi, count),
        np.full(count, 6.0),
    ))


def config_for(obstacle_count=12, scenario_count=2, extra=()):
    config, warnings, errors, notes = parse_grammar(
        list(extra), obstacle_count, scenario_count)
    if errors:
        raise AssertionError(errors)
    return config


class GrammarTests(unittest.TestCase):

    def test_template_parses_and_lists_resolve(self):
        config, warnings, errors, notes = parse_grammar(
            grammar_template(), obstacle_count=20, scenario_count=2)
        self.assertEqual(errors, [])
        resolved = resolved_lists(config, 20, 2)
        self.assertEqual(len(resolved["structure.porosity"]), 20)
        self.assertEqual(resolved["scenario.weights"], [0.5, 0.5])

    def test_template_also_accepts_one_flow_scenario(self):
        config, warnings, errors, notes = parse_grammar(
            grammar_template(), obstacle_count=20, scenario_count=1)
        self.assertEqual(errors, [])
        resolved = resolved_lists(config, 20, 1)
        self.assertEqual(resolved["scenario.weights"], [1.0])
        self.assertEqual(resolved["scenario.duration_h"], [6.0])

    def test_unknown_computational_key_is_rejected_but_notes_are_inert(self):
        config, warnings, errors, notes = parse_grammar([
            "site.dept_m = 10",
            "note.intent = maximize flushing",
            "source.local = field campaign 2026",
        ])
        self.assertTrue(any("unknown key" in error for error in errors))
        self.assertEqual(notes["note.intent"], "maximize flushing")
        self.assertEqual(config["site.depth_m"], 12.0)

    def test_bad_per_obstacle_and_scenario_lengths_are_rejected(self):
        config, warnings, errors, notes = parse_grammar([
            "structure.porosity = 0.5, 0.6",
            "scenario.weights = 1, 2, 3",
        ], obstacle_count=5, scenario_count=2)
        self.assertTrue(any("structure.porosity" in error for error in errors))
        self.assertTrue(any("scenario.weights" in error for error in errors))

    def test_dangerous_rates_are_rejected_and_subzero_water_is_parseable(self):
        config, warnings, errors, notes = parse_grammar([
            "site.temperature_c = -1",
            "species.clearance_q10 = 0",
            "hydrodynamics.min_speed_ratio = 1.2",
            "stocking.dry_tissue_kg_per_obstacle = -0.5",
        ])
        self.assertEqual(config["site.temperature_c"], -1.0)
        self.assertTrue(any("species.clearance_q10" in error
                            for error in errors))
        self.assertTrue(any("hydrodynamics.min_speed_ratio" in error
                            for error in errors))
        self.assertTrue(any("dry_tissue" in error for error in errors))


class GrasshopperSdkContractTests(unittest.TestCase):

    def test_sdk_signature_is_typed_and_return_count_is_stable(self):
        path = Path(__file__).with_name(
            "musselflow_bio_optimizer_gh_sdk.py")
        source = path.read_text(encoding="utf-8")
        self.assertNotIn("import rhinoscriptsyntax", source)
        self.assertIn("FLOW VECTOR ERROR", source)
        self.assertIn("domain plane, so plan-view speed is zero", source)
        tree = ast.parse(source, filename=str(path))
        run_script = next(
            node for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef) and node.name == "RunScript")
        self.assertEqual(
            [argument.arg for argument in run_script.args.args],
            ["self", "run", "obstacles", "domain", "probes",
             "flowVectors", "grammar"])
        self.assertTrue(all(
            argument.annotation is not None
            for argument in run_script.args.args[1:]))
        collection_annotations = {
            argument.arg: ast.unparse(argument.annotation)
            for argument in run_script.args.args
            if argument.arg in ("obstacles", "probes", "flowVectors", "grammar")
        }
        self.assertTrue(all(
            "Grasshopper.DataTree" in annotation
            for annotation in collection_annotations.values()))
        tuple_returns = [
            node for node in ast.walk(run_script)
            if isinstance(node, ast.Return) and
            isinstance(node.value, ast.Tuple)]
        self.assertEqual([len(node.value.elts) for node in tuple_returns], [33])

        empty_function = next(
            node for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef) and
            node.name == "empty_outputs")
        empty_returns = [
            node for node in ast.walk(empty_function)
            if isinstance(node, ast.Return) and
            isinstance(node.value, ast.Tuple)]
        self.assertEqual([len(node.value.elts) for node in empty_returns], [33])


class PhysicalContractTests(unittest.TestCase):

    def test_particle_size_changes_retention(self):
        config = config_for(1, 1)
        self.assertLess(
            float(particle_retention(2.0, config)),
            float(particle_retention(15.0, config)))

    def test_high_food_saturates_clearance_and_reduces_assimilation(self):
        config = config_for(1, 1)
        moderate = feeding_state(5.0, 3.0, 0.15, config)
        high = feeding_state(50.0, 3.0, 0.15, config)
        self.assertLess(high["clearance_activity"],
                        moderate["clearance_activity"])
        self.assertLess(high["assimilation_efficiency"],
                        moderate["assimilation_efficiency"])

    def test_low_oxygen_reduces_feeding_activity(self):
        normal = config_for(1, 1, ["site.boundary_do_mg_l = 9"])
        low = config_for(1, 1, ["site.boundary_do_mg_l = 1"])
        self.assertLess(
            feeding_state(5.0, 3.0, 0.15, low)["oxygen_activity"],
            feeding_state(5.0, 3.0, 0.15, normal)["oxygen_activity"])

    def test_aggregation_protects_clearance_in_current(self):
        small = config_for(1, 1, [
            "stocking.effective_aggregation_size = 3"])
        protected = config_for(1, 1, [
            "stocking.effective_aggregation_size = 20"])
        self.assertLess(current_activity(0.4, small),
                        current_activity(0.4, protected))

    def test_high_tsm_increases_pseudofaeces(self):
        config = config_for(1, 1)
        low = feeding_state(5.0, 1.0, 0.15, config)
        high = feeding_state(5.0, 20.0, 0.15, config)
        self.assertLess(low["pseudofaeces_fraction"],
                        high["pseudofaeces_fraction"])

    def test_organic_mass_partition_closes(self):
        obstacles = layout()
        config = config_for(len(obstacles), 1)
        result = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.15, 0.0]]), config)
        partition = (
            result["assimilated_organic_kg_day"]+
            result["faeces_organic_kg_day"]+
            result["pseudofaeces_organic_kg_day"])
        self.assertAlmostEqual(
            result["particulate_capture_kg_day"] *
            config["site.particulate_organic_fraction"],
            partition, places=12)
        self.assertTrue(np.isfinite(result["scope_for_growth_kj_day"]))
        self.assertGreaterEqual(result["potential_growth_g_dw_day"], 0.0)

    def test_oxygen_solubility_reference_values(self):
        self.assertAlmostEqual(
            oxygen_saturation_mg_l(20.0, 35.0), 7.3950, places=3)
        self.assertAlmostEqual(
            oxygen_saturation_mg_l(10.0, 20.0), 9.934, places=3)

    def test_zero_current_has_no_advective_capture(self):
        obstacles = layout()
        config = config_for(len(obstacles), 1)
        result = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.0, 0.0]]), config)
        self.assertEqual(result["chlorophyll_capture_g_day"], 0.0)
        self.assertEqual(result["particulate_capture_kg_day"], 0.0)
        self.assertGreater(result["mussel_respiration_kg_o2_day"], 0.0)
        self.assertEqual(result["probe_source"], "obstacle_centres")

    def test_capture_never_exceeds_domain_inflow(self):
        obstacles = layout(80)
        # Co-locate sources in streamwise space so this deliberately impossible
        # stress case tries to process more water than crosses the domain.
        obstacles[:, 0] = 20.0
        obstacles[:, 2] = 10.0
        obstacles[:, 3] = 10.0
        obstacles[:, 5] = 12.0
        config = config_for(80, 1, [
            "stocking.mussels_per_obstacle = 1000000",
            "structure.porosity = 0.95",
        ])
        speed = 0.01
        chlorophyll = config["site.chlorophyll_ug_l"]
        domain_inflow_g_day = (
            speed*30.0*config["site.depth_m"]*chlorophyll *
            86400.0/1000.0)
        result = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[speed, 0.0]]), config)
        self.assertLessEqual(
            result["chlorophyll_capture_g_day"],
            domain_inflow_g_day*(1.0+1e-12))
        self.assertLess(result["scenario_mass_balance_scale"][0], 1.0)

    def test_zero_chlorophyll_does_not_erase_other_particulates(self):
        obstacles = layout()
        config = config_for(len(obstacles), 1, [
            "site.chlorophyll_ug_l = 0",
            "site.tsm_mg_l = 4",
        ])
        result = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.15, 0.0]]), config)
        self.assertEqual(result["chlorophyll_capture_g_day"], 0.0)
        self.assertGreater(result["particulate_capture_kg_day"], 0.0)

    def test_salinity_profile_changes_activity_but_is_explicitly_editable(self):
        obstacles = layout()
        normal = config_for(len(obstacles), 1, [
            "site.salinity_psu = 20",
        ])
        closed = config_for(len(obstacles), 1, [
            "site.salinity_psu = 2",
        ])
        normal_result = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.15, 0.0]]), normal)
        closed_result = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.15, 0.0]]), closed)
        self.assertEqual(closed_result["salinity_activity"], 0.0)
        self.assertEqual(closed_result["chlorophyll_capture_g_day"], 0.0)
        self.assertGreater(
            normal_result["chlorophyll_capture_g_day"], 0.0)

    def test_more_solidity_causes_a_stronger_wake(self):
        obstacles = layout()
        probes = np.array([[35.0, 15.0]])
        open_config = config_for(len(obstacles), 1, [
            "structure.porosity = 0.85",
        ])
        dense_config = config_for(len(obstacles), 1, [
            "structure.porosity = 0.30",
        ])
        open_result = evaluate_layout(
            obstacles, rectangle(), probes, np.array([[0.2, 0.0]]),
            open_config)
        dense_result = evaluate_layout(
            obstacles, rectangle(), probes, np.array([[0.2, 0.0]]),
            dense_config)
        self.assertLess(
            dense_result["probe_speed_ratio"][0],
            open_result["probe_speed_ratio"][0])

    def test_flow_vector_length_changes_absolute_transport(self):
        obstacles = layout()
        config = config_for(len(obstacles), 1)
        slow = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.02, 0.0]]), config)
        fast = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.2, 0.0]]), config)
        self.assertGreater(
            np.mean(fast["probe_speed_m_s"]),
            np.mean(slow["probe_speed_m_s"]))
        self.assertNotEqual(
            fast["oxygen"]["final_mg_l"],
            slow["oxygen"]["final_mg_l"])

    def test_rotating_anisotropic_obstacles_changes_flow(self):
        obstacles = layout()
        obstacles[:, 2] = 4.0
        obstacles[:, 3] = 0.25
        probes = np.array([[35.0, 15.0], [35.0, 13.0], [35.0, 17.0]])
        config = config_for(len(obstacles), 1)
        aligned = obstacles.copy()
        aligned[:, 4] = 0.0
        normal = obstacles.copy()
        normal[:, 4] = 0.5*np.pi
        aligned_result = evaluate_layout(
            aligned, rectangle(), probes, np.array([[0.15, 0.0]]), config)
        normal_result = evaluate_layout(
            normal, rectangle(), probes, np.array([[0.15, 0.0]]), config)
        self.assertLess(
            np.mean(normal_result["probe_speed_ratio"]),
            np.mean(aligned_result["probe_speed_ratio"]))

    def test_refined_frontal_fill_changes_wake_and_is_scenario_specific(self):
        obstacles = layout(4)
        probes = np.array([[35.0, 15.0]])
        flows = np.array([[0.15, 0.0], [0.0, 0.15]])
        config = config_for(len(obstacles), 2)
        baseline = evaluate_layout(
            obstacles, rectangle(), probes, flows, config)
        low_fill_profiles = [
            {
                "frontal_area_m2": np.full(len(obstacles), 0.2),
                "frontal_fill": np.full(len(obstacles), 0.05),
            },
            {
                "frontal_area_m2": np.full(len(obstacles), 0.4),
                "frontal_fill": np.full(len(obstacles), 0.10),
            },
        ]
        refined = evaluate_layout(
            obstacles, rectangle(), probes, flows, config,
            hydraulic_profiles=low_fill_profiles)
        self.assertGreater(
            np.mean(refined["probe_speed_ratio"]),
            np.mean(baseline["probe_speed_ratio"]))
        np.testing.assert_allclose(
            refined["scenario_frontal_fill"][0], 0.05)
        np.testing.assert_allclose(
            refined["scenario_frontal_fill"][1], 0.10)

    def test_refined_profiles_must_match_flow_count(self):
        obstacles = layout(4)
        config = config_for(len(obstacles), 2)
        with self.assertRaisesRegex(
                ValueError, "one profile per flow vector"):
            evaluate_layout(
                obstacles, rectangle(), np.empty((0, 2)),
                np.array([[0.15, 0.0], [0.0, 0.15]]), config,
                hydraulic_profiles=[None])

    def test_hard_constraint_dominates_objective(self):
        obstacles = layout()
        unconstrained = config_for(len(obstacles), 1)
        impossible = config_for(len(obstacles), 1, [
            "constraint.min_do_mg_l = 20",
        ])
        feasible_result = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.15, 0.0]]), unconstrained)
        infeasible_result = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.15, 0.0]]), impossible)
        self.assertTrue(feasible_result["feasible"])
        self.assertFalse(infeasible_result["feasible"])
        self.assertGreaterEqual(feasible_result["fitness"], 0.25)
        self.assertLessEqual(infeasible_result["fitness"], 0.25)

    def test_declared_use_envelope_is_a_constraint_unless_overridden(self):
        obstacles = layout()
        guarded = config_for(len(obstacles), 1, [
            "species.valid_flow_max_m_s = 0.5",
            "validation.allow_extrapolation = false",
        ])
        explicit_override = config_for(len(obstacles), 1, [
            "species.valid_flow_max_m_s = 0.5",
            "validation.allow_extrapolation = true",
        ])
        guarded_result = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.8, 0.0]]), guarded)
        override_result = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)),
            np.array([[0.8, 0.0]]), explicit_override)
        self.assertFalse(guarded_result["feasible"])
        self.assertIn(
            "validation_envelope", guarded_result["constraint_margins"])
        self.assertTrue(override_result["feasible"])

    def test_deterministic_and_finite_under_random_stress(self):
        rng = np.random.default_rng(20260728)
        for _ in range(20):
            count = int(rng.integers(1, 80))
            obstacles = np.column_stack((
                rng.uniform(1.0, 39.0, count),
                rng.uniform(1.0, 29.0, count),
                rng.uniform(0.2, 4.0, count),
                rng.uniform(0.1, 2.0, count),
                rng.uniform(-np.pi, np.pi, count),
                rng.uniform(0.5, 10.0, count),
            ))
            probes = rng.uniform([0.0, 0.0], [40.0, 30.0], (40, 2))
            flows = rng.uniform(-0.3, 0.3, (3, 2))
            config = config_for(count, 3)
            first = evaluate_layout(
                obstacles, rectangle(), probes, flows, config)
            second = evaluate_layout(
                obstacles, rectangle(), probes, flows, config)
            self.assertEqual(first["fitness"], second["fitness"])
            self.assertTrue(np.isfinite(first["fitness"]))
            for key in (
                    "probe_speed_m_s", "probe_food_fraction",
                    "obstacle_speed_m_s",
                    "chlorophyll_capture_g_day_by_obstacle"):
                self.assertTrue(np.all(np.isfinite(first[key])), key)
                np.testing.assert_array_equal(first[key], second[key])


class PerformanceCharacterizationTests(unittest.TestCase):

    def test_typical_optimizer_case_remains_reduced_order(self):
        """Guard against accidentally reintroducing a seconds-long grid solve."""
        rng = np.random.default_rng(8)
        count = 100
        obstacles = np.column_stack((
            rng.uniform(1.0, 99.0, count),
            rng.uniform(1.0, 59.0, count),
            rng.uniform(0.5, 3.0, count),
            rng.uniform(0.2, 1.0, count),
            rng.uniform(-np.pi, np.pi, count),
            rng.uniform(2.0, 10.0, count),
        ))
        probes = rng.uniform([0.0, 0.0], [100.0, 60.0], (400, 2))
        flows = np.array([
            [0.15, 0.00],
            [-0.10, 0.00],
            [0.00, 0.08],
            [0.00, -0.08],
        ])
        config = config_for(count, len(flows))
        domain = rectangle(100.0, 60.0)

        # Warm imports/BLAS setup before timing.
        evaluate_layout(obstacles, domain, probes, flows, config)
        elapsed = []
        for unused_index in range(5):
            start = time.perf_counter()
            evaluate_layout(obstacles, domain, probes, flows, config)
            elapsed.append((time.perf_counter()-start)*1000.0)
        # A median rejects desktop scheduling spikes without hiding a
        # systematic regression back to a slow grid solver.
        self.assertLess(float(np.median(elapsed)), 300.0)


class ProbePhysicsDecouplingTests(unittest.TestCase):
    """Diagnostic probes must not change farm physics or farm-wide objectives.

    Architecture rule 3.4 / build Phase 4: adding, moving, or removing user
    probes may only change probe diagnostic outputs (and any explicitly
    probe-defined constraint), never obstacle capture, biodeposition,
    dissolved oxygen, or the farm objective/fitness values.
    """

    def farm_signature(self, result):
        oxygen = result["oxygen"]
        return {
            "fitness": result["fitness"],
            "raw_objective": result["raw_objective"],
            "feasible": result["feasible"],
            "objectives": dict(result["objectives"]),
            "minimum_do": oxygen["minimum_mg_l"],
            "final_do": oxygen["final_mg_l"],
            "deposition": oxygen["weighted_deposition_kg_m2_day"],
            "chlorophyll_capture": result["chlorophyll_capture_g_day"],
            "particulate_capture": result["particulate_capture_kg_day"],
            "cleared": result["effective_cleared_m3_day"],
        }

    def assert_same_farm(self, expected, actual):
        self.assertEqual(expected["feasible"], actual["feasible"])
        for key in ("fitness", "raw_objective", "minimum_do", "final_do",
                    "deposition", "chlorophyll_capture",
                    "particulate_capture", "cleared"):
            self.assertEqual(expected[key], actual[key], key)
        for name, value in expected["objectives"].items():
            self.assertEqual(value, actual["objectives"][name], name)

    def test_probes_do_not_change_farm_physics_or_objectives(self):
        obstacles = layout()
        domain = rectangle()
        flows = np.array([[0.15, 0.0], [-0.10, 0.0]])
        config = config_for(len(obstacles), len(flows))
        rng = np.random.default_rng(7)

        no_probes = evaluate_layout(
            obstacles, domain, np.empty((0, 2)), flows, config)
        few = evaluate_layout(
            obstacles, domain,
            np.array([[35.0, 15.0], [10.0, 20.0], [20.0, 5.0]]),
            flows, config)
        many = evaluate_layout(
            obstacles, domain,
            rng.uniform([0.0, 0.0], [40.0, 30.0], (60, 2)), flows, config)

        base = self.farm_signature(no_probes)
        self.assert_same_farm(base, self.farm_signature(few))
        self.assert_same_farm(base, self.farm_signature(many))

        # The obstacle-centre control set that physics uses is probe-invariant.
        np.testing.assert_array_equal(
            no_probes["obstacle_speed_ratio"], few["obstacle_speed_ratio"])
        np.testing.assert_array_equal(
            no_probes["obstacle_speed_ratio"], many["obstacle_speed_ratio"])

    def test_probe_diagnostics_still_respond_to_probe_placement(self):
        # Decoupling must not silence the probes themselves: a probe in the
        # wake must still read a different speed than one in the free stream.
        obstacles = layout()
        domain = rectangle()
        flows = np.array([[0.15, 0.0]])
        config = config_for(len(obstacles), len(flows))
        upstream = evaluate_layout(
            obstacles, domain, np.array([[2.0, 15.0]]), flows, config)
        wake = evaluate_layout(
            obstacles, domain, np.array([[38.0, 15.0]]), flows, config)
        self.assertNotAlmostEqual(
            float(upstream["probe_speed_m_s"][0]),
            float(wake["probe_speed_m_s"][0]))


class EnsembleEvaluatorTests(unittest.TestCase):
    """Native ensemble forcing: independent states combined by probability,
    with worst-case hard constraints and order-independent aggregates.
    """

    def setup(self):
        obstacles = layout()
        domain = rectangle()
        flows = np.array([[0.15, 0.0], [-0.10, 0.0]])
        config = config_for(len(obstacles), 1)
        probes = np.empty((0, 2))
        return obstacles, domain, flows, config, probes

    def states(self, config, probabilities=(0.6, 0.4)):
        return [
            {"id": "flood", "flow_vector_index": 0,
             "probability": probabilities[0], "config": copy.deepcopy(config)},
            {"id": "ebb", "flow_vector_index": 1,
             "probability": probabilities[1], "config": copy.deepcopy(config)},
        ]

    def test_objectives_are_probability_weighted(self):
        obstacles, domain, flows, config, probes = self.setup()
        flood = evaluate_layout(obstacles, domain, probes, flows[0:1], config)
        ebb = evaluate_layout(obstacles, domain, probes, flows[1:2], config)
        result = evaluate_ensemble(
            obstacles, domain, probes, flows, self.states(config, (0.6, 0.4)))
        for name in flood["objectives"]:
            expected = (
                0.6*flood["objectives"][name] + 0.4*ebb["objectives"][name])
            self.assertAlmostEqual(result["objectives"][name], expected)
        self.assertAlmostEqual(
            result["chlorophyll_capture_g_day"],
            0.6*flood["chlorophyll_capture_g_day"] +
            0.4*ebb["chlorophyll_capture_g_day"])

    def test_aggregate_is_independent_of_state_order(self):
        obstacles, domain, flows, config, probes = self.setup()
        forward = evaluate_ensemble(
            obstacles, domain, probes, flows, self.states(config, (0.6, 0.4)))
        reversed_states = list(reversed(self.states(config, (0.6, 0.4))))
        backward = evaluate_ensemble(
            obstacles, domain, probes, flows, reversed_states)
        self.assertEqual(forward["fitness"], backward["fitness"])
        self.assertEqual(forward["feasible"], backward["feasible"])
        self.assertEqual(forward["minimum_do_mg_l"], backward["minimum_do_mg_l"])
        for name in forward["objectives"]:
            self.assertAlmostEqual(
                forward["objectives"][name], backward["objectives"][name])

    def test_unnormalised_probabilities_are_renormalised(self):
        obstacles, domain, flows, config, probes = self.setup()
        normalised = evaluate_ensemble(
            obstacles, domain, probes, flows, self.states(config, (0.6, 0.4)))
        scaled = evaluate_ensemble(
            obstacles, domain, probes, flows, self.states(config, (6.0, 4.0)))
        self.assertAlmostEqual(normalised["fitness"], scaled["fitness"])

    def test_worst_state_makes_the_ensemble_infeasible(self):
        obstacles, domain, flows, config, probes = self.setup()
        states = self.states(config, (0.6, 0.4))
        # Force the ebb state outside its declared current envelope so its
        # validation_envelope margin goes negative; the ensemble must inherit
        # that worst-case infeasibility even though the flood state is fine.
        states[1]["config"]["species.valid_flow_max_m_s"] = 0.05
        result = evaluate_ensemble(obstacles, domain, probes, flows, states)
        self.assertFalse(result["feasible"])
        self.assertEqual(result["status"], "OUT_OF_ENVELOPE")
        self.assertLess(
            result["constraint_margins"]["validation_envelope"], 0.0)


class GuardedSurrogateTests(unittest.TestCase):

    def test_layout_feature_schema_is_finite_and_stable(self):
        obstacles = layout()
        flows = np.array([[0.15, 0.0], [-0.10, 0.0]])
        config = config_for(len(obstacles), len(flows))
        baseline = evaluate_layout(
            obstacles, rectangle(), np.empty((0, 2)), flows, config)
        row, names = extract_layout_features(
            obstacles, rectangle(), flows, config, baseline)
        self.assertEqual(names, LAYOUT_FEATURE_NAMES)
        self.assertEqual(len(row), len(names))
        self.assertTrue(np.all(np.isfinite(row)))

    def test_residual_model_round_trip_and_domain_guard(self):
        rng = np.random.default_rng(19)
        features = rng.uniform(-1.0, 1.0, (120, 3))
        residuals = np.column_stack((
            0.2*features[:, 0]-0.1*features[:, 1]**2,
            features[:, 0]*features[:, 2]+0.05*features[:, 1],
        ))
        model = fit_residual_ensemble(
            features, residuals,
            ["solidity", "reynolds", "blockage"],
            ["capture_residual", "flow_residual"],
            ensemble_size=6)
        prediction = predict_residual(features[:5], model)
        np.testing.assert_allclose(
            prediction["mean"], residuals[:5], atol=0.08)
        self.assertTrue(np.all(prediction["in_domain"]))

        outside = predict_residual(np.array([[20.0, 0.0, 0.0]]), model)
        self.assertFalse(outside["in_domain"][0])

        with tempfile.TemporaryDirectory() as folder:
            path = folder+"/surrogate.npz"
            save_model(path, model)
            loaded = load_model(path)
            reloaded = predict_residual(features[:5], loaded)
        np.testing.assert_allclose(
            prediction["mean"], reloaded["mean"], rtol=0.0, atol=0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
