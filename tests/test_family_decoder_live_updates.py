import importlib.util
import pathlib
import sys
import types
import unittest


class _TypedList:
    @classmethod
    def __class_getitem__(cls, item):
        return list


class _Interval:
    def __init__(self, low, high):
        self.T0 = low
        self.T1 = high
        self.IsValid = True


def _load_decoder():
    system = types.ModuleType("System")
    system.Collections = types.SimpleNamespace(
        Generic=types.SimpleNamespace(List=_TypedList))
    rhino = types.ModuleType("Rhino")
    rhino.Geometry = types.SimpleNamespace(Interval=_Interval)
    grasshopper = types.ModuleType("Grasshopper")
    grasshopper.Kernel = types.SimpleNamespace(GH_ScriptInstance=object)
    sys.modules["System"] = system
    sys.modules["Rhino"] = rhino
    sys.modules["Grasshopper"] = grasshopper

    source = (
        pathlib.Path(__file__).parents[1]
        / "Miesmuschel"
        / "musselflow_family_decoder_gh_sdk.py")
    spec = importlib.util.spec_from_file_location("family_decoder", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DECODER = _load_decoder()


class FamilyDecoderLiveUpdateTests(unittest.TestCase):
    def setUp(self):
        self.component = DECODER.Script_Instance()
        self.u_genes = [0.31, 0.46, 0.22, 0.67, 0.38, 0.54, 0.29, 0.71]
        self.v_genes = [0.62, 0.27, 0.58, 0.41, 0.49, 0.36, 0.73, 0.18]
        self.rod_genes = [0.34, 0.66, 0.28, 0.72, 0.21, 0.43, 0.57, 0.39]
        self.domain = _Interval(0.8, 1.6)

    def evaluate(self, rods=None, u_genes=None):
        return self.component.RunScript(
            True,
            u_genes or self.u_genes,
            self.v_genes,
            rods or self.rod_genes,
            20,
            self.domain,
            False)

    def test_rod_randomization_changes_only_rod_outputs(self):
        first = self.evaluate()
        randomized = [0.81, 0.19, 0.77, 0.23, 0.69, 0.84, 0.35, 0.91]
        second = self.evaluate(rods=randomized)

        self.assertEqual(first[0], second[0])
        self.assertEqual(first[1], second[1])
        self.assertEqual(first[3], second[3])
        self.assertEqual(first[4], second[4])
        self.assertNotEqual(first[2], second[2])
        self.assertNotEqual(first[5][1], second[5][1])

    def test_every_rod_gene_has_a_continuous_effect(self):
        baseline = self.evaluate()[2]
        for index in range(8):
            changed = list(self.rod_genes)
            changed[index] = min(0.98, changed[index]+0.071)
            with self.subTest(gene=index):
                self.assertNotEqual(baseline, self.evaluate(rods=changed)[2])

    def test_rod_values_remain_inside_physical_domain(self):
        rods = self.evaluate()[2]
        self.assertEqual(40, len(rods))
        self.assertTrue(all(0.8 <= value <= 1.6 for value in rods))

    def test_u_gene_change_moves_positions_but_rod_values_stay_fixed(self):
        first = self.evaluate()
        changed_u = list(self.u_genes)
        changed_u[0] = 0.83
        second = self.evaluate(u_genes=changed_u)

        self.assertNotEqual(first[0], second[0])
        self.assertEqual(first[2], second[2])

    def test_report_identifies_new_live_build(self):
        report = self.evaluate()[5]
        self.assertIn("build 2026-08-25a", report[0])
        self.assertIn("LIVE INPUT", report[1])


if __name__ == "__main__":
    unittest.main()
