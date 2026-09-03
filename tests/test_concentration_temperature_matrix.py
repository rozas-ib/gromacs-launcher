import unittest

from launcher.case_matrix import build_master_combos, normalize_component_ratios
from launcher.concentration_optimizer import load_optimizer_configs, render_ratio_snippet


def optimizer_config(temperature):
    return {
        "concentration_optimizer": {
            "enabled": True,
            "target_group": "salt",
            "target_molarity_mol_l": [0.5, 1.0],
            "tolerance_mol_l": 0.05,
            "reference_count": 100,
            "box_size_nm": 10.0,
            "initial_ratio_name": "ratio_1",
            "force_field_name": "ff_1",
            "temperature": temperature,
        }
    }


class OptimizerTemperatureMatrixTests(unittest.TestCase):
    def test_scalar_temperature_remains_supported(self):
        configs = load_optimizer_configs(optimizer_config(298.15))
        self.assertEqual([(item.target_molarity_mol_l, item.temperature) for item in configs], [
            (0.5, 298.15),
            (1.0, 298.15),
        ])

    def test_molarity_temperature_cartesian_product(self):
        cfg = optimizer_config([298.15, 323.15])
        configs = load_optimizer_configs(cfg)
        self.assertEqual([(item.target_molarity_mol_l, item.temperature) for item in configs], [
            (0.5, 298.15),
            (0.5, 323.15),
            (1.0, 298.15),
            (1.0, 323.15),
        ])
        filtered = load_optimizer_configs(cfg, target_molarity_filter=1.0, temperature_filter=323.15)
        self.assertEqual(len(filtered), 1)
        self.assertEqual(filtered[0].temperature, 323.15)

    def test_generated_ratio_is_copy_pasteable_and_temperature_restricted(self):
        opt_cfg = load_optimizer_configs(optimizer_config([298.15]))[0]
        snippet = render_ratio_snippet({"salt": 1.0, "solvent": 4.25}, opt_cfg)
        self.assertIn('name = "ratio_1_0.5M_T298.15K"', snippet)
        self.assertIn("target_temps = [298.15]", snippet)
        self.assertIn("solvent = 4.25", snippet)


class LauncherTemperatureRestrictionTests(unittest.TestCase):
    def test_ratio_specific_temperatures_prevent_cross_combinations(self):
        cfg = {
            "screening": {
                "target_temps": [298.15, 323.15],
                "component_ratios": [
                    {"name": "low", "salt": 1, "target_temps": [298.15]},
                    {"name": "high", "salt": 2, "target_temps": [323.15]},
                ],
                "force_field": [{"name": "ff", "ion": "ion.itp"}],
            },
            "system_sizing": {"mode": "target_atoms", "target_atoms": 100},
        }
        combos = build_master_combos(cfg, ["salt"], ["ion"])
        self.assertEqual([(combo[0]["name"], combo[2]) for combo in combos], [
            ("low", 298.15),
            ("high", 323.15),
        ])

    def test_ratio_without_override_uses_global_temperatures(self):
        entries = normalize_component_ratios(
            {"component_ratios": [{"salt": 1}]}, ["salt"]
        )
        self.assertIsNone(entries[0]["target_temps"])


if __name__ == "__main__":
    unittest.main()
