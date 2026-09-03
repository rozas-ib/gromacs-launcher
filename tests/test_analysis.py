import csv
import math
import os
import tempfile
import unittest

from launcher.analysis.analyses.molarity import compute_molarity_mol_l, run_molarity_analysis
from launcher.analysis.geometry import parse_gro_box_volume_nm3
from launcher.analysis.models import AnalysisCase
from launcher.analysis.reporting import write_csv


def write_empty_gro(path, box_line):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as gro_file:
        gro_file.write(f"test\n0\n{box_line}\n")


class GeometryTests(unittest.TestCase):
    def test_orthorhombic_volume(self):
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, "box.gro")
            write_empty_gro(path, "2 3 4")
            self.assertEqual(parse_gro_box_volume_nm3(path), 24.0)

    def test_triclinic_volume(self):
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, "box.gro")
            write_empty_gro(path, "2 3 4 0 0 0 0 0 0")
            self.assertEqual(parse_gro_box_volume_nm3(path), 24.0)


class MolarityTests(unittest.TestCase):
    def test_reports_each_replica_and_mean_of_valid_molarities(self):
        with tempfile.TemporaryDirectory() as root:
            write_empty_gro(os.path.join(root, "rep_1", "4_prod", "start.gro"), "5 5 5")
            write_empty_gro(os.path.join(root, "rep_2", "4_prod", "start.gro"), "10 10 10")
            case = AnalysisCase("case", root, 3, {"LiFSI_salt": 150})

            detail_headers, details, summary_headers, summaries = run_molarity_analysis(
                case, {"groups": ["LiFSI_salt"]}
            )

            self.assertIn("Molarity (mol/L)", detail_headers)
            self.assertIn("Mean Molarity (mol/L)", summary_headers)
            self.assertEqual(len(details), 3)
            self.assertEqual(details[2]["Status"], "Error")
            self.assertEqual(summaries[0]["Valid Replicas"], 2)
            self.assertEqual(summaries[0]["Status"], "Incomplete")
            expected = (
                compute_molarity_mol_l(150, 125) + compute_molarity_mol_l(150, 1000)
            ) / 2
            self.assertTrue(math.isclose(summaries[0]["Mean Molarity (mol/L)"], expected))

    def test_csv_always_contains_readable_headers(self):
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, "report.csv")
            headers = ["System", "Volume (nm^3)", "Molarity (mol/L)"]
            write_csv(path, headers, [])
            with open(path, newline="", encoding="utf-8") as csv_file:
                self.assertEqual(next(csv.reader(csv_file)), headers)


if __name__ == "__main__":
    unittest.main()
