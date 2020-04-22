import math
import unittest

from utils import model


class TestModel(unittest.TestCase):

    def test_compare_coordinates(self):
        coordinates1 = model.Coordinates(11.01, 12.02, 13.03)
        coordinates2 = model.Coordinates(11.01, 12.02, 13.03)
        coordinates3 = model.Coordinates(15.01, 14.02, 13.03)

        self.assertEqual(coordinates1, coordinates2)
        self.assertNotEqual(coordinates1, coordinates3)

    def test_norm(self):
        coordinates1 = model.Coordinates(1.0, 0.0, 0.0)
        coordinates2 = model.Coordinates(3.0, 4.0, 0.0)
        coordinates3 = model.Coordinates(-1.0, 0.0, 1.0)

        self.assertEqual(1.0, coordinates1.norm())
        self.assertEqual(5.0, coordinates2.norm())
        self.assertEqual(math.sqrt(2.0), coordinates3.norm())

    def test_distance(self):
        coordinates1 = model.Coordinates(1.0, 0.0, 0.0)
        coordinates2 = model.Coordinates(1.0, 0.0, 0.0)
        distance12 = coordinates1.calculate_distance(coordinates2)

        self.assertEqual(0.0, distance12)

        coordinates3 = model.Coordinates(4.0, 4.0, 0.0)
        distance13 = coordinates1.calculate_distance(coordinates3)
        self.assertEqual(5.0, distance13)

        coordinates4 = model.Coordinates(-3.0, 0.0, 4.0)
        distance14 = coordinates1.calculate_distance(coordinates4)
        distance43 = coordinates4.calculate_distance(coordinates3)
        self.assertGreaterEqual(distance14 + distance43, distance13)

    def test_compare_atoms(self):
        coordinates = model.Coordinates(12.03, 11.01, 10.02)
        atom1 = model.Atom("C", coordinates)
        atom2 = model.Atom("C", coordinates)
        atom3 = model.Atom("Ce", coordinates)

        self.assertEqual(atom1, atom2)
        self.assertNotEqual(atom1, atom3)

    def test_compare_molecules(self):
        coordinates1 = model.Coordinates()
        coordinates2 = model.Coordinates(1.0, 0.0, 0.0)
        coordinates3 = model.Coordinates(0.5, 1.0, 0.0)

        atom1 = model.Atom("H", coordinates1)
        atom2 = model.Atom("H", coordinates2)
        atom3 = model.Atom("H", coordinates3)

        molecule1 = model.Molecule([atom1, atom2], "MOL")
        molecule2 = model.Molecule([atom2, atom1], "MOL")
        molecule3 = model.Molecule([atom1, atom2, atom3], "MOL")

        self.assertEqual(molecule1, molecule2)
        self.assertNotEqual(molecule2, molecule3)

    def test_determine_bonds(self):
        coordinatesH1 = model.Coordinates(0.0, 0.0, 0.95)
        coordinatesH2 = model.Coordinates(0.89567, 0.00, -0.316663)
        coordinatesO1 = model.Coordinates(0.0, 0.0, 0.0)
        coordinatesO2 = model.Coordinates(4.0, 0.0, 0.0)

        atomH1 = model.Atom("H", coordinatesH1)
        atomH2 = model.Atom("H", coordinatesH2)
        atomO1 = model.Atom("O", coordinatesO1)
        atomO2 = model.Atom("O", coordinatesO2)

        molecule = model.Molecule([atomH1, atomH2, atomO1, atomO2], "MOL")
        molecule.determine_bonds(1.50)
        bonds = molecule.bonds
        expected_bonds = ((0, 2), (1, 2))
        self.assertEqual(bonds, expected_bonds)
