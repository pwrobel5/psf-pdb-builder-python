import model
import unittest


class TestModel(unittest.TestCase):

    def test_compare_coordinates(self):
        coordinates1 = model.Coordinates(11.01, 12.02, 13.03)
        coordinates2 = model.Coordinates(11.01, 12.02, 13.03)
        coordinates3 = model.Coordinates(15.01, 14.02, 13.03)

        self.assertEqual(coordinates1, coordinates2)
        self.assertNotEqual(coordinates1, coordinates3)

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
