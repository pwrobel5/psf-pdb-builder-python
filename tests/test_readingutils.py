import unittest

from utils import model, readingutils

PACKMOL_FILE_NAME = "li-ec/li-ec-01.inp"
PACKMOL_TINKER_FILE_NAME = "fsi-tinker/fsi.inp"
PACKMOL_DRUDE_FILE_NAME = "fsi-tinker-drude/fsi.inp"


class TestFileUtils(unittest.TestCase):

    def test_input_reader(self):
        reader = readingutils.InputReader(PACKMOL_FILE_NAME)
        reader.parse_packmol_input()

        expected_output = [
            ("li-ec/li.xyz", 4),
            ("li-ec/ec.xyz", 46),
            ("li-ec/tfsi.xyz", 4)
        ]
        self.assertEqual(expected_output, reader.xyz_data)

    def test_xyz_default_reader(self):
        reader = readingutils.InputReader(PACKMOL_FILE_NAME)
        reader.parse_packmol_input()

        system = reader.read_xyz_data()
        li_atom = model.Atom("Li", model.Coordinates(0.0, 0.0, 0.0), 1.00, 6.997, "kLi")
        molecule = model.Molecule([li_atom], "LI")

        molecules = system.molecules
        self.assertEqual(molecules[0][0], molecule)
        self.assertEqual(system.xyz_file_name, "li-ec/li-ec-01.xyz")

    def test_xyz_tinker_reader(self):
        reader = readingutils.InputReader(PACKMOL_TINKER_FILE_NAME, tinker_format=True)
        reader.parse_packmol_input()

        system = reader.read_xyz_data()
        f_atom = model.Atom("F", model.Coordinates(1.714650, -1.415751, 0.421888), -0.13, 18.9984, "Ff")

        molecule = system.molecules[0][0]
        atoms = molecule.atoms
        self.assertIn(f_atom, atoms)
        self.assertEqual(system.xyz_file_name, "fsi-tinker/fsi-res.xyz")

        bonds = molecule.bonds
        self.assertIn((0, 3), bonds)
        self.assertIn((5, 8), bonds)

        angles = molecule.angles
        self.assertIn((4, 5, 6), angles)
        self.assertIn((1, 0, 4), angles)

        dihedrals = molecule.dihedrals
        self.assertIn((3, 0, 4, 5), dihedrals)
        self.assertIn((1, 0, 4, 5), dihedrals)

    def test_xyz_reader_drude(self):
        reader = readingutils.InputReader(PACKMOL_DRUDE_FILE_NAME, tinker_format=True, include_drude=True)
        reader.parse_packmol_input()

        system = reader.read_xyz_data()
        molecule = system.molecules[0][0]
        atoms = molecule.atoms
        bonds = molecule.bonds
        drude_bonds = molecule.drude_bonds
        angles = molecule.angles
        dihedrals = molecule.dihedrals

        self.assertTrue(atoms[0].drude_atom is not None)

        self.assertIn((0, 1), drude_bonds)
        self.assertIn((0, 8), bonds)

        self.assertIn((4, 0, 6), angles)
        self.assertIn((8, 9, 11), angles)

        self.assertIn((4, 0, 8, 9), dihedrals)
        self.assertIn((0, 8, 9, 13), dihedrals)
