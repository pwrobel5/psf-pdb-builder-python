import readingutils
import unittest
import model

PACKMOL_FILE_NAME = 'li-ec-01.inp'


class TestFileUtils(unittest.TestCase):

    def test_input_reader(self):
        reader = readingutils.InputReader(PACKMOL_FILE_NAME)
        reader.parse_packmol_input()

        expected_output = [
            ('li.xyz', 4),
            ('ec.xyz', 46),
            ('tfsi.xyz', 4)
        ]
        self.assertEqual(expected_output, reader.xyz_data)

    def test_xyz_reader(self):
        reader = readingutils.InputReader(PACKMOL_FILE_NAME)
        reader.parse_packmol_input()

        system = reader.read_xyz_data()
        na_atom = model.Atom("Li", model.Coordinates(0.0, 0.0, 0.0))
        molecule = model.Molecule([na_atom])

        molecules = system.molecules
        self.assertEqual(molecules[0], molecule)
        self.assertEqual(system.xyz_file_name, 'li-ec-01.xyz')
