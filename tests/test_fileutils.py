import fileutils
import unittest
import model

PACKMOL_FILE_NAME = 'na-f2ec-tfsi-1M-10.inp'


class TestFileUtils(unittest.TestCase):

    def test_input_reader(self):
        reader = fileutils.InputReader(PACKMOL_FILE_NAME)
        reader.parse_packmol_input()

        expected_output = [
            ('na.xyz', 47),
            ('f1ec.xyz', 461),
            ('tfsi.xyz', 47)
        ]
        self.assertEqual(expected_output, reader.xyz_data)

    def test_xyz_reader(self):
        reader = fileutils.InputReader(PACKMOL_FILE_NAME)
        reader.parse_packmol_input()

        xyz_data = reader.read_xyz_data()
        na_atom = model.Atom("Na", model.Coordinates(0.0, 0.0, 0.0))
        molecule = model.Molecule([na_atom])

        self.assertEqual(xyz_data[0], molecule)
