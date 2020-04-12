import fileutils
import unittest


class TestFileUtils(unittest.TestCase):

    def test_input_reader(self):
        reader = fileutils.InputReader('na-f2ec-tfsi-1M-10.inp')
        reader.parse_input()

        expected_output = {
            'na.xyz': 47,
            'f1ec.xyz': 461,
            'tfsi.xyz': 47
        }
        self.assertEqual(expected_output, reader.xyz_data)
