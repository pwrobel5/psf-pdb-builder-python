import sys
import readingutils
import savingutils

CORRECT_ARGC_VALUE = 2


def main():
    if len(sys.argv) != CORRECT_ARGC_VALUE:
        input_file_name = input("Enter input file name: ")
    else:
        input_file_name = sys.argv[1]

    reader = readingutils.InputReader(input_file_name)
    reader.parse_packmol_input()
    system = reader.read_xyz_data()

    psf_file_name = input("Enter psf file name: ")
    psf_saver = savingutils.PSFSaver(psf_file_name, system)
    psf_saver.save_to_file()


if __name__ == "__main__":
    main()
