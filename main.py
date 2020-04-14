import sys
import readingutils
import savingutils

INPUT_ARGC_VALUE = 2
PSF_ARGC_VALUE = 3
PDB_ARGC_VALUE = 4


def main():
    if len(sys.argv) < INPUT_ARGC_VALUE:
        input_file_name = input("Enter input file name: ")
    else:
        input_file_name = sys.argv[1]

    reader = readingutils.InputReader(input_file_name)
    reader.parse_packmol_input()
    system = reader.read_xyz_data()

    if len(sys.argv) < PSF_ARGC_VALUE:
        psf_file_name = input("Enter psf file name: ")
    else:
        psf_file_name = sys.argv[2]
    psf_saver = savingutils.PSFSaver(psf_file_name, system)
    psf_saver.save_to_file()

    if len(sys.argv) < PDB_ARGC_VALUE:
        pdb_file_name = input("Enter pdb file name: ")
    else:
        pdb_file_name = sys.argv[3]
    pdb_saver = savingutils.PDBSaver(pdb_file_name, system)
    pdb_saver.save_to_file()


if __name__ == "__main__":
    main()
