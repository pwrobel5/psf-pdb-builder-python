import argparse

from utils import readingutils, savingutils


def main():
    parser = argparse.ArgumentParser(description="Create PSF and PDB file from Packmol/Tinker output")

    parser.add_argument("Input", metavar="input", type=str, help='Path to Packmol input file')
    parser.add_argument("-psf", type=str, help="Output PSF file name")
    parser.add_argument("-pdb", type=str, help="Output PDB file name")
    parser.add_argument("-t", "--tinker", action="store_true", help="Read .xyz in Tinker analyse format")

    args = parser.parse_args()
    input_file_name = args.Input

    reader = readingutils.InputReader(input_file_name, args.tinker)
    reader.parse_packmol_input()
    system = reader.read_xyz_data()

    psf_file_name = args.psf
    if psf_file_name:
        psf_saver = savingutils.PSFSaver(psf_file_name, system)
        psf_saver.save_to_file()

    pdb_file_name = args.pdb
    if pdb_file_name:
        pdb_saver = savingutils.PDBSaver(pdb_file_name, system)
        pdb_saver.save_to_file()

    print("Done")


if __name__ == "__main__":
    main()
