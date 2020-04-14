import model

COORDINATE_LINE_COLUMNS = 4


class InputReader:
    def __init__(self, input_file_name):
        self._input_file_name = input_file_name
        self._xyz_data = []
        self._packmol_output_name = None

    @property
    def xyz_data(self):
        return self._xyz_data

    @property
    def packmol_output_name(self):
        return self._packmol_output_name

    def parse_packmol_input(self):
        input_file = open(self._input_file_name, 'r')

        for line in input_file:
            line = line.split()

            if line[0] == 'output':
                self._packmol_output_name = line[1]
            elif line[0] == 'structure':
                xyz_name = line[1]
                line = next(input_file).split()

                while line[0] != 'end':
                    if line[0] == 'number':
                        self._xyz_data.append((xyz_name, int(line[1])))

                    line = next(input_file).split()

        input_file.close()

    def read_xyz_data(self):
        molecules = []

        for (xyz_file_name, molecule_count) in self._xyz_data:
            xyz_file = open(xyz_file_name, 'r')
            dat_file = open(xyz_file_name.replace('.xyz', '.dat'), 'r')

            residue_name = dat_file.readline().replace('\n', '')
            atoms_number = int(xyz_file.readline())
            xyz_file.readline()  # omit commentary line in XYZ file
            atoms = []

            for xyz_line in xyz_file:
                xyz_line = xyz_line.split()
                if len(xyz_line) == COORDINATE_LINE_COLUMNS:
                    atom_symbol = xyz_line[0]
                    atom_coordinates = model.Coordinates(float(xyz_line[1]), float(xyz_line[2]), float(xyz_line[3]))

                    dat_line = dat_file.readline().split()
                    namd_symbol = dat_line[0]
                    charge = float(dat_line[1])
                    mass = float(dat_line[2])

                    atoms.append(model.Atom(atom_symbol, atom_coordinates, charge, mass, namd_symbol))

            if len(atoms) != atoms_number:
                print(
                    '[WARNING] Difference between declared and read atoms number in file {}, declared: {}, read: {}'.format(
                        xyz_file_name, atoms_number, len(atoms)))

            molecule = model.Molecule(atoms, residue_name)
            if 'BONDS' in dat_file.readline():
                bonds = []
                dat_line = dat_file.readline().split()
                while 'END' not in dat_line:
                    bond = (int(dat_line[0]), int(dat_line[1]))
                    bonds.append(bond)
                    dat_line = dat_file.readline().split()

                molecule.bonds = bonds
            else:
                molecule.determine_bonds()

            molecule.determine_angles()

            molecules.append((molecule, molecule_count))
            xyz_file.close()
            dat_file.close()

        system = model.System(molecules, self._packmol_output_name)
        return system
