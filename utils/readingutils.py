import os

from utils import model


class InputReader:
    COORDINATE_LINE_COLUMNS = 4
    COORDINATE_LINE_COLUMNS_TINKER = 5
    DRUDE_LINE_COLUMNS = 5
    BOND_SECTION_BEGINNING = "Bond Stretching Parameters"
    ANGLE_SECTION_BEGINNING = "Angle Bending Parameters"
    DIHEDRAL_SECTION_BEGINNING = "Torsional Angle Parameters"
    CONN_LINES_TO_OMIT = 3

    def __init__(self, input_file_name, tinker_format=False, include_drude=False):
        self._input_file_name = input_file_name
        self._xyz_data = []
        self._packmol_output_name = None
        self._tinker_format = tinker_format
        self._include_drude = include_drude

    @property
    def xyz_data(self):
        return self._xyz_data

    @property
    def packmol_output_name(self):
        return self._packmol_output_name

    def parse_packmol_input(self):
        input_file = open(self._input_file_name, 'r')
        dirname = os.path.dirname(input_file.name)
        if dirname != "":
            dirname += "/"

        for line in input_file:
            line = line.split()

            if len(line) != 0:
                if line[0] == 'output':
                    self._packmol_output_name = dirname + line[1]
                elif line[0] == 'structure':
                    xyz_name = dirname + line[1]
                    line = next(input_file).split()

                    while line[0] != 'end':
                        if line[0] == 'number':
                            self._xyz_data.append((xyz_name, int(line[1])))

                        line = next(input_file).split()

        input_file.close()

    def read_xyz_data(self):
        if self._tinker_format:
            return self.__read_xyz_data_tinker()
        else:
            return self.__read_xyz_data_default()

    def __read_xyz_data_default(self):
        molecules = []

        for (xyz_file_name, molecule_count) in self._xyz_data:
            xyz_file = open(xyz_file_name, 'r')

            atoms_number = int(xyz_file.readline())
            xyz_file.readline()  # omit commentary line in XYZ file
            atoms = []

            for xyz_line in xyz_file:
                xyz_line = xyz_line.split()
                if len(xyz_line) >= self.COORDINATE_LINE_COLUMNS:
                    atom_symbol = xyz_line[0]
                    atom_coordinates = model.Coordinates(float(xyz_line[1]), float(xyz_line[2]), float(xyz_line[3]))

                    atoms.append(model.Atom(atom_symbol, atom_coordinates))

            xyz_file.close()

            if len(atoms) != atoms_number:
                print(
                    '[WARNING] Difference between declared and read atoms number in file {}, declared: {}, read: {}'.format(
                        xyz_file_name, atoms_number, len(atoms)))

            molecule = model.Molecule(atoms)
            self.__read_dat_data(xyz_file_name.replace('.xyz', '.dat'), molecule)

            molecule.determine_angles()
            molecule.determine_dihedrals()
            molecules.append((molecule, molecule_count))

        system = model.System(molecules, self._packmol_output_name)
        return system

    def __read_xyz_data_tinker(self):
        molecules = []

        for (xyz_file_name, molecule_count) in self._xyz_data:
            xyz_file = open(xyz_file_name, 'r')

            atoms_number = int(xyz_file.readline().split()[0])
            atoms = []

            for xyz_line in xyz_file:
                xyz_line = xyz_line.split()
                if len(xyz_line) >= self.COORDINATE_LINE_COLUMNS_TINKER:
                    atom_symbol = xyz_line[1]
                    atom_coordinates = model.Coordinates(float(xyz_line[2]), float(xyz_line[3]), float(xyz_line[4]))

                    atoms.append(model.Atom(atom_symbol, atom_coordinates))

            xyz_file.close()

            if len(atoms) != atoms_number:
                print(
                    '[WARNING] Difference between declared and read atoms number in file {}, declared: {}, read: {}'.format(
                        xyz_file_name, atoms_number, len(atoms)))

            molecule = model.Molecule(atoms)
            self.__read_dat_data(xyz_file_name.replace('.xyz', '.dat'), molecule)

            if molecule.atoms_number > 1:
                self.__read_conn_data(xyz_file_name.replace('.xyz', '.conn'), molecule)

            molecules.append((molecule, molecule_count))

        system = model.System(molecules, self._packmol_output_name)
        return system

    def __read_dat_data(self, dat_file_name, molecule):
        dat_file = open(dat_file_name, 'r')
        molecule.residue_name = dat_file.readline().replace('\n', '')
        drude_bonds_list = []

        for atom in molecule.atoms:
            dat_line = dat_file.readline().split()
            atom.namd_symbol = dat_line[0]
            atom.charge = float(dat_line[1])
            atom.mass = float(dat_line[2])

            if self._include_drude and len(dat_line) >= self.DRUDE_LINE_COLUMNS:
                molecule.add_drude_atom(atom, float(dat_line[4]), drude_bonds_list)

        molecule.drude_bonds = drude_bonds_list

        if not self._tinker_format:
            if 'BONDS' in dat_file.readline():
                shifts = molecule.shifts
                bonds = []
                dat_line = dat_file.readline().split()
                while 'END' not in dat_line:
                    bond = (self.shift(dat_line[0], shifts), self.shift(dat_line[1], shifts))
                    bonds.append(bond)
                    dat_line = dat_file.readline().split()

                molecule.bonds = bonds
            else:
                molecule.determine_bonds()

        dat_file.close()

    def __read_conn_data(self, conn_file_name, molecule):
        conn_file = open(conn_file_name, 'r')

        self.__jump_to_conn_section(conn_file, self.BOND_SECTION_BEGINNING, self.CONN_LINES_TO_OMIT)
        self.__read_bond_conn_data(conn_file, molecule)

        if molecule.atoms_number > 2:
            self.__jump_to_conn_section(conn_file, self.ANGLE_SECTION_BEGINNING, self.CONN_LINES_TO_OMIT)
            self.__read_angles_conn_data(conn_file, molecule)

        if molecule.atoms_number > 3:
            self.__jump_to_conn_section(conn_file, self.DIHEDRAL_SECTION_BEGINNING, self.CONN_LINES_TO_OMIT)
            self.__read_dihedrals_conn_data(conn_file, molecule)

        conn_file.close()

    @staticmethod
    def __jump_to_conn_section(conn_file, section_beginning, lines_to_omit):
        line = conn_file.readline()
        while section_beginning not in line:
            line = conn_file.readline()

        for _ in range(0, lines_to_omit):
            conn_file.readline()

    @staticmethod
    def __read_bond_conn_data(conn_file, molecule):
        shifts = molecule.shifts
        bonds = []
        line = conn_file.readline().strip()
        while line != "":
            line = line.split()
            bond = (InputReader.shift(line[1], shifts), InputReader.shift(line[2], shifts))
            bonds.append(bond)
            line = conn_file.readline().strip()

        molecule.bonds = bonds

    @staticmethod
    def __read_angles_conn_data(conn_file, molecule):
        shifts = molecule.shifts
        angles = []
        line = conn_file.readline().strip()
        while line != "":
            line = line.split()
            angle = (
                InputReader.shift(line[1], shifts),
                InputReader.shift(line[2], shifts),
                InputReader.shift(line[3], shifts)
            )
            angles.append(angle)
            line = conn_file.readline().strip()

        molecule.angles = angles

    @staticmethod
    def __read_dihedrals_conn_data(conn_file, molecule):
        shifts = molecule.shifts
        dihedrals = []
        line = conn_file.readline().strip()
        while line != "":
            line = line.split()
            dihedral = (
                InputReader.shift(line[1], shifts),
                InputReader.shift(line[2], shifts),
                InputReader.shift(line[3], shifts),
                InputReader.shift(line[4], shifts)
            )
            dihedrals.append(dihedral)
            line = conn_file.readline().strip()

        molecule.dihedrals = dihedrals

    @staticmethod
    def shift(string_value, shifts):
        result = int(string_value) - 1
        result += shifts[result]
        return result
