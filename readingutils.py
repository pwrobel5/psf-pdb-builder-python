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

        for (xyz_file_name, _) in self._xyz_data:
            xyz_file = open(xyz_file_name)

            atoms_number = int(xyz_file.readline())
            xyz_file.readline()  # omit commentary line in XYZ file
            atoms = []

            for line in xyz_file:
                line = line.split()
                if len(line) == COORDINATE_LINE_COLUMNS:
                    atom_symbol = line[0]
                    atom_coordinates = model.Coordinates(float(line[1]), float(line[2]), float(line[3]))
                    atoms.append(model.Atom(atom_symbol, atom_coordinates))

            if len(atoms) != atoms_number:
                print(
                    '[WARNING] Difference between declared and read atoms number in file %s, declared: %d, read: %d' % (
                        xyz_file_name, atoms_number, len(atoms)))

            molecules.append(model.Molecule(atoms))
            xyz_file.close()

        segment = model.Segment(molecules)
        return segment
