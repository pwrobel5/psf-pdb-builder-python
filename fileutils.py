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
