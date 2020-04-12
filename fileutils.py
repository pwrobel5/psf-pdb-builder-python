class InputReader:
    def __init__(self, input_file_name):
        self.input_file_name = input_file_name
        self.xyz_data = {}
        self.packmol_output_name = None

    def parse_input(self):
        input_file = open(self.input_file_name, 'r')

        for line in input_file:
            line = line.split()

            if line[0] == 'output':
                self.packmol_output_name = line[1]
            elif line[0] == 'structure':
                xyz_name = line[1]
                line = next(input_file).split()

                while line[0] != 'end':
                    if line[0] == 'number':
                        self.xyz_data[xyz_name] = int(line[1])

                    line = next(input_file).split()

        input_file.close()
