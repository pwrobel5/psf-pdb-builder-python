import abc
import model


class FileSaver(metaclass=abc.ABCMeta):
    def __init__(self, output_file_name, system):
        self._output_file = open(output_file_name, 'w')
        self._system = system

    @abc.abstractmethod
    def save_to_file(self):
        raise NotImplementedError('Saving not implemented!')


class PSFSaver(FileSaver):
    PSF_HEADER = "PSF CMAP\n\n       1 !NTITLE\n REMARKS written by psf-pdb-builder\n\n"
    ATOM_LINE_FORMAT = "{:>8}    {:3} {:3d} {:3}  {:3}   {:3}  {:.6f}   {:.4f}  0\n"

    def save_to_file(self):
        self._output_file.write(self.PSF_HEADER)
        self.__save_atoms_section()
        self.__save_bonds_section()
        self._output_file.close()
        print('PSF file successfully written')

    def __save_atoms_section(self):
        atom_number = 1
        self._output_file.write("{:>8} !NATOM\n".format(self._system.atoms_number))

        for (molecule, counter) in self._system.molecules:
            residue_id = 1

            for i in range(0, counter):
                for atom in molecule.atoms:
                    self._output_file.write(self.ATOM_LINE_FORMAT.format(atom_number,
                                                                         self._system.segment_id,
                                                                         residue_id,
                                                                         molecule.residue_name,
                                                                         atom.symbol,
                                                                         atom.namd_symbol,
                                                                         atom.charge,
                                                                         atom.mass))
                    atom_number += 1
                residue_id += 1

    def __save_bonds_section(self):
        bonds = {}
        bonds_number = 0

        for (molecule, counter) in self._system.molecules:
            bonds[molecule] = molecule.bonds
            bonds_number += counter * len(bonds[molecule])

        self._output_file.write("{:>8} !NBOND: bonds\n".format(bonds_number))
        base = 1

        self._output_file.write(" ")
        items_in_line = 0

        for (molecule, counter) in self._system.molecules:
            bond_list = bonds[molecule]

            for i in range(0, counter):
                for (first, second) in bond_list:
                    self._output_file.write("{:7d} {:7d} ".format(first + base, second + base))
                    items_in_line = (items_in_line + 1) % 4

                    if items_in_line == 0:
                        self._output_file.write('\n ')
                base += molecule.atoms_number

        if items_in_line != 0:
            self._output_file.write('\n')


class PDBSaver(FileSaver):
    LINE_FORMAT = "ATOM {:6d}  {:3} {:3} {:5d}     {:7.3f} {:7.3f} {:7.3f} {:5.2f} {:5.2f}      {:4} {:3d}\n"

    def save_to_file(self):
        xyz_file = open(self._system.xyz_file_name, 'r')

        occupancy = 0.00
        temperature_factor = 0.00
        atom_number = 1

        # check if number of atoms in .xyz matches number given in system
        atoms_number = int(xyz_file.readline())
        if atoms_number != self._system.atoms_number:
            print('[ERROR] Number of atoms calculated from Packmol input does not match number of atoms in .xyz file')
            xyz_file.close()
            return

        # omit commentary line in .xyz
        xyz_file.readline()

        for (molecule, counter) in self._system.molecules:
            residue_id = 1

            for i in range(0, counter):
                for atom in molecule.atoms:
                    xyz_line = xyz_file.readline().split()

                    if xyz_line[0] != atom.symbol:
                        print('[ERROR] No match between atoms in line {}, symbol is {}, expected {}'.format(xyz_line,
                                                                                                            xyz_line[0],
                                                                                                            atom.symbol))
                        xyz_file.close()
                        self._output_file.close()
                        return

                    coordinates = model.Coordinates(float(xyz_line[1]), float(xyz_line[2]), float(xyz_line[3]))
                    self._output_file.write(self.LINE_FORMAT.format(atom_number,
                                                                    atom.namd_symbol,
                                                                    molecule.residue_name,
                                                                    residue_id,
                                                                    coordinates.x,
                                                                    coordinates.y,
                                                                    coordinates.z,
                                                                    occupancy,
                                                                    temperature_factor,
                                                                    self._system.segment_id,
                                                                    atom_number))

                    atom_number += 1
                residue_id += 1

        self._output_file.write('END\n')

        xyz_file.close()
        self._output_file.close()
        print('PDB file successfully written')
