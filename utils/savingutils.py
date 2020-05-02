import abc

from utils import model


class FileSaver(metaclass=abc.ABCMeta):
    def __init__(self, output_file_name, system):
        self._output_file = open(output_file_name, 'w')
        self._system = system

    @abc.abstractmethod
    def save_to_file(self):
        raise NotImplementedError('Saving not implemented!')


class PSFSaver(FileSaver):
    PSF_HEADER = "PSF CMAP\n\n       1 !NTITLE\n REMARKS written by psf-pdb-builder\n\n"
    ATOM_LINE_FORMAT = "{:>8}    {:3} {:3d} {:3}  {:3}   {:3}  {:.6f}   {:.4f}  0   {:.5f}       {:.5f}\n"

    def save_to_file(self):
        self._output_file.write(self.PSF_HEADER)
        self.__save_atoms_section()
        self.__save_bonds_section()
        self.__save_angles_section()
        self.__save_dihedrals_section()
        self.__save_file_ending()
        self._output_file.close()
        print('PSF file successfully written')

    def __save_atoms_section(self):
        atom_number = 1
        self._output_file.write("{:>8} !NATOM\n".format(self._system.atoms_number_with_drude))

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
                                                                         atom.mass,
                                                                         atom.polarizability,
                                                                         atom.last_parameter))

                    drude_atom = atom.drude_atom
                    if drude_atom is not None:
                        atom_number += 1
                        self._output_file.write(self.ATOM_LINE_FORMAT.format(atom_number,
                                                                             self._system.segment_id,
                                                                             residue_id,
                                                                             molecule.residue_name,
                                                                             drude_atom.symbol,
                                                                             drude_atom.namd_symbol,
                                                                             drude_atom.charge,
                                                                             drude_atom.mass,
                                                                             drude_atom.polarizability,
                                                                             drude_atom.last_parameter))

                    atom_number += 1
                residue_id += 1

    def __save_bonds_section(self):
        bonds = {}
        bonds_number = 0

        for (molecule, counter) in self._system.molecules:
            bonds[molecule] = molecule.bonds + molecule.drude_bonds
            bonds_number += counter * len(bonds[molecule])

        self._output_file.write("{:>8} !NBOND: bonds\n".format(bonds_number))
        base = 1

        self._output_file.write(" ")
        items_in_line = 0

        bond_counter = 0

        for (molecule, counter) in self._system.molecules:
            bond_list = bonds[molecule]

            for i in range(0, counter):
                for (first, second) in bond_list:
                    self._output_file.write("{:7d} {:7d} ".format(first + base, second + base))
                    items_in_line = (items_in_line + 1) % 4
                    bond_counter += 1

                    if items_in_line == 0 and bond_counter < bonds_number:
                        self._output_file.write('\n ')
                base += molecule.atoms_number_with_drude

        self._output_file.write('\n')

    def __save_angles_section(self):
        angles = {}
        angles_number = 0

        for (molecule, counter) in self._system.molecules:
            angles[molecule] = molecule.angles
            angles_number += counter * len(angles[molecule])

        self._output_file.write("{:>8} !NTHETA: angles\n".format(angles_number))
        base = 1

        self._output_file.write(" ")
        items_in_line = 0

        angle_counter = 0

        for (molecule, counter) in self._system.molecules:
            angle_list = angles[molecule]

            for i in range(0, counter):
                for (first, middle, last) in angle_list:
                    self._output_file.write("{:7d} {:7d} {:7d} ".format(first + base, middle + base, last + base))
                    items_in_line = (items_in_line + 1) % 3
                    angle_counter += 1

                    if items_in_line == 0 and angle_counter < angles_number:
                        self._output_file.write('\n ')
                base += molecule.atoms_number_with_drude

        self._output_file.write('\n')

    def __save_dihedrals_section(self):
        dihedrals = {}
        dihedrals_number = 0

        for (molecule, counter) in self._system.molecules:
            dihedrals[molecule] = molecule.dihedrals
            dihedrals_number += counter * len(dihedrals[molecule])

        self._output_file.write("{:>8} !NPHI: dihedrals\n".format(dihedrals_number))
        base = 1

        self._output_file.write(" ")
        items_in_line = 0

        dihedral_counter = 0

        for (molecule, counter) in self._system.molecules:
            dihedral_list = dihedrals[molecule]

            for i in range(0, counter):
                for (first, second, third, fourth) in dihedral_list:
                    self._output_file.write(
                        "{:7d} {:7d} {:7d} {:7d} ".format(first + base, second + base, third + base, fourth + base))
                    items_in_line = (items_in_line + 1) % 2
                    dihedral_counter += 1

                    if items_in_line == 0 and dihedral_counter < dihedrals_number:
                        self._output_file.write('\n ')
                base += molecule.atoms_number_with_drude

        self._output_file.write('\n')

    def __save_file_ending(self):
        self._output_file.write("{:>8} !NIMPHI: impropers\n".format(0))
        self._output_file.write("{:>8} !NDON: donors\n".format(0))
        self._output_file.write("{:>8} !NACC: acceptors\n".format(0))
        self._output_file.write("{:>8} !NNB\n".format(0))
        self._output_file.write("{:>8} !NGRP\n".format(0))
        self._output_file.write("{:>8}  {:>8} !NUMLP NUMLPH\n".format(0, 0))
        self._output_file.write("{:>8} !NCRTERM: cross-terms\n".format(0))
        self._output_file.write("{:>8} !NUMANISO\n".format(0))
        self._output_file.write("END\n")


class PDBSaver(FileSaver):
    FIRST_LINE = "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n"
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

        # write header line in PDB
        self._output_file.write(self.FIRST_LINE)

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

                    drude_atom = atom.drude_atom
                    if drude_atom is not None:
                        atom_number += 1
                        self._output_file.write(self.LINE_FORMAT.format(atom_number,
                                                                        drude_atom.namd_symbol,
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
