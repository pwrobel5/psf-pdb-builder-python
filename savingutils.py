import abc


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
        self._output_file.close()

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
