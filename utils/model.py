from collections import Counter
from functools import reduce
from math import sqrt


class Coordinates:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self._x = x
        self._y = y
        self._z = z

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    def norm(self):
        return sqrt(self._x ** 2 + self._y ** 2 + self._z ** 2)

    def calculate_distance(self, other):
        distance = Coordinates(self.x - other.x,
                               self.y - other.y,
                               self.z - other.z)
        return distance.norm()

    def __eq__(self, other):
        if not isinstance(other, Coordinates):
            return NotImplemented

        return self._x == other._x and self._y == other._y and self._z == other._z

    def __hash__(self):
        return hash((self._x, self._y, self._z))

    def __str__(self):
        return str(self._x) + ' ' + str(self._y) + ' ' + str(self._z)


class Atom:
    DRUDE_CHARGE_FACTOR = -0.0548768646057431

    def __init__(self, symbol, coordinates, charge=0.0, mass=0.0, namd_symbol=None):
        self._symbol = symbol
        self._coordinates = coordinates
        self._charge = charge
        self._mass = mass
        self._drude_atom = None

        if namd_symbol is None:
            self._namd_symbol = symbol
        else:
            self._namd_symbol = namd_symbol

    @property
    def symbol(self):
        return self._symbol

    @property
    def coordinates(self):
        return self._coordinates

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, charge):
        self._charge = charge

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, mass):
        self._mass = mass

    @property
    def namd_symbol(self):
        return self._namd_symbol

    @namd_symbol.setter
    def namd_symbol(self, namd_symbol):
        self._namd_symbol = namd_symbol

    @property
    def drude_atom(self):
        return self._drude_atom

    @drude_atom.setter
    def drude_atom(self, drude_atom):
        self._drude_atom = drude_atom

    def __eq__(self, other):
        if not isinstance(other, Atom):
            return NotImplemented

        return self._symbol == other._symbol and \
               self._coordinates == other._coordinates and \
               self._charge == other._charge and \
               self._mass == other._mass and \
               self._namd_symbol == other._namd_symbol

    def create_drude_atom(self, polarizability, bond_const=1000, mass=0.4):
        drude_namd_symbol = "D" + self._namd_symbol[:2]
        drude_charge = self.DRUDE_CHARGE_FACTOR * sqrt(polarizability * bond_const)
        drude_atom = Atom(self._symbol, self._coordinates, drude_charge, mass, drude_namd_symbol)

        self._mass -= mass
        self._charge -= drude_charge
        self._drude_atom = drude_atom

    def __hash__(self):
        return hash((self._symbol, self._coordinates, self._charge, self._mass, self._namd_symbol))

    def __str__(self):
        return self._namd_symbol + '(' + self._symbol + ') ' + str(self._coordinates) + ' charge: ' + str(
            self._charge) + ', mass: ' + str(self._mass)


class Molecule:
    def __init__(self, atoms, residue_name="MOL"):
        self._atoms = atoms
        self._residue_name = residue_name
        self._bonds = None
        self._drude_bonds = None
        self._angles = None
        self._dihedrals = None
        self._shifts = {i: 0 for i in range(0, self.atoms_number)}

    @property
    def atoms_number(self):
        return len(self._atoms)

    @property
    def atoms_number_with_drude(self):
        return reduce(lambda a, b: a + (1 if b.drude_atom is None else 2), self._atoms, 0)

    @property
    def atoms(self):
        return self._atoms

    @property
    def residue_name(self):
        return self._residue_name

    @residue_name.setter
    def residue_name(self, residue_name):
        self._residue_name = residue_name

    @property
    def bonds(self):
        return self._bonds

    @bonds.setter
    def bonds(self, bonds):
        self._bonds = tuple(bonds)

    @property
    def drude_bonds(self):
        return self._drude_bonds

    @drude_bonds.setter
    def drude_bonds(self, drude_bonds):
        self._drude_bonds = tuple(drude_bonds)

    @property
    def angles(self):
        return self._angles

    @angles.setter
    def angles(self, angles):
        self._angles = tuple(angles)

    @property
    def dihedrals(self):
        return self._dihedrals

    @dihedrals.setter
    def dihedrals(self, dihedrals):
        self._dihedrals = tuple(dihedrals)

    @property
    def shifts(self):
        return self._shifts

    def add_drude_atom(self, atom, polarizability, drude_bonds_list):
        atom.create_drude_atom(polarizability)
        index = self._atoms.index(atom)
        shifted_index = index + self.shifts[index]
        drude_bonds_list.append((shifted_index, shifted_index + 1))

        for i in range(index + 1, self.atoms_number):
            self._shifts[i] += 1

    def __eq__(self, other):
        if not isinstance(other, Molecule):
            return NotImplemented

        return Counter(self._atoms) == Counter(other._atoms) and self._residue_name == other._residue_name

    def __hash__(self):
        return hash((tuple(self._atoms), self._residue_name))

    def __str__(self):
        result = 'Molecule: ' + self._residue_name + '\n'

        for atom in self._atoms:
            result += str(atom) + '\n'

        return result

    def determine_bonds(self, threshold=1.70):
        result = []

        for atom in self._atoms[:-1]:
            atom_index = self._atoms.index(atom)
            for neighbour in self._atoms[atom_index + 1:]:
                if atom.coordinates.calculate_distance(neighbour.coordinates) <= threshold:
                    atom_index += self._shifts[atom_index]
                    neighbour_index = self._atoms.index(neighbour)
                    neighbour_index += self._shifts[neighbour_index]
                    result.append((atom_index, neighbour_index))
        self._bonds = tuple(result)

    def determine_angles(self):
        if self.bonds is None:
            return

        result = []
        for i in range(0, self.atoms_number):
            filtered_bonds = filter(lambda x: i in x, self._bonds)
            left_sides = filter(lambda x: x[1] == i, filtered_bonds)
            right_sides = filter(lambda x: x[0] == i, filtered_bonds)

            for left in left_sides:
                for right in right_sides:
                    result.append(left + (right[1],))

        self._angles = tuple(result)

    def determine_dihedrals(self):
        if self.angles is None:
            return

        result = []
        for i in range(0, self.atoms_number):
            left_sides = list(filter(lambda x: x[2] == i, self._angles))
            right_sides = list(filter(lambda x: x[0] == i, self._bonds))

            for left in left_sides:
                for right in right_sides:
                    result.append(left + (right[1],))

        self._dihedrals = tuple(result)


class System:
    def __init__(self, molecules, xyz_file_name, segment_id="IL"):
        self._molecules = molecules
        self._xyz_file_name = xyz_file_name
        self._segment_id = segment_id

    @property
    def molecules(self):
        return self._molecules

    @property
    def xyz_file_name(self):
        return self._xyz_file_name

    @property
    def segment_id(self):
        return self._segment_id

    @property
    def atoms_number(self):
        result = 0
        for (molecule, molecule_count) in self._molecules:
            result += molecule.atoms_number * molecule_count

        return result

    @property
    def atoms_number_with_drude(self):
        result = 0
        for (molecule, molecule_count) in self._molecules:
            result += molecule.atoms_number_with_drude * molecule_count

        return result
