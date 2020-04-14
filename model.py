from collections import Counter
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
    def __init__(self, symbol, coordinates, charge=0.0, mass=0.0, namd_symbol=None):
        self._symbol = symbol
        self._coordinates = coordinates
        self._charge = charge
        self._mass = mass

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

    @property
    def mass(self):
        return self._mass

    @property
    def namd_symbol(self):
        return self._namd_symbol

    def __eq__(self, other):
        if not isinstance(other, Atom):
            return NotImplemented

        return self._symbol == other._symbol and \
               self._coordinates == other._coordinates and \
               self._charge == other._charge and \
               self._mass == other._mass and \
               self._namd_symbol == other._namd_symbol

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
        self._angles = None

    @property
    def atoms_number(self):
        return len(self._atoms)

    @property
    def atoms(self):
        return self._atoms

    @property
    def residue_name(self):
        return self._residue_name

    @property
    def bonds(self):
        return self._bonds

    @bonds.setter
    def bonds(self, bonds):
        self._bonds = tuple(bonds)

    @property
    def angles(self):
        return self._angles

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
                    neighbour_index = self._atoms.index(neighbour)
                    result.append((atom_index, neighbour_index))
        self._bonds = tuple(result)

    def determine_angles(self):
        if self.bonds is None:
            return

        result = []
        for i in range(0, self.atoms_number):
            filtered_bonds = list(filter(lambda x: i in x, self._bonds))
            left_sides = list(filter(lambda x: x[1] == i, filtered_bonds))
            right_sides = list(filter(lambda x: x[0] == i, filtered_bonds))

            for left in left_sides:
                for right in right_sides:
                    result.append(left + (right[1],))

        self._angles = tuple(result)


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
