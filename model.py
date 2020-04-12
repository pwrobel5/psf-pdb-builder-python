from collections import Counter


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

    def __eq__(self, other):
        if not isinstance(other, Coordinates):
            return NotImplemented

        return self._x == other._x and self._y == other._y and self._z == other._z

    def __hash__(self):
        return hash((self._x, self._y, self._z))

    def __str__(self):
        return str(self._x) + ' ' + str(self._y) + ' ' + str(self._z)


class Atom:
    def __init__(self, symbol, coordinates, charge=0.0, mass=0.0, namd_symbol=""):
        self._symbol = symbol
        self._coordinates = coordinates
        self._charge = charge
        self._mass = mass
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

    def __eq__(self, other):
        if not isinstance(other, Molecule):
            return NotImplemented

        return Counter(self._atoms) == Counter(other._atoms) and self._residue_name == other._residue_name

    def __hash__(self):
        return hash((self._atoms, self._residue_name))

    def __str__(self):
        result = 'Molecule: ' + self._residue_name + '\n'

        for atom in self._atoms:
            result += str(atom) + '\n'

        return result
