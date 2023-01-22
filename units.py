import numpy as np
import scm.plams as plams
from yutility import dictfunc

si_prefixes = {-24: 'y',
               -21: 'z',
               -18: 'a',
               -15: 'f',
               -12: 'p',
               -9: 'n',
               -6: 'u',
               -3: 'm',
               # -2:  'c',
               # -1:  'd',
               0: '',
               # 1:   'da',
               # 2:   'h',
               3: 'k',
               6: 'M',
               9: 'G',
               12: 'T',
               15: 'P',
               18: 'E',
               21: 'Z',
               24: 'Y'}

si_prefixes_inv = dictfunc.invert(si_prefixes)


def get_si_prefix(x, only_major=True):
    # order = np.floor(np.log10(x))
    # splits = str(x).split('.')
    science = f'{x:e}'
    order = int(science.split('e')[1])
    if order not in si_prefixes:
        order -= order % 3
    if only_major:
        order -= order % 3
    return si_prefixes.get(order, ''), 10**-order


def get_scientific(x):
    order = np.floor(np.log10(x))
    return f'·10^{int(order)}', 10**-order


def convert_time(t):
    # converts time in seconds to a better unit
    if t / 60 < 1:
        return 's', 1
    if t / 60 / 60 < 1:
        return 'm', 1 / 60
    return 'h', 1 / 60 / 60


######################
def hartree2kcalmol(v):
    return plams.Units.convert(v, 'Hartree', 'kcal/mol')


h2k = hartree2kcalmol


def hartree2eV(v):
    return plams.Units.convert(v, 'Hartree', 'eV')


h2e = hartree2eV


#####################
def bohr2angstrom(v):
    return plams.Units.convert(v, 'bohr', 'angstrom')


b2a = bohr2angstrom


class Unit:
    def __init__(self, unit='', inverse=False):
        self.unit = unit
        self.best_unit = unit
        self.inverse = inverse
        self.conversion_factors = {}

    def string(self, val, rounding=2, use_si=True):
        val, un = self.convert(val, use_si=use_si)
        return f'{val:.{rounding}f} {un}'

    def __repr__(self):
        return f'{type(self).__name__}:{self.unit}'

    def convert_si(self, val):
        # converts a value to readable value
        # ex. convert(100, kind='SI') == 100 u
        #     convert(1000, kind='SI') == 1 ku
        return get_si_prefix(val)

    def convert_scientific(self, val):
        return get_scientific(val)

    def convert_unit(self, val, unit):
        if unit not in self.conversion_factors:
            raise Exception(
                f'Unit {unit} not in conversion_factors. Please add it in __init__.')
        if self.unit not in self.conversion_factors:
            raise Exception(
                f'Unit {self.unit} not in conversion_factors. Please add it in __init__.')
        if self.inverse:
            return val * \
                self.conversion_factors[self.unit] / \
                self.conversion_factors[unit]
        else:
            return val / \
                self.conversion_factors[self.unit] * \
                self.conversion_factors[unit]

    def convert(self, val, unit=None, use_si=False, si_only_major=True):
        if unit is None:
            unit = self.unit
        if val is None:
            return
        sunit_si, self.unit = self._parse_unit(self.unit)
        tunit_si, tunit = self._parse_unit(unit)
        sfac, tfac = 10**si_prefixes_inv[sunit_si], 10**si_prefixes_inv[tunit_si]
        val = sfac * val
        tval = self.convert_unit(val, tunit)
        if use_si:
            pre, fac = self.convert_si(tval)
        else:
            pre, fac = tunit_si, 1 / tfac
        # print(fac)
        if self.inverse:
            return tval * fac, f'({pre}{tunit})^-1'
        return tval * fac, f'{pre}{tunit}'

    def _parse_unit(self, unit):
        units_ = list(self.conversion_factors.keys())
        match = [unit_ for unit_ in units_ if unit.endswith(unit_)][0]
        if len(match) == len(unit):
            si = ''
        else:
            si = unit[:-len(match)]
        return si, match


class Energy(Unit):
    def __init__(self, unit='Ha', best_unit='kcal/mol', inverse=False):
        self.unit = unit
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            'Ha': 1,
            'ha': 1,
            'hartree': 1,
            'Hartree': 1,
            'cal/mol': 627509.6080305927,
            'J/mol': 2625500.2,
            'eV': 27.211399,
            'cm^-1': 219474.63
        }


class Length(Unit):
    def __init__(self, unit='Å', best_unit='Å', inverse=False):
        self.unit = unit
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            'm': 1e-10,
            'Å': 1
        }


class Frequency(Unit):
    def __init__(self, unit='cm^-1', best_unit='cm^-1', inverse=False):
        self.unit = unit
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            'cm^-1': 1,
            'Hz': 29979245800,
        }


class VibIntensity(Unit):
    def __init__(self, unit='km/mol', best_unit='km/mol', inverse=False):
        self.unit = unit
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            'km/mol': 1,
        }


class Time(Unit):
    def __init__(self, unit='s', best_unit='s', inverse=False):
        self.unit = unit
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            's': 1,
            'm': 1 / 60,
            'h': 1 / 60 / 60,
            'd': 1 / 60 / 60 / 24
        }

    def convert(self, val, unit=None, use_si=False, si_only_major=True):
        if unit is None:
            val = self.convert_unit(val, 's')
            self.unit = 's'
            if val < 60:
                tunit = 's'
            elif val < 60 * 60:
                tunit = 'm'
            elif val < 60 * 60 * 24:
                tunit = 'h'
            else:
                tunit = 'd'
            val = self.convert_unit(val, tunit)
            self.unit = tunit
        return super().convert(val, unit=unit, use_si=use_si, si_only_major=si_only_major)


class UnitLess(Unit):
    def __init__(self, unit='', best_unit='', inverse=False):
        self.unit = ''
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            '': 1,
            'mol': 1 / 6.0221415e23}


class Binary(Unit):
    def __init__(self, unit='B', best_unit='B', inverse=False):
        self.unit = unit
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            'B': 1,
            'b': 8,
        }


class Electrons(Unit):
    def __init__(self, unit='e^-', best_unit='e^-', inverse=False):
        self.unit = unit
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            'e^-': 1,
        }


class Electrons2(Unit):
    def __init__(self, unit='(e^-)^2', best_unit='(e^-)^2', inverse=False):
        self.unit = unit
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            '(e^-)^2': 1,
        }


class Electrons2Energy(Unit):
    def __init__(self, unit='(e^-)^2Ha', best_unit='(e^-)^2eV', inverse=False):
        self.unit = unit
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            '(e^-)^2Ha': 1,
            '(e^-)^2eV': 27.211399,
        }


class Electrons2EnergyInv(Unit):
    def __init__(self, unit='(e^-)^2/Ha',
                 best_unit='(e^-)^2/eV', inverse=False):
        self.unit = unit
        self.best_unit = best_unit
        self.inverse = inverse
        self.conversion_factors = {
            '(e^-)^2/Ha': 1,
            '(e^-)^2/eV': 1 / 27.211399,
        }

def from_string(string):
    subclassname, unitname = string.split(':')
    for subclass in Unit.__subclasses__():
        if subclass.__name__ == subclassname:
            return subclass(unitname)


if __name__ == '__main__':
    # print(get_si_prefix(1))
    # print(get_scientific(1e2))
    # print(get_si_prefix(1e6+1))
    # print(get_si_prefix(1e9))

    # print(Time('h').convert_unit(90, 'm'))
    # print(Length('m').convert_unit(1, 'm'))
    # print(Binary('Mb').convert(8, 'MB', use_si=False))

    print(Energy('Ha', inverse=True).convert(1, 'eV'))
    unit = Energy('Ha', inverse=True)
    print(unit)

    unit = from_string('Energy:Ha')
    print(type(unit).__name__, unit.unit)