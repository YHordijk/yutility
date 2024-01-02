from scm import plams
from TCutility import formula, log, results


class Job:
    '''This is the base Job class used to build more advanced classes such as AMSJob and ORCAJob.
    The base class contains an empty DotDict object that holds the settings. It also provides __enter__ and __exit__ methods to make use of context manager syntax.'''
    def __init__(self):
        self.settings = results.Result()
        self.sbatch = results.Result()
        self._molecule = None
        self.name = None
        self.rundir = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.run()

    def __repr__(self):
        s = '#! /bin/bash\n\n'
        if self.sbatch:
            for key, val in self.sbatch.items():
                key.replace('_', '-')
                if len(key) > 1:
                    s += f'#SBATCH --{key}={val}\n'
                else:
                    s += f'#SBATCH -{key} {val}\n'
            s += '\n'

        return s

    def run(self):
        NotImplemented

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, mol):
        assert isinstance(mol, (str, plams.Molecule)), f'Argument should be a plams.Molecule object or a path'
        
        if isinstance(mol, plams.Molecule):
            self._molecule = mol

        elif isinstance(mol, str):
            self._molecule = plams.Molecule(mol)


class ADFJob(Job):
    def __init__(self):
        super().__init__()
        self._functional = None
        self._basis_set = None

    def __setattr__(self, key, val):
        if key.lower() == 'basis_set':
            self.settings.input.adf.basis.type = val
            self.settings.input.adf.basis.core = 'None'
            return

        if key.lower() == 'task':
            self.settings.input.ams.task = val
            return

        if key.lower() == 'vibrations':
            if val is True:
                self.settings.input.ams.Properties.NormalModes = 'Yes'
                self.settings.input.ams.Properties.PESPointCharacter = 'Yes'
                self.settings.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = tolerance
                self.settings.input.ams.NormalModes.ReScanFreqRange = '-10000000.0 10.0'
            return

        if key.lower() == 'charge':
            self.settings.input.ams.System.Charge = val
            return

        if key.lower() == 'spin_polarization':
            self.settings.input.adf.SpinPolarization = val
            if val != 0:
                self.settings.input.adf.Unrestricted = 'Yes'
            return

        super().__setattr__(key, val)

    @property
    def functional(self):
        return self._functional

    @functional.setter
    def functional(self, functional):
        # before adding the new functional we should clear any previous functional settings
        self.settings.input.adf.pop('XC', None)

        functional = functional.strip()
        self._functional = functional

        # first handle dispersion corrections
        disp_map = {
            '-D4(EEQ)': 'GRIMME4',
            '-D3(BJ)': 'GRIMME3 BJDAMP',
            '-D3BJ': 'GRIMME3 BJDAMP',
            '-D3': 'GRIMME3',
            '-dDsC': 'dDsC',
            '-dUFF': 'UFF',
            '-MBD@rsSC': 'MBD',
            '-D': 'DEFAULT'
        }

        # look through all possible dispersion options
        for name, adf_name in disp_map.items():
            # check if the user requests the dispersion correction
            if functional.endswith(name):
                # correctly set the dispersion correction
                self.settings.input.adf.XC.Dispersion = adf_name
                # remove the correction from the functional name
                # we would use str.removesuffix, but that is only since python 3.9, so we just use slicing instead
                functional = functional[:-len(name)]

        # handle the actual functional part
        functional_map = {
            'LDA': ['VWN', 'PW92'],
            'GGA': ['BLYP', 'BP86', 'GAM', 'HTBS', 'KT1', 'KT2', 'mPW', 'mPBE', 'N12', 'OLYP', 'OPBE', 'PBE', 'PBEsol', 'PW91', 'revPBE', 'RPBE', 'BEE'],
            'Hybrid': ['B3LYP', 'B1LYP', 'B1PW91', 'B3LYP*', 'BHandH', 'BHandHLYP', 'KMLYP', 'MPW1PW', 'MPW1K', 'O3LYP', 'OPBE0', 'PBE0', 'S12h', 'X3LYP', 'HTBS'],
            'MetaGGA': ['M06L', 'MN15-L', 'MVS', 'SCAN', 'revTPSS', 'SSB', 'TASKxc', 'TPSS', 'r2SCAN-3c'],
            'LibXC': ['rSCAN', 'revSCAN', 'r2SCAN'] + ['LCY-BLYP', 'LCY-BP86', 'LCY-PBE', 'CAM-B3LYP', 'CAMY-B3LYP', 'HSE03', 'HSE06', 'M11', 'MN12-SX', 'N12-SX', 'WB97', 'WB97X'] + ['revSCAN0'],
            'DoubleHybrid': ['rev-DOD-PBEP86-D4', 'rev-DOD-BLYP-D4', 'rev-DOD-PBE-D4', 'B2PLYP', 'B2GPPLYP'],
            'MetaHybrid': ['MN15', 'M06', 'M06-2X', 'M06-HF', 'TPSSH'],
            'model': ['SAOP'],
        }

        for category, functionals in functional_map.items():
            if functional in functionals:
                self.settings.input.adf.XC[category] = functional
                return

        # define some preset functionals
        if functional == 'BMK':
            self.settings.input.adf.XC.LibXC = 'HYB_MGGA_X_BMK GGA_C_BMK'
            return

        # LDA is the standard functional
        if functional == 'LDA':
            return

        log.warn(f'XC-functional {functional} not defined. Defaulting to using LibXC.')
        self.settings.input.adf.XC.LibXC = functional


if __name__ == '__main__':
    from pprint import pprint

    with ADFJob() as job:
        job.molecule = '/Users/yumanhordijk/PhD/ychem/calculations/00f8a789a08ae1a976ebd9eb09652290720e60a029bc13b84dc573b2e19a8a94/input_mol.xyz'
        job.sbatch.p = 'tc'
        job.sbatch.ntasks_per_node = 15

        job.functional = 'BMK'
        job.functional = 'PW92'

        job.charge = 10
        job.spin_polarization = 1

    pprint(job)
