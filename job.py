from scm import plams
from TCutility import log, results, formula
import os


gr = plams.GridRunner(parallel=True, maxjobs=0, grid='auto')


class Job:
    '''This is the base Job class used to build more advanced classes such as AMSJob and ORCAJob.
    The base class contains an empty DotDict object that holds the settings. It also provides __enter__ and __exit__ methods to make use of context manager syntax.'''
    def __init__(self):
        self.settings = results.Result()
        self._sbatch = None
        self._molecule = None
        self.name = 'calc'
        self.rundir = 'tmp'

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.run()

    def __repr__(self):
        s = '#!/bin/bash\n'
        self._sbatch.prune()
        if self._sbatch:
            for key, val in self._sbatch.items():
                key = key.replace('_', '-')
                if len(key) > 1:
                    s += f'#SBATCH --{key}={val}\n'
                else:
                    s += f'#SBATCH -{key} {val}\n'
            s += '\n'
        return s

    def sbatch(self, **kwargs):
        self._sbatch = results.Result(kwargs)

    def run(self):
        NotImplemented

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, mol):
        assert isinstance(mol, (str, plams.Molecule)), f'Argument should be a plams.Molecule object or a path, not {type(mol)}'
        
        if isinstance(mol, plams.Molecule):
            self._molecule = mol

        elif isinstance(mol, str):
            self._molecule = plams.Molecule(mol)
            log.info(f'Succesfully loaded molecule {formula.molecule(self._molecule)} from path.')


class ADFJob(Job):
    def __init__(self):
        super().__init__()
        self.functional('LDA')
        self.basis_set('DZ')
        self.single_point()

    def __str__(self):
        return f'{self._task}({self._functional}/{self._basis_set}), running in {os.path.join(os.path.abspath(self.rundir), self.name)}'

    def basis_set(self, typ, core='None'):
        self._basis_set = typ
        self.settings.input.adf.basis.type = typ
        self.settings.input.adf.basis.core = core

    def single_point(self):
        self._task = 'SP'
        self.settings.input.ams.pop('TransitionStateSearch', None)
        self.settings.input.ams.task = 'SinglePoint'

    def transition_state(self, distances=None, angles=None, coordinates=None, dihedrals=None, ModeToFollow=1):
        self._task = 'TS'
        self.settings.input.ams.task = 'TransitionStateSearch'
        self.settings.input.ams.TransitionStateSearch.ReactionCoordinate.ModeToFollow = ModeToFollow

        if distances is not None:
            self.settings.input.ams.TransitionStateSearch.ReactionCoordinate.Distance = [" ".join([str(x) for x in dist]) for dist in distances]
        if angles is not None:
            self.settings.input.ams.TransitionStateSearch.ReactionCoordinate.Angle = [" ".join([str(x) for x in ang]) for ang in angles]
        if coordinates is not None:
            self.settings.input.ams.TransitionStateSearch.ReactionCoordinate.Coordinate = [" ".join([str(x) for x in coord]) for coord in coordinates]
        if dihedrals is not None:
            self.settings.input.ams.TransitionStateSearch.ReactionCoordinate.Dihedral = [" ".join([str(x) for x in dihedral]) for dihedral in dihedrals]

    def optimization(self):
        self._task = 'GO'
        self.settings.input.ams.pop('TransitionStateSearch', None)
        self.settings.input.ams.task = 'GeometryOptimization'

    def vibrations(self, enable=True, PESPointCharacter=True, NegativeFrequenciesTolerance=-5, ReScanFreqRange='-10000000.0 10.0'):
        self.settings.input.ams.Properties.NormalModes = 'Yes' if enable else 'No'
        self.settings.input.ams.Properties.PESPointCharacter = 'Yes' if enable else 'No'
        self.settings.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = NegativeFrequenciesTolerance
        self.settings.input.ams.NormalModes.ReScanFreqRange = ReScanFreqRange

    def charge(self, val):
        self.settings.input.ams.System.Charge = val

    def spin_polarization(self, val):
        self.settings.input.adf.SpinPolarization = val
        if val != 0:
            self.settings.input.adf.Unrestricted = 'Yes'

    def multiplicity(self, val):
        self.settings.input.adf.SpinPolarization = (val - 1)//2
        if val != 0:
            self.settings.input.adf.Unrestricted = 'Yes'

    def unrestricted(self, val):
        self.settings.input.adf.Unrestricted = 'Yes' if val else 'No'

    def quality(self, val='Normal'):
        self.settings.input.adf.NumericalQuality = val

    def functional(self, val):
        # before adding the new functional we should clear any previous functional settings
        self.settings.input.adf.pop('XC', None)

        functional = val.strip()
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

    def solvent(self, name=None, eps=None, rad=None, use_klamt=False):
        self.settings.input.adf.Solvation.Surf = 'Delley'
        solv_string = ''
        if name:
            solv_string += f'name={name} '
        else:
            self.settings.input.adf.Solvation.Solv = f'eps={eps} rad={rad} '
        if use_klamt:
            solv_string += 'cav0=0.0 cav1=0.0'
        else:
            solv_string += 'cav0=0.0 cav1=0.0067639'
        self.settings.input.adf.Solvation.Solv = solv_string

        self.settings.input.adf.Solvation.Charged = 'method=CONJ corr'
        self.settings.input.adf.Solvation['C-Mat'] = 'POT'
        self.settings.input.adf.Solvation.SCF = 'VAR ALL'
        self.settings.input.adf.Solvation.CSMRSP = None

        if use_klamt:
            radii = {
                'H': 1.30,
                'C': 2.00,
                'N': 1.83,
                'O': 1.72,
                'F': 1.72,
                'Si': 2.48,
                'P': 2.13,
                'S': 2.16,
                'Cl': 2.05,
                'Br': 2.16,
                'I': 2.32
            }
            self.settings.input.adf.solvation.radii = radii

    def run(self):
        plams.init()
        sett = self.settings.as_plams_settings()
        sett.keep = ['-', 't21.*', 'CreateAtoms.out', '*.dill']
        job = plams.AMSJob(name=self.name, molecule=self.molecule, settings=sett)
        job.run(jobrunner=gr, queue='tc', n=32, J=self.name)
        plams.finish()


class OrcaJob(Job):
    def __init__(self):
        super().__init__()
        self.settings.main = {'LARGEPRINT'}
        self._charge = 0
        self._multiplicity = 1
        self.memory = None
        self.processes = None

    def __casefold_main(self):
        self.settings.main = {key.casefold() for key in self.settings.main}

    def __remove_task(self):
        self.__casefold_main()
        [self.settings.main.discard(task) for task in ['sp', 'opt', 'tsopt', 'neb-ts']]

    def single_point(self):
        self.__remove_task()
        self.settings.main.add('sp')

    def transition_state(self):
        self.__remove_task()
        self.settings.main.add('optts')

    def optimization(self):
        self.__remove_task()
        self.settings.main.add('opt')

    def vibrations(self, enable=True, numerical=False):
        self.__casefold_main()
        self.settings.main.discard('numfreq')
        self.settings.main.discard('freq')
        if enable:
            if numerical:
                self.settings.main.append('numfreq')
            else:
                self.settings.main.append('freq')

    def charge(self, val):
        self._charge = val

    def spin_polarization(self, val):
        self._multiplicity = 2 * val + 1

    def multiplicity(self, val):
        self._multiplicity = val

    # def solvent(self, name=None, eps=None, rad=None, model='COSMO'):
    #     assert model in ['CPCM', 'SMD', 'COSMO'], f'Solvation model must be one of [CPCM, SMD, COSMO], not {model}.'

    #     self.__casefold_main()
    #     self.settings.main.discard('cpcm')
    #     self.settings.main.discard('smd')

    def get_memory_usage(self):
        mem = self.memory or self._sbatch.mem or None

        ntasks = self.processes
        if ntasks is None:
            if self._sbatch.n:
                ntasks = self._sbatch.n
            if self._sbatch.ntasks:
                ntasks = self._sbatch.ntasks
            if self._sbatch.ntasks_per_node:
                ntasks = self._sbatch.ntasks_per_node * self._sbatch.get('N', 1) * self._sbatch.get('nodes', 1)

        return mem, ntasks

    @property
    def input(self):
        # set the correct memory usage and processes
        natoms = len(self.molecule)
        mem, ntasks = self.get_memory_usage()
        if ntasks and mem:
            ntasks = min(ntasks, (natoms - 1) * 3)
            self.settings.PAL.nprocs = ntasks
            self.settings.maxcore = int(mem / ntasks * 0.75)
        else:
            log.warn('MaxCore and nprocs not specified. Please use SBATCH settings or set job.processes and job.memory.')

        ret = ''
        for key in self.settings.main:
            ret += f'!{key}\n'
        ret += '\n'

        for option, block in self.settings.items():
            if option == 'main':
                continue

            if isinstance(block, results.Result):
                ret += f'%{option}\n'

                for key, val in block.items():
                    ret += f'    {key} {val}\n'

                ret += 'END\n\n'
            else:
                ret += f'%{option} {block}\n'

        ret += '\n'

        ret += f'* xyz {self._charge} {self._multiplicity}\n'
        for atom in self.molecule:
            ret += f'    {atom.symbol:2} {atom.x: >13f} {atom.y: >13f} {atom.z: >13f}\n'
        ret += '*\n'

        return ret

    def write(self):
        os.makedirs(os.path.join(self.rundir, self.name), exist_ok=True)
        with open(os.path.join(self.rundir, self.name, 'ORCA.inp'), 'w+') as inp:
            inp.write(self.input)

        with open(os.path.join(self.rundir, self.name, f'{self.name}.sh'), 'w+') as run:
            run.write(repr(self))


if __name__ == '__main__':
    with ADFJob() as job:
        job.molecule = r"D:\Users\Yuman\Desktop\PhD\TCutility\test\fixtures\chloromethane_sn2_ts\ts sn2.results\output.xyz"
        job.sbatch(p='tc', ntasks_per_node=15)

        job.functional('BM12K')
        job.charge(0)
        job.spin_polarization(1)
        job.transition_state()
        job.optimization()
        job.solvent('Ethanol')
