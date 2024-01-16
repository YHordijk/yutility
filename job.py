from scm import plams
from tcutility import log, results, formula, slurm
import subprocess as sp
import os
from typing import Union

j = os.path.join


class Job:
    '''This is the base Job class used to build more advanced classes such as AMSJob and ORCAJob.
    The base class contains an empty DotDict object that holds the settings. It also provides __enter__ and __exit__ methods to make use of context manager syntax.'''
    def __init__(self, test_mode=False):
        self.settings = results.Result()
        self._sbatch = results.Result()
        self._molecule = None
        self.name = 'calc'
        self.rundir = 'tmp'
        self.test_mode = test_mode
        self._preambles = []
        self._postambles = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self.can_skip():
            log.info(f'Skipping calculation {j(self.rundir, self.name)}, it is already finished or currently pending or running.')
            return
        self.run()

    def can_skip(self):
        res = results.read(j(self.rundir, self.name))
        return not res.status.fatal

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
        for key, value in kwargs.items():
            if key == 'dependency' and 'dependency' in self._sbatch:
                value = self._sbatch['dependency'] + ',' + value
            self._sbatch[key] = value

    def get_sbatch_command(self):
        self._sbatch.prune()
        c = 'sbatch '
        for key, val in self._sbatch.items():
            key = key.replace('_', '-')
            if len(key) > 1:
                c += f'--{key}={val} '
            else:
                c += f'-{key} {val} '
        return c

    def run(self):
        NotImplemented

    def add_preamble(self, line):
        self._preambles.append(line)

    def add_postamble(self, line):
        self._postambles.append(line)

    def dependency(self, otherjob: 'Job'):
        if hasattr(otherjob, 'slurm_job_id'):
            self.sbatch(dependency=f'afterok:{otherjob.slurm_job_id}')

    @property
    def workdir(self):
        return j(os.path.abspath(self.rundir), self.name)


def get_available_functionals():
    '''
    Function that returns a dictionary of all available XC-functionals.
    '''
    def set_dispersion(func):
        disp_map = {
            '-D4': 'GRIMME4',
            '-D3(BJ)': 'GRIMME3 BJDAMP',
            '-D3BJ': 'GRIMME3 BJDAMP',
            '-D3': 'GRIMME3',
            '-dDsC': 'dDsC',
            '-dUFF': 'UFF',
            '-MBD': 'MBD',
            '-MBD@rsSC': 'MBD',
            '-D': 'DEFAULT'
        }

        # set some default values for useful parameters
        func.name_no_disp = func.name
        func.dispersion = None
        func.dispersion_name = None

        # go through every dispersion correction and check if we are using it
        for disp_suffix, disp_name in disp_map.items():
            if func.name.endswith(disp_suffix):
                # dispersion will be the suffix without the -
                func.dispersion = disp_suffix[1:]
                func.dispersion_name = disp_name

                # check if the functional already includes the dispersion correction
                if func.includes_disp:
                    break

                # else we set the name of the functional without dispersion for later
                func.name_no_disp = func.name[:-len(disp_suffix)]
                # get the dispersion settings for ADF. Try to get custom values if they were provided.
                func.adf_settings.XC.Dispersion = func.disp_params or disp_name
                break

    def set_functional(func):
        # set the functional settings for ADF
        # first go through some special functionals that require special settings
        if func.name_no_disp == 'BMK':
            func.adf_settings.XC.LibXC = 'HYB_MGGA_X_BMK GGA_C_BMK'
            return

        if func.name in ['LCY-BLYP', 'LCY-BP86', 'LCY-PBE']:
            func.adf_settings.XC.GGA = func.name.split('-')[1]
            func.adf_settings.XC.RANGESEP = ''
            func.adf_settings.XC.xcfun = ''
            return

        if func.name in ['CAMY-B3LYP']:
            func.adf_settings.XC.Hybrid = 'CAMY-B3LYP'
            func.adf_settings.XC.RANGESEP = ''
            func.adf_settings.XC.xcfun = ''
            return

        if func.name == 'GGA:SSB-D':
            func.adf_settings.XC.GGA = 'SSB-D'
            return

        if func.name == 'MetaGGA:SSB-D':
            func.adf_settings.XC.MetaGGA = 'SSB-D'
            return

        if func.name == 'MP2':
            func.adf_settings.XC.MP2 = ''
            return

        if func.name in ['SOS-MP2', 'SCS-MP2']:
            func.adf_settings.XC.MP2 = ''
            func.adf_settings.XC.EmpiricalScaling = func.name[:-4]
            return

        # the normal functionals are defined based on their category, or selected from libxc
        if func.use_libxc:
            func.adf_settings.XC.LibXC = func.name_no_disp
        else:
            func.adf_settings.XC[func.category] = func.name_no_disp

    # gather all data about available functionals
    functionals = results.Result()  # store all info in this dict

    with open(j(os.path.split(__file__)[0], 'available_functionals.txt')) as file:
        lines = file.readlines()

    for line in lines:
        if not line.strip():
            continue

        # functional names are given starting with -
        # category names without -
        if not line.startswith('- '):
            curr_category = line.strip()
            continue

        # store data about the func in a dict
        func = results.Result()
        func.category = curr_category

        # separate the functional name from the line
        functional_name = line[2:].split('!')[0].split(',')[0].strip()
        func.name = functional_name

        # check if custom params were given for dispersion
        if 'GRIMME' in line:
            func.disp_params = line.split('!')[0].split(',')[1].strip().strip("'")

        func.use_libxc = '!libxc' in line
        func.includes_disp = '!includesdisp' in line

        set_dispersion(func)
        set_functional(func)

        functionals[functional_name] = func

    return functionals

functionals = get_available_functionals()

# from pprint import pprint
# print('BMK-D3(BJ)' in functionals)
# pprint(functionals['BMK-D3(BJ)'])
# pprint(get_available_functionals()['HartreeFock-D4'])
# pprint(get_available_functionals()['MetaGGA:SSB-D'])


class ADFJob(Job):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._functional = None
        self.basis_set('TZ2P')
        self.quality('Good')
        self.SCF_converge(1e-8)
        self.single_point()
        self.solvent('vacuum')

    def __str__(self):
        return f'{self._task}({self._functional}/{self._basis_set}), running in {os.path.join(os.path.abspath(self.rundir), self.name)}'

    def basis_set(self, typ: str = 'TZ2P', core: str = 'None'):
        '''
        Set the basis-set type and frozen core approximation for this calculation.
        Note: If the selected functional is the r2SCAN-3c functional, then the basis-set will be set to mTZ2P.

        Args:
            typ: the type of basis-set to use. Default is TZ2P.
            core: the size of the frozen core approximation. Default is None.
        '''
        if self._functional == 'r2SCAN-3c' and typ != 'mTZ2P':
            log.warn(f'Basis set {typ} is not allowed with r2SCAN-3c, switching to mTZ2P.')
            typ = 'mTZ2P'
        self._basis_set = typ
        self.settings.input.adf.basis.type = typ
        self.settings.input.adf.basis.core = core

    def single_point(self):
        '''
        Set the task of the job to single point.
        '''
        self._task = 'SP'
        self.settings.input.ams.task = 'SinglePoint'

    def transition_state(self, distances: list = None, angles: list = None, dihedrals: list = None, ModeToFollow: int = 1):
        '''
        Set the task of the job to transition state search. Optionally you can give some TS coordinates to accelerate convergence.
        By default also calculates the normal modes after convergence.

        Args:
            distances: sequence of tuples or lists containing [atom_index1, atom_index2, factor]. Atom indices start at 1. 
            angles: sequence of tuples or lists containing [atom_index1, atom_index2, atom_index3, factor]. Atom indices start at 1.
            dihedrals: sequence of tuples or lists containing [atom_index1, atom_index2, atom_index3, atom_index4, factor]. Atom indices start at 1.
            ModeToFollow: the vibrational mode to follow during optimization.
        '''
        self._task = 'TS'
        self.settings.input.ams.task = 'TransitionStateSearch'

        self.settings.input.ams.TransitionStateSearch.ReactionCoordinate.ModeToFollow = ModeToFollow

        if distances is not None:
            self.settings.input.ams.TransitionStateSearch.ReactionCoordinate.Distance = [" ".join([str(x) for x in dist]) for dist in distances]
        if angles is not None:
            self.settings.input.ams.TransitionStateSearch.ReactionCoordinate.Angle = [" ".join([str(x) for x in ang]) for ang in angles]
        if dihedrals is not None:
            self.settings.input.ams.TransitionStateSearch.ReactionCoordinate.Dihedral = [" ".join([str(x) for x in dihedral]) for dihedral in dihedrals]

        # for TS searches we quickly calculate the hessian with DFTB
        self.settings.input.ams.GeometryOptimization.InitialHessian.Type = 'CalculateWithFasterEngine'
        self.vibrations(True)  # also calculate vibrations by default

    def optimization(self):
        '''
        Set the task of the job to transition state search. By default also calculates the normal modes after convergence.
        '''
        self._task = 'GO'
        self.settings.input.ams.task = 'GeometryOptimization'
        self.vibrations(True)

    def vibrations(self, enable: bool = True, PESPointCharacter: bool = True, NegativeFrequenciesTolerance: float = -5, ReScanFreqRange: tuple[float ,float] = [-10000000.0, 10.0]):
        '''
        Set the calculation of vibrational modes. 

        Args:
            enable: whether to calculate the vibrational modes.
            PESPointCharacter: whether to report the PES character in the output.
            NegativeFrequenciesTolerance: the tolerance for negative modes. Modes above this value will not be counted as imaginary. Use this option when you experience a lot of numerical noise.
            ReScanFreqRange: the rescan range. Any mode that has a frequency in this range will be refined.
        '''
        self.settings.input.ams.Properties.NormalModes = 'Yes' if enable else 'No'
        self.settings.input.ams.Properties.PESPointCharacter = 'Yes' if enable else 'No'
        self.settings.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = NegativeFrequenciesTolerance
        self.settings.input.ams.NormalModes.ReScanFreqRange = ' '.join([str(x) for x in ReScanFreqRange])

    def charge(self, val: int):
        '''
        Set the charge of the system.
        '''
        self.settings.input.ams.System.Charge = val

    def spin_polarization(self, val: int):
        '''
        Set the spin-polarization of the system. If the value is not zero the calculation will also be unrestricted.
        '''
        self.settings.input.adf.SpinPolarization = val
        if val != 0:
            self.settings.input.adf.Unrestricted = 'Yes'

    def multiplicity(self, val: int):
        '''
        Set the multiplicity of the system. If the value is not one the calculation will also be unrestricted.
        We use the following values:
            1: singlet
            2: doublet
            3: triplet
            ...
        The multiplicity is equal to 2*S+1 for spin-polarization of S.
        '''
        self.settings.input.adf.SpinPolarization = (val - 1)//2
        if val != 1:
            self.settings.input.adf.Unrestricted = 'Yes'

    def unrestricted(self, val: bool):
        '''
        Whether the calculation should be unrestricted.
        '''
        self.settings.input.adf.Unrestricted = 'Yes' if val else 'No'

    def quality(self, val: str = 'Good'):
        '''
        Set the numerical quality of the calculation. Defaults to Good.

        Args:
            val: the numerical quality value to set to. This is the same as the ones used in the ADF GUI.
        '''
        self.settings.input.adf.NumericalQuality = val

    def SCF_converge(self, val: float = 1e-8):
        '''
        Set the convergence criterion for the SCF procedure. ADF default is 1e-6, but we use a stricter 1e-8 criterion.

        Args:
            val: the value of the criterion, default value = 1e-8.
        '''
        self.settings.input.adf.SCF.converge = val

    def functional(self, val: str):
        '''
        Set the functional to be used by the calculation. This also sets the dispersion if it is specified in the functional name.
        Note: setting the functional to r2SCAN-3c will automatically set the basis-set to mTZ2P.

        Args:
            val: the value to set the functional to. The value can be the same as the ones used in the ADF GUI. There are also extra functionals described in ADFJob.custom_disp_params. The keys of this dictionary are also valid functions. For a full list of functionals please see the ADFJob.available_functionals class-method.
        '''
        # before adding the new functional we should clear any previous functional settings
        self.settings.input.adf.pop('XC', None)

        functional = val.strip()
        functional = functional.replace('-D4(EEQ)', '-D4')  # D4(EEQ) and D4 are the same, unlike with D3 and D3(BJ)
        self._functional = functional

        if functional == 'r2SCAN-3c' and self._basis_set != 'mTZ2P':
            log.warn(f'Switching basis set from {self._basis_set} to mTZ2P for r2SCAN-3c.')
            self.basis_set('mTZ2P')

        if functional == 'SSB-D':
            log.error(f'There are two functionals called SSB-D, please use "GGA:SSB-D" or "MetaGGA:SSB-D".')
            return

        if functional not in functionals:
            log.warn(f'XC-functional {functional} not found. Please warn your local developer. Adding functional as LibXC.')
        else:
            func = functionals[functional]
            self.settings.input.adf.update(func.adf_settings)

        self.settings.input.adf.XC.LibXC = functional


    def relativity(self, level: str = 'Scalar'):
        '''
        Set the treatment of relativistic effects for this calculation. By default it is set to Scalar.

        Args:
            level: the level to set. Can be the same as the values in the ADF GUI.
        '''
        self.settings.input.adf.relativity.level = level

    def solvent(self, name: str = None, eps: float = None, rad: float = None, use_klamt: bool = False):
        if name:
            self._solvent = name
        else:
            self._solvent = f'COSMO(eps={eps} rad={rad})'

        if name == 'vacuum':
            self.settings.input.adf.pop('Solvation', None)
            return

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

    def molecule(self, mol: Union[str, plams.Molecule, plams.Atom, list[plams.Atom]]):
        '''
        Read a molecule in various formats.

        Args:
            mol: the molecule to read, can be a path (str). If the path exists already we read it. If it does not exist yet, it will be read in later. mol can also be a plams.Molecule object or a single or a list of plams.Atom objects.
        '''
        if isinstance(mol, plams.Molecule):
            self._molecule = mol

        elif isinstance(mol, str) and os.path.exists(mol):
            self._molecule = plams.Molecule(mol)
            # log.info(f'Succesfully loaded molecule {formula.molecule(self._molecule)} from path.')

        elif isinstance(mol, str):
            self._molecule = None
            self.settings.input.ams.system.GeometryFile = mol
            # log.info(f'Could not find molecule in file {mol}, will load it from the filename, so it should exist when the job starts.')

        elif isinstance(mol, list) and isinstance(mol[0], plams.Atom):
            self._molecule = plams.Molecule()
            [self._molecule.add_atom(atom) for atom in mol]

        elif isinstance(mol, plams.Atom):
            self._molecule = plams.Molecule()
            self._molecule.add_atom(mol)

    def run(self):
        '''
        Run this calculation. This method will write the input and run files and then calls sbatch with the provided settings.
        You can run this manually, or preferably use the context manager syntax for these jobs.
        '''
        os.makedirs(self.rundir, exist_ok=True)

        sett = self.settings.as_plams_settings()
        sett.keep = ['-', 't21.*', 't12.*', 'CreateAtoms.out', '$JN.dill']
        job = plams.AMSJob(name=self.name, molecule=self._molecule, settings=sett)

        os.makedirs(self.workdir, exist_ok=True)
        with open(j(self.workdir, f'{self.name}.in'), 'w+') as inpf:
            inpf.write(job.get_input())

        with open(j(self.workdir, f'{self.name}.run'), 'w+') as runf:
            runf.write('#!/bin/sh\n\n')
            runf.write('\n'.join(self._preambles) + '\n\n')
            runf.write(job.get_runscript())
            runf.write('\n'.join(self._postambles))

        cmd = self.get_sbatch_command() + f'-D {self.workdir} -J {self.rundir}/{self.name} -o ams.out {self.name}.run'
        with open(j(self.workdir, 'submit'), 'w+') as cmd_file:
            cmd_file.write(cmd)

        if not self.test_mode:
            with open(os.devnull, 'wb') as devnull:
                sp.run(cmd.split(), stdout=devnull, stderr=sp.STDOUT)

        # set the slurm job id for this calculation
        self.slurm_job_id = slurm.workdir_info(self.workdir).id
        # self.slurm_job_id = None

    @property
    def output_mol_path(self):
        return j(self.workdir, 'output.xyz')


class ADFFragmentJob(ADFJob):
    def __init__(self, *args, **kwargs):
        self.childjobs = {}
        super().__init__(*args, **kwargs)
        self.name = 'complex'

    def add_fragment(self, mol, name=None):
        if isinstance(mol, list) and isinstance(mol[0], plams.Atom):
            mol_ = plams.Molecule()
            [mol_.add_atom(atom) for atom in mol]
            mol = mol_

        name = name or f'fragment{len(self.childjobs) + 1}'
        self.childjobs[name] = ADFJob()
        self.childjobs[name].molecule(mol)
        setattr(self, name, self.childjobs[name])

        if self._molecule is None:
            self._molecule = self.childjobs[name]._molecule.copy()
        else:
            self._molecule = self._molecule + self.childjobs[name]._molecule.copy()

    def run(self):
        mol_str = " + ".join([formula.molecule(child._molecule) for child in self.childjobs.values()])
        log.flow(f'ADFFragmentJob [{mol_str}]', ['start'])
        # obtain some system wide properties of the molecules
        charge = sum([child.settings.input.ams.System.charge or 0 for child in self.childjobs.values()])
        unrestricted = any([(child.settings.input.adf.Unrestricted or 'no').lower() == 'yes' for child in self.childjobs.values()])
        spinpol = sum([child.settings.input.adf.SpinPolarization or 0 for child in self.childjobs.values()])
        log.flow(f'Level:             {self._functional}/{self._basis_set}')
        log.flow(f'Solvent:           {self._solvent}')
        log.flow(f'Charge:            {charge}', ['straight'])
        log.flow(f'Unrestricted:      {unrestricted}', ['straight'])
        log.flow(f'Spin-Polarization: {spinpol}', ['straight'])
        log.flow()
        # this job and all its children should have the same value for unrestricted
        [child.unrestricted(unrestricted) for child in self.childjobs.values()]

        # we now update the child settings with the parent settings
        # this is because we have to propagate settings such as functionals, basis sets etc.
        sett = self.settings.as_plams_settings()  # first create a plams settings object
        # same for the children
        child_setts = {name: child.settings.as_plams_settings() for name, child in self.childjobs.items()}
        # update the children using the parent settings
        [child_sett.update(sett) for child_sett in child_setts.values()]
        # same for sbatch settings
        [child.sbatch(**self._sbatch) for child in self.childjobs.values()]

        # now set the charge, spinpol, unrestricted for the parent 
        self.charge(charge)
        self.spin_polarization(spinpol)
        self.unrestricted(unrestricted)
        if unrestricted:
            self.settings.input.adf.UnrestrictedFragments = 'Yes'

        # now we are going to run each child job
        for i, (childname, child) in enumerate(self.childjobs.items(), start=1):
            log.flow(f'Child job ({i}/{len(self.childjobs)}) {childname} [{formula.molecule(child._molecule)}]', ['split'])
            log.flow(f'Charge:            {child.settings.input.ams.System.charge or 0}', ['straight', 'straight'])
            log.flow(f'Spin-Polarization: {child.settings.input.adf.SpinPolarization or 0}', ['straight', 'straight'])
            # the child name will be prepended with SP showing that it is the singlepoint calculation
            child.name = f'frag_{childname}'
            child.rundir = self.rundir

            # add the path to the child adf.rkf file as a dependency to the parent job
            self.settings.input.adf.fragments[childname] = j(child.workdir, 'adf.rkf')

            if child.can_skip():
                log.flow(log.Emojis.warning + f' Already ran, skipping', ['straight', 'end'])
                log.flow()
                # log.info(f'Child calculation {j(child.rundir, child.name)} already finished.')
                continue

            log.flow(log.Emojis.good + f' Submitting', ['straight', 'end'])
            # recast the plams.Settings object into a Result object as that is what run expects
            child.settings = results.Result(child_setts[childname])
            child.run()
            self.dependency(child)

            log.flow(f'SlurmID:  {child.slurm_job_id}', ['straight', 'skip', 'end'])
            log.flow(f'Work dir: {child.workdir}', ['straight', 'skip', 'end'])
            log.flow()

        # in the parent job the atoms should have the region and adf.f defined as options
        depend_atoms = []
        # for each atom we check which child it came from
        for atom in self._molecule:
            for childname, child in self.childjobs.items():
                for childatom in child._molecule:
                    # we check by looking at the symbol and coordinates of the atom
                    if (atom.symbol, atom.x, atom.y, atom.z) == (childatom.symbol, childatom.x, childatom.y, childatom.z):
                        # now write the symbol and coords as a string with the correct suffix
                        depend_atoms.append(f'\t\t{atom.symbol} {atom.x} {atom.y} {atom.z} region={childname} adf.f={childname}')
        # write the atoms block as a string with new line characters
        self.settings.input.ams.system.atoms = ('\n' + '\n'.join(depend_atoms) + '\n\tEnd').expandtabs(4)
        # set the _molecule to None, otherwise it will overwrite the atoms block
        self._molecule = None
        # run this job
        log.flow(log.Emojis.good + f' Submitting parent job', ['split'])
        super().run()
        log.flow(f'SlurmID: {self.slurm_job_id}', ['straight', 'end'])
        log.flow()

        # also do the calculation with SCF cycles set to 1
        self.settings.input.adf.SCF.Iterations = 1
        self.settings.input.adf.print = 'FMatSFO'  # by default print the fock matrix for each SCF cycle
        self.settings.input.adf.AllPoints = 'Yes'
        self.settings.input.adf.FullFock = 'Yes'
        self.name = self.name + '_SCF1'
        log.flow(log.Emojis.good + f' Submitting extra job with 1 SCF cycle', ['split'])

        super().run()
        log.flow(f'SlurmID: {self.slurm_job_id}', ['straight', 'end'])
        log.flow()
        log.flow(log.Emojis.finish + ' Done, bye!', ['startinv'])


class OrcaJob(Job):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
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


class NMRJob(Job):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pre_nmr_job = ADFJob()
        self.pre_nmr_job.functional('SAOP')
        self.pre_nmr_job.basis_set('TZ2P')
        self.pre_nmr_job.settings.input.adf.save = 'TAPE10'
        self.nics_points = []

    def get_runscript(self):
        preamble = '\n'.join(self._preambles)
        postamble = '\n'.join(self._postambles)
        return f'''#!/bin/sh
{preamble}
$AMSBIN/nmr {j(self.workdir, f"{self.name}.in")} < /dev/null
{postamble}
        '''

    def get_input(self):
        ghost_block = ''
        if len(self.nics_points) > 0:
            ghost_block += 'Ghosts\n'

        for point in self.nics_points:
            ghost_block += f'    {point[0]} {point[1]} {point[2]}\n'

        if len(self.nics_points) > 0:
            ghost_block += 'SubEnd\n'

        r = f'''NMR
    out all
    {ghost_block}
End
        '''

        return r

    def add_nics_point(self, p):
        self.nics_points.append(p)

    def molecule(self, *args, **kwargs):
        self.pre_nmr_job.molecule(*args, **kwargs)

    def run(self):
        os.makedirs(self.workdir, exist_ok=True)

        self.pre_nmr_job.name = self.name + '_pre'
        self.pre_nmr_job.rundir = self.rundir
        self.pre_nmr_job.sbatch(**self._sbatch)
        self.pre_nmr_job.add_postamble(f'cp {j(self.pre_nmr_job.workdir, "adf.rkf")} {j(self.workdir, "TAPE21")}')
        self.pre_nmr_job.add_postamble(f'cp {j(self.pre_nmr_job.workdir, "TAPE10")} {j(self.workdir, "TAPE10")}')
        self.pre_nmr_job.run()

        self.dependency(self.pre_nmr_job)
        self.add_postamble(f'mv {j(self.workdir, "TAPE21")} {j(self.workdir, "adf.rkf")}')

        with open(j(self.workdir, f'{self.name}.in'), 'w+') as inpf:
            inpf.write(self.get_input())

        with open(j(self.workdir, f'{self.name}.run'), 'w+') as runf:
            runf.write(self.get_runscript())

        cmd = self.get_sbatch_command() + f'-D {self.workdir} -J {self.rundir}/{self.name} -o ams.out {self.name}.run'
        with open(j(self.workdir, 'submit'), 'w+') as cmd_file:
            cmd_file.write(cmd)

        if not self.test_mode:
            with open(os.devnull, 'wb') as devnull:
                sp.run(cmd.split(), stdout=devnull, stderr=sp.STDOUT)

        # set the slurm job id for this calculation
        self.slurm_job_id = slurm.workdir_info(self.workdir).id


if __name__ == '__main__':
    for i, func in enumerate(get_available_functionals()):
        try:
            with ADFJob() as job:
                job.molecule('./test/xyz/H2O.xyz')
                job.rundir = 'tmp/functional_test'
                job.name = f'{i}.{func}'
                job.sbatch(p='tc', ntasks_per_node=32)
                # job.optimization()
                job.functional(func)
                job.basis_set('TZ2P')
                job.add_preamble('module load ams/2023.101')
        except Exception as e:
            print(e)

    # with NMRJob() as job:
    #     job.molecule('./test/xyz/H2O.xyz')
    #     job.rundir = 'tmp/NMR'
    #     job.name = 'H2O'
    #     job.sbatch(p='tc', ntasks_per_node=15)


    # with ADFFragmentJob() as job:
    #     mol = plams.Molecule('./test/xyz/NH3BH3.xyz')
    #     job.rundir = 'tmp/NH3BH3/EDA'
    #     job.sbatch(p='tc', ntasks_per_node=15)
    #     job.functional('r2SCAN')
    #     job.basis_set('TZ2P')
    #     job.add_fragment(mol.atoms[:4], 'Donor')
    #     job.add_fragment(mol.atoms[4:], 'Acceptor')

    # with ADFFragmentJob() as job:
    #     mol = plams.Molecule('./test/xyz/propyne.xyz')
    #     job.rundir = 'tmp/propyne/EDA'
    #     job.sbatch(p='tc', ntasks_per_node=15)
    #     job.functional('BP86')
    #     job.basis_set('TZ2P')
    #     job.quality('VeryGood')
    #     job.add_fragment(mol.atoms[:-1], 'Alkyne')
    #     job.add_fragment(mol.atoms[-1], 'Hydrogen')
    #     job.Alkyne.spin_polarization(-1)
    #     job.Hydrogen.spin_polarization(1)

    # with ADFFragmentJob() as job:
    #     mol = plams.Molecule('./test/xyz/SN2_TS.xyz')
    #     job.add_fragment(mol.atoms[:8], 'EtCl')
    #     job.add_fragment(mol.atoms[8:], 'Phenolate')
    #     job.Phenolate.charge(-1)

    #     job.rundir = 'tmp/SN2/EDA'
    #     job.sbatch(p='tc', ntasks_per_node=15)
    #     job.functional('OLYP')
    #     job.basis_set('DZP')

    # with ADFFragmentJob() as job:
    #     mol = plams.Molecule('./test/xyz/radadd.xyz')
    #     job.add_fragment(mol.atoms[:15], 'Substrate')
    #     job.add_fragment(mol.atoms[15:], 'Radical')
    #     job.Radical.spin_polarization(1)

    #     job.rundir = 'tmp/RA'
    #     job.sbatch(p='tc', ntasks_per_node=15)
    #     job.functional('BLYP-D3(BJ)')
    #     job.basis_set('TZ2P')

    # with ADFFragmentJob() as job:
    #     mol = plams.Molecule('./test/xyz/NaCl.xyz')
    #     job.add_fragment(mol.atoms[0], 'Cl')
    #     job.add_fragment(mol.atoms[1], 'Na')
    #     job.Na.charge(1)
    #     job.Cl.charge(-1)

    #     job.rundir = 'tmp/NaCl'
    #     job.sbatch(p='tc', ntasks_per_node=15)
    #     job.functional('BLYP-D3(BJ)')
    #     job.basis_set('TZ2P')
    #     job.quality('Good')

    with ADFJob() as opt_job:
        opt_job.molecule('./test/xyz/SN2_TS.xyz')
        opt_job.charge(-1)

        opt_job.rundir = 'tmp/SN2'
        opt_job.name = 'TS_OPT'
        opt_job.sbatch(p='tc', ntasks_per_node=15)
        opt_job.functional('BMK-D3(BJ)')
        opt_job.basis_set('DZP')
        opt_job.optimization()
        exit()

    # with ADFJob() as sp_job:
    #     sp_job.dependency(opt_job)  # this job will only run when opt_job finishes
    #     sp_job.molecule(opt_job.output_mol_path)  # we can take the output.xyz from the opt_job's workdir
    #     sp_job.charge(-1)

    #     sp_job.rundir = 'tmp/SN2'
    #     sp_job.name = 'SP_M062X'
    #     sp_job.sbatch(p='tc', ntasks_per_node=15)
    #     # now we can use a higher level of theory
    #     sp_job.functional('M06-2X')
    #     sp_job.basis_set('TZ2P')
