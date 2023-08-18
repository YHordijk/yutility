import scm.plams as plams
from yutility import units, log, settings, pathfunc, plotfunc
import os
import numpy as np
import subprocess
import shutil
import matplotlib.pyplot as plt

j = os.path.join

DEFAULT_RUN_PATH = './tmp'


def init(path, folder):
    plams.init(path=path, folder=folder)


def workdir():
    return plams.config.default_jobmanager.workdir


def run(mol, sett, name='calc', folder=None, path=DEFAULT_RUN_PATH, do_init=True, skip_already_run=False, run_kwargs={}):
    with log.NoPrint():
        if skip_already_run and os.path.exists(j(path, folder)):
            log.warn(f'Calculation in\n{j(path, folder)}\nalready run, skipping ...')
            return
        os.makedirs(path, exist_ok=True)
        if do_init:
            init(path, folder)

        job = plams.AMSJob(molecule=mol, settings=sett, name=name)
        results = job.run(**run_kwargs)
        if not check_success(job) is True:
            try:
                job.pickle()
            except BaseException:
                pass
        if do_init:
            plams.finish()
        return results


def optimize(mol, sett, name='optimize', folder=None, path=DEFAULT_RUN_PATH, do_init=True):
    with log.NoPrint():
        os.makedirs(path, exist_ok=True)
        if do_init:
            init(path, folder)

        job = plams.AMSJob(molecule=mol, settings=sett, name=name)
        results = job.run()
        if not check_success(job) is True:
            try:
                job.pickle()
            except BaseException:
                pass
        plams.finish()
        return results


def pes_scan(mol, sett, name='pes_scan', folder=None, path=DEFAULT_RUN_PATH, do_init=True):
    with log.NoPrint():
        os.makedirs(path, exist_ok=True)
        if do_init:
            init(path, folder)
        job = plams.AMSJob(molecule=mol, settings=sett, name=name)
        job.run()
        if not check_success(job) is True:
            try:
                job.pickle()
            except BaseException:
                pass

        kf = plams.KFFile(os.path.join(job.path, 'ams.rkf'))
        tkf = plams.RKFTrajectoryFile(os.path.join(job.path, 'ams.rkf'))
        energies = kf.read('PESScan', 'PES')
        idx = kf.read('PESScan', 'HistoryIndices')

        ts = tkf.get_plamsmol().copy()
        product = ts.copy()

        ts_idx = idx[np.argmax(energies)]
        tkf.read_frame(ts_idx-1, molecule=ts)
        product_idx = idx[-1]
        tkf.read_frame(product_idx-1, molecule=product)
        plams.finish()
    
        return {'transitionstate': ts, 'product': product}


def crest(mol, folder=None, path=DEFAULT_RUN_PATH, frozen_idxs=None, crest_path='crest', xtb_path='xtb', spinpol=0, charge=0, temp=400, mdlen='x1'):
    with log.NoPrint():
        rundir = pathfunc.next_path_version(j(path, folder))
        os.makedirs(rundir, exist_ok=True)
        mol_path = j(rundir, 'coords.xyz')
        mol.write(mol_path)

        with open(j(rundir, 'run.sh'), 'w+', newline='') as runsh:
            mol_path = str(mol_path).replace('\\', '/').replace('D:', '/mnt/d')
            runsh.write(f'cd {str(os.path.dirname(mol_path))}\n')
            if frozen_idxs is not None:
                runsh.write(f'{crest_path} coords.xyz -constrain {",".join([str(i) for i in frozen_idxs])}\n')
                constr = '-cinp .xcontrol.sample'
            else:
                constr = ''
            runsh.write(f'{crest_path} coords.xyz -xnam "{xtb_path}" --noreftopo -rthr 1 -c {charge} -u {spinpol} {constr} -tnmd {temp} -mdlen {mdlen}\n')
        runshp = j(rundir, 'run.sh').replace('\\', '/').replace('D:', '/mnt/d')
        with open(j(rundir, 'crest.log'), 'w+', newline='') as outfile:
            subprocess.call(['bash', f'{runshp}'], stdout=outfile)

        # Read in the conformers that were generated
        with open(j(rundir, 'crest_conformers.xyz')) as confs:
            lines = [line.strip() for line in confs.readlines()]

        # separate conformers into individual molecules
        natoms = int(lines[0])
        nmols = len(lines)//(natoms+2)
        mols = []
        for i in range(nmols):
            a, b = i*(natoms+2), (i+1)*(natoms+2), 
            mollines = lines[a:b]
            with open(j(rundir, 'tmp.xyz'), 'w+') as tmp:
                [tmp.write(line + '\n') for line in mollines]
            mol = plams.Molecule(j(rundir, 'tmp.xyz'))
            mols.append(mol)
        return mols


def get_calc_time(job):
    kf = plams.KFFile(os.path.join(job.path, 'ams.rkf'))
    time = kf.read('General', 'ElapsedTime')
    conv = units.convert_time(time)
    time, unit = time * conv[1], conv[0]
    return time, unit


def check_success(job):
    errf = j(job.path, f'{job.name}.err')
    with open(errf) as err:
        lines = [l.strip() for l in err.readlines()]
        if any(errtype in lines for errtype in ['ERROR: Geometry optimization failed! (Did not converge.)']):
            return False
        if any(successtype in lines for successtype in ['NORMAL TERMINATION', 'NORMAL TERMINATION with warnings']):
            return True
        else:
            return [l for l in lines if l.startswith('ERROR: ')]


def get_job_status(path):
    direc, name = os.path.split(path)
    if not os.path.exists(path):
        return 'notfound', None

    logp = j(path, 'ams.log')
    xyzp = j(path, 'output.xyz')
    reason = None
    success = False
    if os.path.exists(logp):
        with open(logp) as err:
            lines = [' '.join(l.strip().split()[2:]) for l in err.readlines()]
            if any(errtype in lines for errtype in ['ERROR: Geometry optimization failed! (Did not converge.)']):
                success = False
                reason = 'Geometry optimization failed.'
            elif any(successtype in lines for successtype in ['NORMAL TERMINATION', 'NORMAL TERMINATION with warnings']):
                success = True
    else:
        success = False

    if success and os.path.exists(xyzp):
        return 'completed', reason
    elif not success and os.path.exists(xyzp):
        return 'failed', reason
    else:
        return 'canceled', reason


class NMRResults:
    def __init__(self, kffile):
        self.kf = plams.KFReader(kffile)
        self.default_standards = {
            'H': 31.7,
            'C': 181.1,
        }

    def read(self, section, variable):
        return self.kf.read(section, variable)

    @property
    def elements(self):
        return self.read('Geometry', 'atomtype').split()

    @property
    def ordered_shifts(self):
        '''
        Returns chemical shifts following the input order of atoms
        '''
        return self.read('Properties', 'NMR Shieldings InputOrder')

    @property
    def chemical_shifts(self):
        peaks = {}
        for element in self.elements:
            peaks[element] = self.chemical_shift_by_element(element)
        return peaks 

    def chemical_shift_by_element(self, element):
        return self.chemical_shift_by_element_and_atom(element).values()

    def chemical_shift_by_element_and_atom(self, element):
        nnuc = self.read('Geometry', 'nnuc')
        atom_order_index = self.read('Geometry', 'atom order index')[nnuc:]
        element_order_index = self.read('Geometry', 'fragment and atomtype index')[nnuc:]
        peaks = self.read('Properties', 'NMR Shieldings InputOrder')

        elements_ordered = [self.elements[element_order_index[atom_order_index[i]-1]-1] for i in range(nnuc)]
        return {i: peaks[i] for i in range(nnuc) if elements_ordered[i] == element}

    # method to evaluate the NMR spectrum at given chemical shift(s)
    def __call__(self, *args, **kwargs):
        return self.evaluate(*args, **kwargs)

    def evaluate(self, shifts=None, element='H', width=0.1, standard=None, lineshape='lorentz', **kwargs):
        if element not in self.elements:
            return None, None

        lineshape_func = {
            'lorentz': lambda peak: 1/(1 + 4*(shifts - peak)**2/width**2),
            'gaussian': lambda peak: np.exp(-.6931471806 * (4*(shifts - peak)**2/width**2))
        }

        # default reference standard chemical shifts for some nuclei
        if standard is None:
            # default reference chemical shift is 0
            standard = self.default_standards.get(element, 0)

        peaks = [standard - peak for peak in self.chemical_shift_by_element(element)]

        # if shifts were not given we will find them ourselves
        if shifts is None:
            # we find the lowest and highest and offset them to prevent errors for elements with only 1 atom
            low = min(peaks) - width * 2
            high = max(peaks) + width * 2
            # move window 10% to low and high field and create a large space
            # partitioning could be improved by sampling more points close to peaks
            shifts = np.linspace(low - (high - low)*.1, high + (high - low)*.1, 10000)

        shifts = np.asarray(shifts)

        # iteratively add lines to spectrum
        spectrum = np.zeros_like(shifts)
        for peak in peaks:
            spectrum += lineshape_func[lineshape](peak)

        return shifts, spectrum

    def draw_spectrum(self, *args, **kwargs):
        x, y = self.evaluate(*args, **kwargs)
        if x is None:
            return

        element = kwargs.pop('element', 'H')
        kwargs.pop('shifts', None)
        kwargs.pop('width', None)
        standard = kwargs.pop('standard', self.default_standards.get(element, 0))
        kwargs.pop('lineshape', None)
        plt.plot(x, y, **kwargs)
        annotations = []
        for atom, peak in self.chemical_shift_by_element_and_atom(element).items():
            plt.vlines(standard - peak, 0, max(y) * 1.1, colors='grey', linestyles='dashed')
            annot = plt.gca().annotate(f'{element}{atom+1}', xycoords='data', va='top', ha='left', xy=(standard - peak, self.evaluate(standard - peak)[1]), textcoords='offset pixels', xytext=(10, 50))
            annotations.append(annot)
        plt.xlabel(rf'{element}-NMR $\delta$ (PPM)')
        plt.ylabel('Intensity')
        plotfunc.auto_texts(annotations)


def nmr(mol, dft_settings=None, folder=None, path=DEFAULT_RUN_PATH, do_init=True):
    with log.NoPrint():
        os.makedirs(path, exist_ok=True)
        if do_init:
            init(path, folder)

        # first run a dft calculation
        if dft_settings is None:
            # SAOP is supposed to be the best model potential for NMR
            dft_settings = settings.default('SAOP/TZ2P/Good')
        # we require TAPE10 file from dft
        dft_settings.input.adf.save = 'TAPE10'

        job = plams.AMSJob(molecule=mol, settings=dft_settings, name='pre_nmr')
        job.run()

        # copy files to new folder, NMR jobs need both TAPE21 and TAPE10 in the rundirectory
        # in this case, TAPE10 was written by the DFT job and TAPE21 is the adf.rkf file
        os.makedirs(j(workdir(), 'nmr'), exist_ok=True)
        shutil.copy2(j(workdir(), 'pre_nmr', 'adf.rkf'), j(workdir(), 'nmr', 'TAPE21'))
        shutil.copy2(j(workdir(), 'pre_nmr', 'TAPE10'), j(workdir(), 'nmr', 'TAPE10'))

        # generate a runscript
        runshp = j(workdir(), 'nmr', 'nmr.run')
        with open(runshp, 'w+') as infile:
            infile.write(f'cd {j(workdir(), "nmr")}\n')
            infile.write('"$AMSBIN/nmr" << eor\n')
            infile.write('NMR\n')
            infile.write('    out\n')
            infile.write('    Atoms ' + ' '.join([str(i+1) for i in range(len(mol.atoms))]) + '\n')
            infile.write('    u1k best\n')
            infile.write('    calc all\n')
            infile.write('End\n')
            infile.write('eor\n')
            infile.write('mv TAPE21 adf.rkf\n')

        # run the job
        with open(j(workdir(), 'nmr', 'nmr.out'), 'w+', newline='') as outfile:
            subprocess.call(['bash', f'{runshp}'], stdout=outfile)

        plams.finish()

        return NMRResults(j(workdir(), 'nmr', 'adf.rkf'))


def orbital_cub(rkf_path, name, orbtype='SCF', symlabel='All', overwrite=False):
    from yutility import volume
    
    def get_cub_file():
        for file in os.listdir(workdir):
            if file.startswith(cubprefix) and file.endswith('.cub'):
                return j(workdir, file)

    workdir = os.path.split(rkf_path)[0]
    cubprefix = str(name)

    # first check if the job is already done
    if not overwrite:
        if get_cub_file():
            return volume.CubeFile(get_cub_file())

    with open(j(os.path.split(rkf_path)[0], 'densf.in'), 'w+') as infile:
        infile.write(f'cd {os.path.split(rkf_path)[0]}\n')
        infile.write('"$AMSBIN/densf" << eor\n')
        infile.write(f'ADFFile {os.path.split(rkf_path)[1]}\n')
        infile.write(f'Orbitals {orbtype}\n')
        infile.write(f'    {symlabel} {name}\n')
        infile.write('END\n')
        infile.write(f'CUBOUTPUT {cubprefix}\n')
        infile.write('eor\n')

    with open(j(os.path.split(rkf_path)[0], 'densf.out'), 'w+', newline='') as outfile:
        subprocess.call(['bash', f'{j(os.path.split(rkf_path)[0], "densf.in")}'], stdout=outfile)

    return volume.CubeFile(get_cub_file())


# class ADFFragmentResults(plams.ADFFragmentResults):
#     def __init__(self, **kwargs):
#         super().__init__(**kwargs)



class ADFFragmentJob(plams.MultiJob):
    _result_type = plams.ADFFragmentResults

    def __init__(self, fragment1=None, fragment2=None, full_settings=None, frag1_settings=None, frag2_settings=None, fragment_names=None, **kwargs):
        plams.MultiJob.__init__(self, **kwargs)
        self.fragment1 = fragment1.copy() if isinstance(fragment1, plams.Molecule) else fragment1
        self.fragment2 = fragment2.copy() if isinstance(fragment2, plams.Molecule) else fragment2
        self.fragment_names = fragment_names or ['frag1', 'frag2']
        self.full_settings = full_settings or plams.Settings()
        self.frag1_settings = frag1_settings or plams.Settings()
        self.frag2_settings = frag2_settings or plams.Settings()

    def prerun(self):
        self.f1 = plams.AMSJob(name=self.fragment_names[0], molecule=self.fragment1, settings=self.settings + self.frag1_settings)
        self.f2 = plams.AMSJob(name=self.fragment_names[1], molecule=self.fragment2, settings=self.settings + self.frag2_settings)

        for at in self.fragment1:
            at.properties.suffix = f'adf.f={self.fragment_names[0]}'
        for at in self.fragment2:
            at.properties.suffix = f'adf.f={self.fragment_names[1]}'

        self.full = plams.AMSJob(name='full',
                                 molecule=self.fragment1 + self.fragment2,
                                 settings=self.settings + self.full_settings)

        self.full.settings.input.adf.fragments[self.fragment_names[0]] = (self.f1, 'adf')
        self.full.settings.input.adf.fragments[self.fragment_names[1]] = (self.f2, 'adf')

        self.children = [self.f1, self.f2, self.full]


def EDA(mol1, mol2, settings=None, full_settings=None, frag1_settings=None, frag2_settings=None, fragment_names=None, folder=None, path=DEFAULT_RUN_PATH, do_init=True):
    with log.NoPrint():
        os.makedirs(path, exist_ok=True)
        if do_init:
            init(path, folder)

        job = ADFFragmentJob(settings=settings, fragment1=mol1, fragment2=mol2, 
                             full_settings=full_settings, frag1_settings=frag1_settings, 
                             frag2_settings=frag2_settings, name='EDA', fragment_names=fragment_names)
        result = job.run()

        return result


if __name__ == '__main__':
    from ychem.results import reaction2
    from ReactionRunner import generators

    class SubstrateCatalyst(generators.generic):
        reaction_name = 'SubstrateCatalyst'
        expected_stationary_points = ['substrate', 'catalyst', 'substrate_cat_complex']

        def _generate(self):
            self.define_sp('substrate_cat_complex', ['substrate', 'catalyst'])


            if self.check_hash and self._hash_collision(self.expected_stationary_points):
                for sp in self.expected_stationary_points:
                    self.load(sp)

            else:
                self.combine_from_sphere('substrate_cat_complex', 'substrate', 'catalyst', radius=5, frozen_bonds=True, ntries=20,
                                         description='Creating Substrate-Catalyst-Complex')

            self.save()
            if self.delete_tmp:
                self.clean_tmp()


    # d = '/Users/yumanhordijk/PhD/ychem/calculations2/0b1794d72ee3b1eed65d7c6e50cf9deb7ff567a663d19e27675df55a084bf3a3'
    # mol = plams.Molecule(j(d, 'substrate', 'geometry', 'output.xyz'))
    # substrate_spectrum = nmr(mol, dft_settings=settings.default('Cheap'), path=j(d, 'substrate'), folder='NMR')
    # print(substrate_spectrum.chemical_shifts)
    # print(chemical_shift_by_element)

    # mol = plams.Molecule(j(d, 'substrate_cat_complex', 'geometry', 'output.xyz'))
    # substrate_catalyst_spectrum = nmr(mol, dft_settings=settings.default('Cheap'), path=j(d, 'substrate_cat_complex'), folder='NMR')
    # print(substrate_catalyst_spectrum.chemical_shifts)
    
    mol = plams.Molecule('acroleine.xyz')
    substrate_spectrum = nmr(mol, dft_settings=settings.default('SAOP/TZ2P/Good'), path=j(d, 'substrate'), folder='NMR')
    substrate_cat_spectra = {}

    for cat in ['I2', 'SnCl4', 'ZnCl2', 'TiCl4', 'BF3', 'AlCl3']:
        reactants = {
            'substrate': 'acrolein_E',
            'catalyst': 'catalyst',
        }
        substituents = {
            'substrate': {
                'R1': 'Me',
                'R2': 'Me'
            },
            'catalyst': {
                'Rcat': cat
            }
        }
        rg = SubstrateCatalyst(reactants=reactants, substituents=substituents)
        mols = rg.generate()
        mol = mols['substrate_cat_complex']
        substrate_cat_spectra[cat] = nmr(mol, dft_settings=settings.default('Cheap'), path=j(d, 'substrate'), folder='NMR')

    rxn = reaction2.Reaction(main_dir='/Users/yumanhordijk/PhD/ychem/calculations2/0b1794d72ee3b1eed65d7c6e50cf9deb7ff567a663d19e27675df55a084bf3a3')
    substrate_spectrum = NMRResults(j(rxn.main_dir, 'substrate', 'NMR', 'nmr', 'adf.rkf'))
    substrate_catalyst_spectrum = NMRResults(j(rxn.main_dir, 'substrate_cat_complex', 'NMR', 'nmr', 'adf.rkf'))
    substrate_spectrum.draw_spectrum(element='H', label='Substrate')
    substrate_catalyst_spectrum.draw_spectrum(element='H', label='Substrate-Catalyst-Complex')
    plt.legend()
    plt.show()
