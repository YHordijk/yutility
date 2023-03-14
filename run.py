import scm.plams as plams
from yutility import units, config, log
import os
import sys
import numpy as np
import subprocess

j = os.path.join

DEFAULT_RUN_PATH = './tmp'


def init(path, folder):
    plams.init(path=path, folder=folder)


def workdir():
    return plams.config.default_jobmanager.workdir


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


def crest(mol, folder=None, path=DEFAULT_RUN_PATH, frozen_idxs=None, crest_path='crest', xtb_path='xtb', spinpol=0, charge=0, temp=298.15, mdlen='x1'):
    with log.NoPrint():
        os.makedirs(j(path, folder + '_crest'), exist_ok=True)
        mol_path = j(path, folder + '_crest', 'coords.xyz')
        mol.write(mol_path)

        with open(j(path, folder + '_crest', 'run.sh'), 'w+', newline='') as runsh:
            mol_path = str(mol_path).replace('\\', '/').replace('D:', '/mnt/d')
            runsh.write(f'cd {str(os.path.dirname(mol_path))}\n')
            if frozen_idxs is not None:
                runsh.write(f'{crest_path} coords.xyz -constrain {",".join([str(i) for i in frozen_idxs])}\n')
                constr = '-cinp .xcontrol.sample'
            else:
                constr = ''
            runsh.write(f'{crest_path} coords.xyz -xnam "{xtb_path}" --noreftopo -rthr 1 -c {charge} -u {spinpol} {constr} -tnmd {temp} -mdlen {mdlen}\n')
        runshp = j(path, folder + '_crest', 'run.sh').replace('\\', '/').replace('D:', '/mnt/d')
        with open(j(path, folder + '_crest', 'crest.log'), 'w+', newline='') as outfile:
            subprocess.call(['bash', f'{runshp}'], stdout=outfile)

        # Read in the conformers that were generated
        with open(j(path, folder + '_crest', 'crest_conformers.xyz')) as confs:
            lines = [line.strip() for line in confs.readlines()]

        # separate conformers into individual molecules
        natoms = int(lines[0])
        nmols = len(lines)//(natoms+2)
        mols = []
        for i in range(nmols):
            a, b = i*(natoms+2), (i+1)*(natoms+2), 
            mollines = lines[a:b]
            with open(j(path, folder + '_crest', 'tmp.xyz'), 'w+') as tmp:
                [tmp.write(line + '\n') for line in mollines]
            mol = plams.Molecule(j(path, folder + '_crest', 'tmp.xyz'))
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



# class ADFFragmentJob(plams.MultiJob):
#     _result_type = plams.ADFFragmentResults

#     def __init__(self, fragment1=None, fragment2=None, full_settings=None, frag1_settings=None, frag2_settings=None, fragment_names=None, **kwargs):
#         plams.MultiJob.__init__(self, **kwargs)
#         self.fragment1 = fragment1.copy() if isinstance(fragment1, plams.Molecule) else fragment1
#         self.fragment2 = fragment2.copy() if isinstance(fragment2, plams.Molecule) else fragment2
#         self.fragment_names = fragment_names or ['frag1', 'frag2']
#         self.full_settings = full_settings or plams.Settings()
#         self.frag1_settings = frag1_settings or plams.Settings()
#         self.frag2_settings = frag2_settings or plams.Settings()

#     def prerun(self):
#         self.f1 = plams.AMSJob(name=self.fragment_names[0], molecule=self.fragment1, settings=self.settings + self.frag1_settings)
#         self.f2 = plams.AMSJob(name=self.fragment_names[1], molecule=self.fragment2, settings=self.settings + self.frag2_settings)

#         for at in self.fragment1:
#             at.properties.suffix = 'adf.f=subsystem1'
#         for at in self.fragment2:
#             at.properties.suffix = 'adf.f=subsystem2'

#         self.full = plams.AMSJob(name='full',
#                                  molecule=self.fragment1 + self.fragment2,
#                                  settings=self.settings + self.full_settings)

#         self.full.settings.input.adf.fragments.subsystem1 = (self.f1, 'adf')
#         self.full.settings.input.adf.fragments.subsystem2 = (self.f2, 'adf')

#         self.children = [self.f1, self.f2, self.full]

# def fragment(mol1, mol2, sett)
