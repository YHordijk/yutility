from yutility import settings, run, molecule
from TCutility import results
import os 
from scm import plams

j = os.path.join


def run_full(xyzfile):
    mol = molecule.load(j(xyzdir, xyzfile))
    preset = mol['flags'].get('level', 'OLYP/TZ2P/Good')
    sett = settings.default(preset)
    sett = settings.charge(mol['flags']['charge'], sett)
    sett = settings.spin_polarization(mol['flags']['spinpol'], sett)
    sett = settings.vibrations(sett)
    if 'solvent' in mol['flags']:
        sett = settings.solvation(mol['flags']['solvent'], sett)

    run.run(mol['molecule'], sett, path=j('calculations', xyzfile.removesuffix('.xyz')), folder='dft_freq', skip_already_run=True)


def run_quick(xyzfile, k=3):
    mol = molecule.load(j(xyzdir, xyzfile))
    preset = mol['flags'].get('level', 'OLYP/TZ2P/Good')
    sett = settings.default(preset)
    sett = settings.charge(mol['flags']['charge'], sett)
    sett = settings.spin_polarization(mol['flags']['spinpol'], sett)
    sett = settings.vibrations(sett)
    if 'solvent' in mol['flags']:
        sett = settings.solvation(mol['flags']['solvent'], sett)

    run.quick_SP_check(mol['molecule'], sett, k=k, path=j('calculations', xyzfile.removesuffix('.xyz')), folder=f'quick_freq_{k}', skip_already_run=True)


def run_ppc(xyzfile):
    mol = molecule.load(j(xyzdir, xyzfile))
    preset = mol['flags'].get('level', 'OLYP/TZ2P/Good')
    sett = settings.default(preset)
    sett = settings.charge(mol['flags']['charge'], sett)
    sett = settings.spin_polarization(mol['flags']['spinpol'], sett)
    if 'solvent' in mol['flags']:
        sett = settings.solvation(mol['flags']['solvent'], sett)
    sett.input.ams.properties.PESPointCharacter = 'Yes'

    run.run(mol['molecule'], sett, path=j('calculations', xyzfile.removesuffix('.xyz')), folder='PPC', skip_already_run=True)
    

xyzdir = j('..', 'xyz')
for xyzfile in os.listdir(xyzdir):
    if xyzfile == '.DS_Store':
        continue

    # if xyzfile != '1C.xyz':
    #     continue
        
    if os.path.exists(j('calculations', xyzfile.removesuffix('.xyz'))):
        continue

    os.makedirs(j('calculations', xyzfile.removesuffix('.xyz')), exist_ok=True)
    run_quick(xyzfile, 1)
    run_quick(xyzfile, 2)
    run_quick(xyzfile, 3)
    run_full(xyzfile)
    run_ppc(xyzfile)
