from yutility import settings, run, molecule
from TCutility import results
import os 
from scm import plams

j = os.path.join

timing_ratios = []
errors = []
natoms = []
xyzdir = j('..', 'xyz')
for xyzfile in os.listdir(xyzdir):
    if xyzfile == '.DS_Store':
        continue

    # initialize the calculations
    os.makedirs('calculations', exist_ok=True)
    # plams.init(path='calculations', folder=xyzfile.removesuffix('.xyz') + '_PPC')

    # load the molecule. The comment line contains the expected type of SP, either minimum or transitionstate
    mol = molecule.load(j(xyzdir, xyzfile))
    preset = mol['flags'].get('level', 'OLYP/TZ2P/Good')

    # first calculate the full vibrational spectrum at the dft level
    dft_settings = settings.default(preset)
    dft_settings = settings.charge(mol['flags']['charge'], dft_settings)
    dft_settings = settings.spin_polarization(mol['flags']['spinpol'], dft_settings)
    dft_settings = settings.vibrations(dft_settings)
    run.run(mol['molecule'], dft_settings, name='dft_freq', path='calculations', folder=xyzfile.removesuffix('.xyz'), skip_already_run=True)
    
    # then perform the quick_SP_check
    dft_settings = settings.default(preset)
    dft_settings = settings.charge(mol['flags']['charge'], dft_settings)
    dft_settings = settings.spin_polarization(mol['flags']['spinpol'], dft_settings)
    run.quick_SP_check(mol['molecule'], dft_settings, path='calculations', folder=xyzfile.removesuffix('.xyz'), skip_already_run=True)

    # then perform the pespointcharacter from ADF
    dft_settings = settings.default(preset)
    dft_settings = settings.charge(mol['flags']['charge'], dft_settings)
    dft_settings = settings.spin_polarization(mol['flags']['spinpol'], dft_settings)
    dft_settings.input.ams.properties.PESPointCharacter = 'Yes'
    run.run(mol['molecule'], dft_settings, path='calculations', folder=xyzfile.removesuffix('.xyz') + '_PPC', skip_already_run=True)
    # refine_res = results.read(j('calculations', xyzfile.removesuffix('.xyz'), 'refine'))
    break

    # plams.finish()

#   # analyse the performance boost and errors
#   timing_ratio = opt_res.timing.total/refine_res.timing.total
#   timing_ratios.append(timing_ratio)

#   natoms.append(len(mol['molecule']))
#   print(opt_res.properties.vibrations.frequencies[0], refine_res.properties.vibrations.frequencies[0])
#   errors.append(opt_res.properties.vibrations.frequencies[0] - refine_res.properties.vibrations.frequencies[0])

# print(timing_ratios)
# print(natoms)
# print(errors)