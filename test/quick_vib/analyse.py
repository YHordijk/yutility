from yutility import settings, run, molecule
from TCutility import results
from ychem.ML import stats
import os 
from scm import plams

j = os.path.join


xyzdir = j('..', 'xyz')
# for xyzfile in os.listdir(xyzdir):
#     if xyzfile == '.DS_Store':
#         continue
#     # initialize the calculations
#     os.makedirs('calculations', exist_ok=True)
#     plams.init(path='calculations', folder=xyzfile.removesuffix('.xyz'))

#     # load the molecule. The comment line contains the expected type of SP, either minimum or transitionstate
#     mol = molecule.load(j(xyzdir, xyzfile))
#     # first optimize and calculate the vibrations at the dft level
#     # dft_settings = settings.default('OLYP/TZ2P/Good')
#     # if mol['comment'] == 'minimum':
#     #     dft_settings = settings.optimization(dft_settings)
#     # else:
#     #     dft_settings = settings.transition_state(sett=dft_settings)

#     # dft_settings = settings.charge(mol['flags']['charge'], dft_settings)
#     # run.run(mol['molecule'], dft_settings, name='dft_opt', path='calculations', folder=xyzfile.removesuffix('.xyz'), skip_already_run=True)
    
#     dft_settings = settings.default('OLYP/TZ2P/Good')
#     dft_settings = settings.charge(mol['flags']['charge'], dft_settings)
#     dft_settings = settings.vibrations(dft_settings)
#     run.run(mol['molecule'], dft_settings, name='dft_freq', path='calculations', folder=xyzfile.removesuffix('.xyz'), skip_already_run=True)
#     # fullfreq_res = results.read(j('calculations', xyzfile.removesuffix('.xyz'), 'dft_opt'))

#     # then perform the quick_SP_check
#     dft_settings = settings.default('OLYP/TZ2P/Good')
#     dft_settings = settings.charge(mol['flags']['charge'], dft_settings)
#     run.quick_SP_check(mol['molecule'], dft_settings, path='calculations', folder=xyzfile.removesuffix('.xyz'), skip_already_run=True, do_init=False)
#     # refine_res = results.read(j('calculations', xyzfile.removesuffix('.xyz'), 'refine'))

#     plams.finish()


import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

nfreqs = 1
timing_ratios = []
timing_refine = []
timing_full = []
timing_PPC = []
errors = []
natoms = []
nsuccess = 0
nfailed = 0
character = []
freqs_ref = []
freqs = []
character_pred = []
modes_ref = []
modes = []
for xyzfile in os.listdir(xyzdir):
    if xyzfile == '.DS_Store':
        continue
    mol = molecule.load(j(xyzdir, xyzfile))
    if not os.path.exists(j('calculations', xyzfile.removesuffix('.xyz'))):
        continue
    if not os.path.exists(j('calculations', xyzfile.removesuffix('.xyz') + '_PPC')):
        continue
    fullfreq_res = results.read(j('calculations', xyzfile.removesuffix('.xyz'), 'dft_freq'))
    # refine_res = results.read(j('calculations', xyzfile.removesuffix('.xyz') + '_2', 'refine'))
    PPC_res = results.read(j('calculations', xyzfile.removesuffix('.xyz') + '_PPC', 'calc'))
    refine_res = results.read(j('calculations', xyzfile.removesuffix('.xyz'), 'refine'))

    if not refine_res.properties.vibrations.frequencies:
        nfailed += 1
        continue

    if fullfreq_res.properties.vibrations.character == refine_res.properties.vibrations.character:
        nsuccess += 1
    else:
        nfailed += 1

    character.append(fullfreq_res.properties.vibrations.character)
    character_pred.append(refine_res.properties.vibrations.character)

    # analyse the performance boost and errors
    timing_ratio = fullfreq_res.timing.total/refine_res.timing.total
    timing_full.append(fullfreq_res.timing.total/60)
    timing_refine.append(refine_res.timing.total/60)
    timing_PPC.append(PPC_res.timing.total/60)
    timing_ratios.append(timing_ratio)


    freqs_ref.extend(fullfreq_res.properties.vibrations.frequencies[:nfreqs])
    freqs.extend(refine_res.properties.vibrations.frequencies[:nfreqs])
    modes_ref.extend(fullfreq_res.properties.vibrations.modes[:nfreqs])
    modes.extend(refine_res.properties.vibrations.modes[:nfreqs])
    natoms.append(len(mol['molecule']))
    errors.append(abs((fullfreq_res.properties.vibrations.frequencies[0] - refine_res.properties.vibrations.frequencies[0])/fullfreq_res.properties.vibrations.frequencies[0] * 100))

timing_ratios = np.array(timing_ratios)
print(timing_ratios)
natoms = np.array(natoms)
errors = np.array(errors)

print(f'Success rate: {nsuccess/(nsuccess + nfailed):.0%}')

curve = lambda x, n, c, a: a*(x+c)**n
plt.figure()
plt.title('Total Calculation Time for Both Methods\n' + r'Fitted to $t(N_{atom}) = a(N_{atom}+c)^n$')
popt = curve_fit(curve, natoms, timing_refine, p0=[2, 0, 7e-4])[0]
plt.plot(np.linspace(min(natoms), max(natoms), 100), curve(np.linspace(min(natoms), max(natoms), 100), *popt), linestyle='dashed')
plt.scatter(natoms, timing_refine, label=rf'$t_{{Quick}} = {popt[2]:.1E}(N_{{atom}}{popt[1]:+.1f})^{{{popt[0]:.1f}}}$')

popt = curve_fit(curve, natoms, timing_full, p0=[3, 0, 5e-6])[0]
plt.plot(np.linspace(min(natoms), max(natoms), 100), curve(np.linspace(min(natoms), max(natoms), 100), *popt), linestyle='dashed')
plt.scatter(natoms, timing_full, label=rf'$t_{{Full}} = {popt[2]:.1E}(N_{{atom}}{popt[1]:+.1f})^{{{popt[0]:.1f}}}$')

popt = curve_fit(curve, natoms, timing_PPC, p0=[3, 0, 5e-6])[0]
plt.plot(np.linspace(min(natoms), max(natoms), 100), curve(np.linspace(min(natoms), max(natoms), 100), *popt), linestyle='dashed')
plt.scatter(natoms, timing_PPC, label=rf'$t_{{PPC}} = {popt[2]:.1E}(N_{{atom}}{popt[1]:+.1f})^{{{popt[0]:.1f}}}$')

plt.xlabel('Number of atoms')
plt.ylabel(r'$t [min]$')
plt.legend()

plt.figure()
plt.title(r'Time Gain Ratio ($\gamma$)')
popt = curve_fit(lambda x, a, b: a*x + b, natoms, timing_ratios, p0=[1, 0])[0]
plt.plot(np.linspace(min(natoms), max(natoms), 2), popt[0] * np.linspace(min(natoms), max(natoms), 2) + popt[1], 
         linestyle='dashed', label=rf'$\gamma = {popt[0]:.1f}N_{{atom}} + {popt[1]:.1f}$')
plt.scatter(natoms, timing_ratios)
plt.xlabel('Number of atoms')
plt.ylabel(r'$\gamma$')
plt.legend()

plt.figure()
plt.scatter(freqs_ref, freqs)
plt.xlabel(r'$\bar{\nu}_{Full}$ [$cm^{-1}$]')
plt.ylabel(r'$\bar{\nu}_{Quick}$ [$cm^{-1}$]')

xlim = plt.xlim()
ylim = plt.ylim()

rmin, rmax = min(xlim[0], ylim[0]), max(xlim[1], ylim[1])
plt.plot([rmin, rmax], [rmin, rmax], linestyle='dashed', c='k')
plt.plot([0, 0], [rmin, rmax], c='k', alpha=.3)
plt.plot([rmin, rmax], [0, 0], c='k', alpha=.3)

xlim = plt.xlim(xlim)
ylim = plt.ylim(ylim)

plt.figure()
modes = [mode/sum(np.array(mode)**2) for mode in modes]
modes_ref = [mode/sum(np.array(mode)**2) for mode in modes_ref]

mode_errors = [stats.MAE(mode, mode_ref) for mode, mode_ref in zip(modes, modes_ref)]
natoms3 = []
[natoms3.extend([atom]*3) for atom in natoms]

plt.scatter([natom for natom, freq in zip(natoms3, freqs_ref) if freq >= 0], [err for err, freq in zip(mode_errors, freqs_ref) if freq >= 0])
plt.scatter([natom for natom, freq in zip(natoms3, freqs_ref) if freq < 0], [err for err, freq in zip(mode_errors, freqs_ref) if freq < 0])
plt.xlabel(r'Number of atoms')
plt.ylabel(r'Mode error (MAE)')

xlim = plt.xlim()
ylim = plt.ylim()

rmin, rmax = min(xlim[0], ylim[0]), max(xlim[1], ylim[1])
plt.plot([rmin, rmax], [rmin, rmax], linestyle='dashed', c='k')
plt.plot([0, 0], [rmin, rmax], c='k', alpha=.3)
plt.plot([rmin, rmax], [0, 0], c='k', alpha=.3)

xlim = plt.xlim(xlim)
ylim = plt.ylim(ylim)


# plt.figure()
# C = [[0, 0], [0, 0]]
# for i, char1 in enumerate(['minimum', 'transitionstate']):
#     for k, char2 in enumerate(['minimum', 'transitionstate']):
#         C[i][k] = sum([character[l] == char1 and character_pred[l] == char2 for l in range(len(character))])
# plt.imshow(C)
# plt.xlabel('True label')
# plt.ylabel('Predicted label')


plt.show()
