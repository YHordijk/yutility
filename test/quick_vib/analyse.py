from yutility import settings, run, molecule, dictfunc
from TCutility import results
from ychem.ML import stats
import os 
from scm import plams
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

j = os.path.join


xyzdir = j('..', 'xyz')

data = dictfunc.DotDict()

# read the raw data first
# we check every atom in the xyz directory
for xyzfile in os.listdir(xyzdir):
    if xyzfile == '.DS_Store':
        continue

    molname = xyzfile.removesuffix('.xyz')

    if not os.path.exists(j('calculations', molname)):
        continue

    # get the number of atoms
    mol = molecule.load(j(xyzdir, xyzfile))
    data.natoms[molname] = mol['natoms']

    # for each method (PPC, quick_1, full, etc.) we read the data
    for method in os.listdir(j('calculations', molname)):
        if method.startswith('quick_freq_'):
            calcdir = j('calculations', molname, method, 'refine')
        else:
            calcdir = j('calculations', molname, method, 'calc')

        # we save timing, frequencies and character
        res = results.read(calcdir)
        data.timing[method][molname] = res.timing.total/60  # in minutes
        data.frequencies[method][molname] = res.properties.vibrations.frequencies
        data.character[method][molname] = res.properties.vibrations.character

print(data.frequencies['quick_freq_3'])
# now we distill some information from the raw data
methods = set(data.timing.keys())
molnames = list(set([key for method in methods for key in data.timing[method].keys()]))
natoms = np.array([data.natoms[molname] for molname in molnames], dtype=float)

# first get the timing_ratios
for method in methods:
    data.timing_ratios[method] = []
    data.timing_sorted[method] = []

    for molname in molnames:
        full_method_timing = data.timing['dft_freq'][molname]
        timing = data.timing[method][molname]
        data.timing_ratios[method].append(full_method_timing/timing)
        data.timing_sorted[method].append(timing)


methods_fancy = {
    'dft_freq': r'Full',
    'quick_freq_1': r'Quick(1)',
    'quick_freq_2': r'Quick(2)',
    'quick_freq_3': r'Quick(3)',
    'PPC': r'PPC',
}

# get the accuracies
for method in methods:
    nsuccess = 0
    nfailed = 0
    for molname in molnames:
        full_method_character = data.character['dft_freq'][molname]
        character = data.character[method][molname]
        if character == full_method_character:
            nsuccess += 1
        else:
            nfailed += 1
    data.accuracy[method] = nsuccess / (nsuccess + nfailed)


curve = lambda x, n, c, a: a*(x+c)**n
plt.figure()
plt.title('Total Calculation Time\n' + r'Fitted to $t(N_{atom}) = a(N_{atom}+c)^n$')
for method in methods:
    print(natoms, data.timing_sorted[method])
    popt = curve_fit(curve, natoms, data.timing_sorted[method], p0=[2, 0, 7e-4])[0]
    plt.plot(np.linspace(min(natoms), max(natoms), 100), curve(np.linspace(min(natoms), max(natoms), 100), *popt), linestyle='dashed')
    plt.scatter(natoms, data.timing_sorted[method], label=rf'$t_{{{methods_fancy.get(method)}}} = {popt[2]:.1E}(N_{{atom}}{popt[1]:+.1f})^{{{popt[0]:.1f}}}$', s=10)
plt.xlabel('Number of atoms')
plt.ylabel(r'$t [min]$')
plt.legend()
plt.savefig('figures/timings.png')
plt.close()


plt.figure()
plt.title(r'Time Gain Ratio ($\gamma$)')
for method in methods:
    if method == 'dft_freq':
        continue

    popt = curve_fit(lambda x, a, b: a*x + b, natoms, data.timing_ratios[method], p0=[1, 0])[0]
    plt.plot(np.linspace(min(natoms), max(natoms), 2), popt[0] * np.linspace(min(natoms), max(natoms), 2) + popt[1], 
             linestyle='dashed', label=rf'$\gamma_{{{methods_fancy.get(method)}}} = {popt[0]:.1f}N_{{atom}} + {popt[1]:.1f}$')
    plt.scatter(natoms, data.timing_ratios[method])
plt.xlabel('Number of atoms')
plt.ylabel(r'$\gamma$')
plt.legend()
plt.savefig('figures/timing_ratios.png')
plt.close()


plt.figure()
freqs = np.array(freqs) 
freqs_ref = np.array(freqs_ref)
plt.scatter(freqs_ref, freqs)
plt.xlabel(r'$\bar{\nu}_{Full}$ [$cm^{-1}$]')
plt.ylabel(r'$\bar{\nu}_{Quick}$ [$cm^{-1}$]')
MRE = np.mean(abs(freqs - freqs_ref)/freqs_ref)
plt.title(f'Quick Method Predicted frequencies\nMRE = {MRE:.1%}, acc. = {success_rate:.1%}')

xlim = plt.xlim()
ylim = plt.ylim()

rmin, rmax = min(xlim[0], ylim[0]), max(xlim[1], ylim[1])
plt.plot([rmin, rmax], [rmin, rmax], linestyle='dashed', c='k')
plt.plot([0, 0], [rmin, rmax], c='k', alpha=.3)
plt.plot([rmin, rmax], [0, 0], c='k', alpha=.3)
plt.gca().add_patch(Rectangle((rmin, 0), abs(rmin), rmax, facecolor='red', alpha=.2))
plt.gca().add_patch(Rectangle((0, rmin), rmax, abs(rmin), facecolor='red', alpha=.2))

xlim = plt.xlim(xlim)
ylim = plt.ylim(ylim)


plt.figure()
freqs_DFTB = np.array(freqs_DFTB) 
freqs_DFTB_ref = np.array(freqs_DFTB_ref)
plt.scatter(freqs_DFTB, freqs_DFTB_ref)
plt.xlabel(r'$\bar{\nu}_{Full}$ [$cm^{-1}$]')
plt.ylabel(r'$\bar{\nu}_{DFTB}$ [$cm^{-1}$]')
MRE = np.mean(abs(freqs_DFTB - freqs_DFTB_ref)/freqs_DFTB_ref)
plt.title(f'DFTB Predicted frequencies\nMRE = {MRE:.1%}, acc. = {success_rate_DFTB:.1%}')

xlim = plt.xlim()
ylim = plt.ylim()

rmin, rmax = min(xlim[0], ylim[0]), max(xlim[1], ylim[1])
plt.plot([rmin, rmax], [rmin, rmax], linestyle='dashed', c='k')
plt.plot([0, 0], [rmin, rmax], c='k', alpha=.3)
plt.plot([rmin, rmax], [0, 0], c='k', alpha=.3)
plt.gca().add_patch(Rectangle((rmin, 0), abs(rmin), rmax, facecolor='red', alpha=.2))
plt.gca().add_patch(Rectangle((0, rmin), rmax, abs(rmin), facecolor='red', alpha=.2))

xlim = plt.xlim(xlim)
ylim = plt.ylim(ylim)


plt.figure()
freqs_PPC = np.array(freqs_PPC) 
freqs_PPC_ref = np.array(freqs_PPC_ref)
print(freqs_PPC.size, freqs_PPC_ref.size)
plt.scatter(freqs_PPC, freqs_PPC_ref)
plt.xlabel(r'$\bar{\nu}_{Full}$ [$cm^{-1}$]')
plt.ylabel(r'$\bar{\nu}_{PPC}$ [$cm^{-1}$]')
MRE = np.mean(abs(freqs_PPC - freqs_PPC_ref)/freqs_PPC_ref)
plt.title(f'PPC Method Predicted frequencies\nMRE = {MRE:.1%}, acc. = {success_rate_PPC:.1%}')

xlim = plt.xlim()
ylim = plt.ylim()

rmin, rmax = min(xlim[0], ylim[0]), max(xlim[1], ylim[1])
plt.plot([rmin, rmax], [rmin, rmax], linestyle='dashed', c='k')
plt.plot([0, 0], [rmin, rmax], c='k', alpha=.3)
plt.plot([rmin, rmax], [0, 0], c='k', alpha=.3)
plt.gca().add_patch(Rectangle((rmin, 0), abs(rmin), rmax, facecolor='red', alpha=.2))
plt.gca().add_patch(Rectangle((0, rmin), rmax, abs(rmin), facecolor='red', alpha=.2))

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
