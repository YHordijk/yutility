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


# now we distill some information from the raw data
methods = set(data.timing.keys())
molnames = list(set([key for method in methods for key in data.timing[method].keys()]))
natoms = [data.natoms[molname] for molname in molnames]

# first get the timing_ratios
for method in methods:
    data.timing_ratios[method] = []
    for molname in molnames:
        try:
            full_method_timing = data.timing['dft_freq'][molname]
            timing = data.timing[method][molname]
            data.timing_ratios[method].append(full_method_timing/timing)
        except:
            data.timing_ratios[method].append(None)

# get the accuracies
for method in methods:
    nsuccess = 0
    nfailed = 0
    for molname in molnames:
        try:
            full_method_character = data.character['dft_freq'][molname]
            character = data.character[method][molname]
            if character == full_method_character:
                nsuccess += 1
            else:
                nfailed += 1
        except:
            raise
    data.accuracy[method] = nsuccess / (nsuccess + nfailed)

print(data.timing_ratios)
print(data.accuracy)

    # if not os.path.exists(j('calculations', xyzfile.removesuffix('.xyz'))):
    #     continue
    # if not os.path.exists(j('calculations', xyzfile.removesuffix('.xyz') + '_PPC')):
    #     continue


    # fullfreq_res = results.read(j('calculations', xyzfile.removesuffix('.xyz'), 'dft_freq'))
    # DFTB_res = results.read(j('calculations', xyzfile.removesuffix('.xyz'), 'low_freq'))
    # PPC_res = results.read(j('calculations', xyzfile.removesuffix('.xyz') + '_PPC', 'calc'))
    # refine_res = results.read(j('calculations', xyzfile.removesuffix('.xyz'), 'refine'))
    # # refine_res = results.read(j('calculations', xyzfile.removesuffix('.xyz') + '_2', 'refine'))


    # if not refine_res.properties.vibrations.frequencies:
    #     nfailed += 1
    # #     continue

    # if fullfreq_res.properties.vibrations.character == refine_res.properties.vibrations.character:
    #     nsuccess += 1
    # else:
    #     nfailed += 1


    # if fullfreq_res.properties.vibrations.character == DFTB_res.properties.vibrations.character:
    #     nsuccess_DFTB += 1
    # else:
    #     nfailed_DFTB += 1


    # if fullfreq_res.properties.vibrations.character == PPC_res.properties.vibrations.character:
    #     nsuccess_PPC += 1
    # else:
    #     nfailed_PPC += 1

    # character.append(fullfreq_res.properties.vibrations.character)
    # character_pred.append(refine_res.properties.vibrations.character)

    # # analyse the performance boost and errors
    # timing_ratio = fullfreq_res.timing.total/refine_res.timing.total
    # timing_full.append(fullfreq_res.timing.total/60)
    # timing_refine.append(refine_res.timing.total/60)
    # timing_PPC.append(PPC_res.timing.total/60)
    # timing_DFTB.append(DFTB_res.timing.total/60)
    # timing_ratios.append(timing_ratio)
    # timing_ratios_DFTB.append(fullfreq_res.timing.total/DFTB_res.timing.total)
    # timing_ratios_PPC.append(fullfreq_res.timing.total/PPC_res.timing.total)


    # freqs_ref.extend(fullfreq_res.properties.vibrations.frequencies[:nfreqs])
    # if PPC_res.properties.vibrations.frequencies:
    #     nfreqs_PPC = min(len(PPC_res.properties.vibrations.frequencies), nfreqs)
    #     freqs_PPC.extend(PPC_res.properties.vibrations.frequencies[:nfreqs_PPC])
    #     freqs_PPC_ref.extend(fullfreq_res.properties.vibrations.frequencies[:nfreqs_PPC])

    # freqs_DFTB.extend(DFTB_res.properties.vibrations.frequencies[:3])
    # freqs_DFTB_ref.extend(fullfreq_res.properties.vibrations.frequencies[:3])
    # freqs.extend(refine_res.properties.vibrations.frequencies[:nfreqs])
    # modes_ref.extend(fullfreq_res.properties.vibrations.modes[:nfreqs])
    # modes.extend(refine_res.properties.vibrations.modes[:nfreqs])
    # natoms.append(len(mol['molecule']))
    # errors.append(abs((fullfreq_res.properties.vibrations.frequencies[0] - refine_res.properties.vibrations.frequencies[0])/fullfreq_res.properties.vibrations.frequencies[0] * 100))

print(max(timing_ratios))
print(max(timing_ratios_DFTB))
# print(max(timing_ratios_PPC))
print([tim for tim, nat in zip(timing_ratios_PPC, natoms) if nat == 47])
timing_ratios = np.array(timing_ratios)
natoms = np.array(natoms)
errors = np.array(errors)
success_rate = nsuccess/(nsuccess + nfailed)
success_rate_DFTB = nsuccess_DFTB/(nsuccess_DFTB + nfailed_DFTB)
success_rate_PPC = nsuccess_PPC/(nsuccess_PPC + nfailed_PPC)
print(f'Success rate: {nsuccess/(nsuccess + nfailed):.0%}')
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
plt.title(r'Time Gain Ratio ($\gamma_{DFTB}$)')
popt = curve_fit(curve, natoms, np.array(timing_full)/np.array(timing_DFTB), p0=[2, 0, 7e-4])[0]
plt.plot(np.linspace(min(natoms), max(natoms), 100), curve(np.linspace(min(natoms), max(natoms), 100), *popt), 
         linestyle='dashed', label=rf'$\gamma_{{DFTB}} = {popt[2]:.1f}(N_{{atom}}{popt[1]:+.1f})^{{{popt[0]:.1f}}}$')
plt.scatter(natoms, np.array(timing_full)/np.array(timing_DFTB))
plt.xlabel('Number of atoms')
plt.ylabel(r'$\gamma$')
plt.legend()


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
