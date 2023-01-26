import scm.plams as plams


def default(preset='BLYP-D3(BJ)/TZ2P/Good'):
    '''
    Settings are given in the form:
        Functional[-Dispersion]/Basisset/Numquality
    -Dispersion is optional, omitting it will simply turn off dispersion corrections
    '''
    sett = plams.Settings()
    sett.input.ams.task = 'SinglePoint'

    if preset == 'Cheap':
        sett.input.adf.Relativity.Level = 'Scalar'
        sett.input.adf.Basis.Core = 'None'
        return sett
    if preset == 'DFTB':
        sett.input.ams.task = 'SinglePoint'
        sett.input.DFTB
        sett.input.DFTB.model = 'GFN1-xTB'
        return sett

    decode = preset.split('/')

    functional, basis_size, numerical_quality = decode
    sett.input.adf.basis.Type = basis_size
    sett.input.adf.basis.Core = 'None'

    sett.input.adf.Relativity.Level = 'Scalar'

    # dispersion corrections
    if functional.endswith('-D4(EEQ)'):
        sett.input.adf.XC.Dispersion = 'GRIMME4'
        functional = functional[:-8]
    elif functional.endswith('-D3(BJ)'):
        sett.input.adf.XC.Dispersion = 'GRIMME3 BJDAMP'
        functional = functional[:-7]
    elif functional.endswith('-D3'):
        sett.input.adf.XC.Dispersion = 'GRIMME3'
        functional = functional[:-3]
    elif functional.endswith('-dDsC'):
        sett.input.adf.XC.Dispersion = 'dDsC'
        functional = functional[:-5]
    elif functional.endswith('-dUFF'):
        sett.input.adf.XC.Dispersion = 'UFF'
        functional = functional[:-5]
    elif functional.endswith('-D'):
        sett.input.adf.XC.Dispersion = 'Default'
        functional = functional[:-2]

    # GGA functionals
    if any(fun in functional for fun in [
           'BLYP', 'BP86', 'OLYP', 'OPBE', 'PBE', 'PBEsol']):
        sett.input.adf.XC.GGA = functional
    # Hybrid functionals
    elif any(fun in functional for fun in ['B3LYP', 'B1LYP']):
        sett.input.adf.XC.Hybrid = functional
    # MetaGGA
    elif any(fun in functional for fun in ['M06L', 'MN15-L', 'SCAN']):
        sett.input.adf.XC.MetaGGA = functional
    # LDA is standard functional
    elif functional == 'LDA':
        pass

    # numerical quality
    sett.input.adf.NumericalQuality = numerical_quality
    return sett


def single_point(sett=None, preset='BLYP-D3(BJ)/TZ2P/Good'):
    if sett is None:
        sett = default(preset)

    sett.input.ams.task = 'SinglePoint'
    return sett


def optimization(sett=None, preset='BLYP-D3(BJ)/TZ2P/Good'):
    if sett is None:
        sett = default(preset)

    sett.input.ams.task = 'GeometryOptimization'
    return sett


def transition_state(a, b, f=1.0, sett=None, preset='BLYP-D3(BJ)/TZ2P/Good'):
    if sett is None:
        sett = default(preset)

    sett.input.ams.task = 'TransitionStateSearch'
    try:
        sett.input.ams.TransitionStateSearch.ReactionCoordinate.Distance.append(
            f'{a} {b} {f}')
    except BaseException:
        sett.input.ams.TransitionStateSearch.ReactionCoordinate.Distance = [
            f'{a} {b} {f}']
    return sett


def pes_scan(constraint, n=15, sett=None, preset='BLYP-D3(BJ)/TZ2P/Good'):
    if sett is None:
        sett = default(preset)

    sett.input.ams.task = 'PESScan'
    sett.input.ams.PESScan.ScanCoordinate.nPoints = n
    if type(constraint[0]) in [list, tuple]:
        sett.input.ams.PESScan.ScanCoordinate.Distance = [
            ' '.join([str(x) for x in const]) for const in constraint]
    else:
        sett.input.ams.PESScan.ScanCoordinate.Distance = [
            ' '.join([str(x) for x in constraint])]
    return sett


def vibrations(sett=None, preset='BLYP-D3(BJ)/TZ2P/Good'):
    if sett is None:
        sett = default(preset)

    sett.input.ams.Properties.NormalModes = 'Yes'
    sett.input.ams.Properties.PESPointCharacter = 'Yes'
    sett.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = '-50'
    sett.input.ams.NormalModes.ReScanFreqRange = '-10000000.0 10.0'
    return sett


def distance_constraint(*args, sett=None, preset='BLYP-D3(BJ)/TZ2P/Good'):
    if sett is None:
        sett = default(preset)
    args = [str(arg) for arg in args]
    if isinstance(sett.input.ams.Constraints.Distance, list):
        sett.input.ams.Constraints.Distance.append(' '.join(args))
    else:
        sett.input.ams.Constraints.Distance = [' '.join(args)]
    return sett


def frozen_atoms(idxs, sett=None, preset='BLYP-D3(BJ)/TZ2P/Good'):
    for idx in idxs:
        if type(sett.input.ams.Constraints.Atom) is list:
            sett.input.ams.Constraints.Atom.append(str(idx))
        else:
            sett.input.ams.Constraints.Atom = [str(idx)]
    return sett



def charge(ch=0, sett=None, preset='BLYP-D3(BJ)/TZ2P/Good'):
    if sett is None:
        sett = default(preset)

    sett.input.ams.System.Charge = ch
    return sett


def spin_polarization(sp=0., sett=None, preset='BLYP-D3(BJ)/TZ2P/Good'):
    if sett is None:
        sett = default(preset)

    sett.input.adf.SpinPolarization = sp
    if sp != 0.:
        sett.input.adf.Unrestricted = 'Yes'
    return sett


def constant_bond_dists(mol, sett=None, preset='BLYP-D3(BJ)/TZ2P/Good'):
    if sett is None:
        sett = default(preset)

    for bond in mol.bonds:
        a1i, a2i = mol.atoms.index(bond.atom1) + \
            1, mol.atoms.index(bond.atom2) + 1
        dist = bond.length()
        distance_constraint(a1i, a2i, dist, sett=sett)

    return sett


if __name__ == '__main__':
    s = default('LDA/SZ/Basic')
    print(s)
    # optimization(s)
    # print(s)
    # vibrations(s)
    # print(s)
    # charge(s, 10)
    # print(s)
    # spin_polarization(s, 5000)
    # print(s)
