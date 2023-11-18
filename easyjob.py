from scm import plams
from TCutility import log


class Job:
    def __init__(self, **kwargs):
        self.settings = plams.settings()
        self.molecule = plams.Molecule()
        self.name = 'calc'
        self.path = 'tmp'
        self.functional = None
        self.basis = None
        self.numerical_quality = None
        self.task = 'SinglePoint'
        self.charge = 0
        self.spin_polarization = 0
        self.unrestricted = False
        self.solvent = None
        self.vibrations = False

        for key, value in kwargs.items():
            setattr(self, key, value)

    def build_settings(self):


    # def set_functional(self, functional: str, use_libxc: bool = False):
    #     '''
    #     Set the functional for this calculation.

    #     Args:
    #         functional: name of the functional to use.
    #         use_libxc: whether the functional is a member of LibXC. This will set self.settings.input.adf.XC.LibXC = functional.
    #     '''
    #     self.functional = functional
    #     # LDA functionals
    #     if functional in ['VWN', 'PW92']:
    #         self.settings.input.adf.XC.LDA = functional

    #     # GGA functionals
    #     elif functional in ['BLYP', 'BP86', 'GAM', 'HTBS', 'KT1', 'KT2', 'mPW', 'mPBE', 'N12', 'OLYP', 'OPBE', 'PBE', 'PBEsol', 'PW91', 'revPBE', 'RPBE', 'BEE']:
    #         self.settings.input.adf.XC.GGA = functional

    #     # Hybrid functionals
    #     # BH&H and BHandH should be equivalent
    #     elif functional in ['B3LYP', 'B1LYP', 'B1PW91', 'B3LYP*', 'BHandH', 'BHandHLYP', 'KMLYP', 'MPW1PW', 'MPW1K', 'O3LYP', 'OPBE0', 'PBE0', 'S12h', 'X3LYP', 'HTBS']:
    #         self.settings.input.adf.XC.Hybrid = functional

    #     # MetaGGA
    #     elif functional in ['M06L', 'MN15-L', 'MVS', 'SCAN', 'revTPSS', 'SSB', 'TASKxc', 'TPSS', 'r2SCAN-3c']:
    #         self.settings.input.adf.XC.MetaGGA = functional
    #     elif functional in ['rSCAN', 'revSCAN', 'r2SCAN']:
    #         self.settings.input.adf.XC.LibXC = functional

    #     # range separated
    #     elif functional in ['LCY-BLYP', 'LCY-BP86', 'LCY-PBE', 'CAM-B3LYP', 'CAMY-B3LYP', 'HSE03', 'HSE06', 'M11', 'MN12-SX', 'N12-SX', 'WB97', 'WB97X']:
    #         self.settings.input.adf.XC.LibXC = functional

    #     # DoubleHybrid
    #     elif functional in ['rev-DOD-PBEP86', 'rev-DOD-BLYP', 'rev-DOD-PBE', 'B2PLYP', 'B2GPPLYP']:
    #         self.settings.input.adf.XC.DoubleHybrid = functional

    #     # MetaHybrid
    #     elif functional in ['MN15', 'M06', 'M06-2X', 'M06-HF', 'TPSSH']:
    #         self.settings.input.adf.XC.MetaHybrid = functional
    #     elif functional in ['revSCAN0']:
    #         self.settings.input.adf.XC.LibXC = functional

    #     # libxc
    #     elif functional == 'BMK':
    #         self.settings.input.adf.XC.LibXC = 'HYB_MGGA_X_BMK GGA_C_BMK'

    #     # LDA is standard functional
    #     elif functional == 'LDA':
    #         pass

    #     # SAOP is a special model functional
    #     elif functional == 'SAOP':
    #         self.settings.input.adf.XC.model = 'SAOP'

    #     else:
    #         if not use_libxc:
    #             raise ValueError(f'XC-functional {functional} not defined, if you meant to use a libXC functional, please call this method with "use_libxc=True"')
    #         else:
    #             self.settings.input.adf.XC.LibXC = functional

    #     if functional == 'r2SCAN-3c' and self.basis != 'mTZ2P' and self.numerical_quality != 'good':
    #         log.warn(f'Use of r2SCAN-3c/{self.basis}/{self.numerical_quality} is not recommended.\nUse r2SCAN-3c/mTZ2P/good instead')

    # def set_basis(self, basis: str):
    #     self.basis = basis
    #     self.settings.input.adf.basis.Type = basis

    # def set_numerical_quality(self, numerical_quality: str):
    #     self.numerical_quality = numerical_quality
    #     self.settings.input.adf.NumericalQuality = numerical_quality

    # def set_charge(self, charge: int):
    #     self.settings.input.ams.System.Charge = charge

    # def set_spinpol(self, spinpol: int):
    #     self.settings.input.ams.SpinPolarization = spinpol
    #     if spinpol != 0.:
    #         self.settings.input.input.adf.Unrestricted = 'Yes'
    


