"""Physical constants in seconds, centimeters and MeV."""

# Set the base units
MEV = 1.0
KEV = 1e-3
CM = 1.0
ANGSTROM = 1e-8
SEC = 1.0

# Dimensionless constants
ALPHA = 1/137
K_PI = 3.14159

# Fundamental constants
C_LIGHT = 2.99792458e10 * (CM / SEC)  # Speed of light
HBAR = 6.58212e-22 * MEV * SEC  # Reduced Planck's constant (hbar)
HBARC = 1.97236e-11 * MEV * CM  # hbar * c
HC = 1.239e-10 * MEV * CM  # h c
KB = 8.617333262145e-11  # Boltzmann's constant in MeV K^-1
MU0 = 1.25663706212e-6  # Permeability of vacuum in kg / (A * s^2)
CHARGE_COULOMBS = 1.602176634e-19  # elementary charge in coulombs

# Scientific constants
AVOGADRO = 6.022e23

# Particle masses
M_E = 0.511 * MEV  # electron
M_P = 938.0 * MEV  # proton
M_MU = 105.658369 * MEV  # muon
M_TAU = 1.77699e3 * MEV  # tau
M_U = 2.2 * MEV  # up quark
M_D = 4.7 * MEV  # down quark
M_C = 1.280e3 * MEV  # charm quark
M_PI = 139.57 * MEV  # pi+ / pi-
M_PI0 = 134.98 * MEV  # pi0
M_ETA = 547.862 * MEV  # eta
M_N = 939.56 * MEV  # neutron
M_A = 1014 * MEV  # axial mass
M_K = 493.67 * MEV  # charged kaon

# Conversion between units
METER_BY_MEV = 6.58212e-22 * 2.998e8  # MeV*m
MEV_PER_KG = 5.6095887e29  # MeV/kg
MEV_PER_HZ = 6.58e-22  # 1 s^-1 = 6.58e-22 MeV
CM_PER_ANG = 1e-8  # cm to angstroms
KEV_CM = HBARC / KEV / ANGSTROM
MEV2_CM2 = (METER_BY_MEV * 100)**2
S_PER_DAY = 3600*24
HBARC_KEV_ANG = 1.97

# SM parameters
SSW = 0.2312  # Sine-squared of the Weinberg angle
G_F = 1.16638e-11  # MeV^-2
CABIBBO = 0.9743  # Cosine of the Cabibbo angle
R_E = 2.81794e-13  # classical electron radius in cm
V_UD = 0.9737  # ckm ud
V_US = 0.2245  # ckm us
F_PI = 130.2 * MEV  # pion decay constant
F_K = 155.7 * MEV  # kaon decay constant

KAON_TO_PION_XS = 0.75

BR_PI_MUNU = 0.999877
BR_K_MUNU = 0.6344

# Lifetimes
KAON_LIFETIME = 1.238e-8 * SEC
KSHORT_LIFETIME = 8.954e-11 * SEC
KLONG_LIFETIME = 5.116e-8 * SEC
PION_LIFETIME = 2.6e-8 * SEC
PI0_LIFETIME = 8.5e-17 * SEC

# Total widths in MeV
PION_WIDTH = 2.524e-14 * MEV
PI0_WIDTH = 1.0
KAON_WIDTH = 5.317e-14 * MEV

