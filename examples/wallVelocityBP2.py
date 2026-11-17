"""
Example where we compute the wall velocity using various approximations.
The model considered here corresponds to BP2 in 2609.xxxxx.
"""

from HydrodynamicsFOPT import Hydrodynamics, EntropyNLTE, EntropyBallistic, EntropyInterpolated

# WallGo result
vwWG = 0.544

# Hydrodynamics parameters corresponding to BP2
Tn = 98.39873456136364
vn = 2.13167139e+02
alN = 0.01325992878519393
psi = 0.9745633990796925
cb2 = 0.32467118701022657
cs2 = 0.33286025097450783
gstar = 107.75
wallWidth = 0.032

# Creating the Hydrodynamics object
hydro = Hydrodynamics(alN, cb2, cs2, psi, Tn)

# Properties of the out-of-equilibrium species (topL, topR and W)
# Masses
mTopL = 173*vn/246
mTopR = 173*vn/246
mW = 80.34*vn/246
# W thermal mass
mWT = 0.6256669389784285*Tn

# Degrees of freedom
dofTopL = 6
dofTopR = 6
dofW = 9

# Statistics (1 for fermion, -1 for boson)
statisticTopL = 1
statisticTopR = 1
statisticW = -1

# Relaxation time tau (in units of 1/Tn)
tauTopL = 19.58409902
tauTopR = 59.16350693
tauW = 127.88493832

# Creating the Entropy objects
entropyNLTE = EntropyNLTE(massesSymmetricPhase=[0, 0, mWT],
                          massesBrokenPhase=[mTopL, mTopR, (mW**2+mWT**2)**0.5],
                          dofs=[dofTopL, dofTopR, dofW],
                          statistics=[statisticTopL, statisticTopR, statisticW],
                          tau=[tauTopL, tauTopR, tauW],
                          L=wallWidth,
                          Tn=Tn,
                          gstar=gstar,
                          cs2=cs2)
entropyBall = EntropyBallistic(massesSymmetricPhase=[0, 0, mWT],
                               massesBrokenPhase=[mTopL, mTopR, (mW**2+mWT**2)**0.5],
                               dofs=[dofTopL, dofTopR, dofW],
                               statistics=[statisticTopL, statisticTopR, statisticW],
                               Tn=Tn,
                               gstar=gstar,
                               cs2=cs2)
entropyInterp = EntropyInterpolated(massesSymmetricPhase=[0, 0, mWT],
                                    massesBrokenPhase=[mTopL, mTopR, (mW**2+mWT**2)**0.5],
                                    dofs=[dofTopL, dofTopR, dofW],
                                    statistics=[statisticTopL, statisticTopR, statisticW],
                                    tau=[tauTopL, tauTopR, tauW],
                                    L=wallWidth,
                                    Tn=Tn,
                                    gstar=gstar,
                                    cs2=cs2,
                                    interpolateTotal=False)

# Compute all the wall velocities
print('Wall velocities:')
print('LTE:', hydro.findVwLTE())
print('NLTE:', hydro.findVwLTE(entropyNLTE))
print('Ballistic:', hydro.findVwLTE(entropyBall))
print('Interpolated:', hydro.findVwLTE(entropyInterp))
print('WallGo:', vwWG)