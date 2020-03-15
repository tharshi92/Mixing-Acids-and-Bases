import numpy as np

# Setup problem
acidK = 1e-5
waterK = 1e-14

initBufferConc = 1e-2
initBufferVol = 0e+2
initAcidConc = 1e-2
initAcidVol = 1e+2

# Calculate intermediate Variables
totalVol = initAcidVol + initBufferVol
eqAcidAnionConc = initBufferConc * initBufferVol / totalVol

# Create cubic coeffcients
quadraticTerm = acidK + eqAcidAnionConc
linearTerm = -(waterK + acidK * (initAcidConc * initAcidVol / totalVol))
constant = - acidK * waterK

# Setup the cubic and solve for each base volume
polyCoeff = [1.0, quadraticTerm, linearTerm, constant]
eqH3OConc = np.max(np.roots(polyCoeff))

# Calculate Important Quantities
eqOHConc = waterK / eqH3OConc
eqAConcs = eqAcidAnionConc + eqH3OConc - eqOHConc
eqHAConc = eqH3OConc * eqAConcs / acidK

pH = - np.log10(eqH3OConc)
pOH = - np.log10(eqOHConc)
pHA = - np.log10(eqHAConc)
pA = - np.log10(eqAConcs)

# Return Result
print("At equilibrium: pH = {0:.2f}, pOH = {1:.2f}, pHA = {2:.2f}, pA = {3:.2f}".format(pH, pOH, pHA, pA))
