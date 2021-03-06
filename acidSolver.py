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
constant = -acidK * waterK

# Setup the cubic and solve for each base volume
polyCoeff = [1.0, quadraticTerm, linearTerm, constant]
eqH3OConc = np.max(np.roots(polyCoeff))

# Calculate Important Quantities
eqOHConc = waterK / eqH3OConc
eqAConc = eqAcidAnionConc + eqH3OConc - eqOHConc
eqHAConc = eqH3OConc * eqAConc / acidK

pH = -np.log10(eqH3OConc)
pOH = -np.log10(eqOHConc)
pA = -np.log10(eqAConc)
pHA = -np.log10(eqHAConc)

# Return Result
print(
	"At equilibrium: [H3O+] = {0:.2E}, [OH-] = {1:.2E}, [A-] = {2:.2E}, [HA] = {3:.2Ef}".
	format(eqH3OConc, eqOHConc, eqAConc, eqHAConc))
print(
	"At equilibrium: pH = {0:.2f}, pOH = {1:.2f}, pA = {2:.2f}, pHA = {3:.2f}".
	format(pH, pOH, pA, pHA))

