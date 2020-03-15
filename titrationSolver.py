import numpy as np
import matplotlib.pyplot as plt

# Graphics parameters

numPoints = 301

# Setup problem
acidConstant = 1e-5
waterConstant = 1e-14

initBufferConc = 1e-2
initBufferVol = 0e+2
initAcidConc = 1e-2
initAcidVol = 1e+2
baseConc = 1e-2

# Create variable base volume array
minVol = 0.0
maxVol = 2.0 * initAcidConc * initAcidVol / baseConc
baseVols = np.linspace(minVol, maxVol, numPoints)

# Calculate intermediate Variables
totalVols = baseVols + initAcidVol + initBufferVol
eqBaseCationConcs = baseConc * baseVols / totalVols
eqAcidAnionConcs = initBufferConc * initBufferVol / totalVols

# Create cubic coeffcients
quadraticTerms = acidConstant + eqBaseCationConcs + eqAcidAnionConcs
linearTerms = -(waterConstant + acidConstant * (
	(initAcidConc * initAcidVol / totalVols) - eqBaseCationConcs))
constant = -acidConstant * waterConstant

# Solve for every base volume
eqH3OConcs = np.linspace(minVol, maxVol, numPoints)

for i in range(numPoints):
	# Setup the cubic and solve for each base volume
	polyCoeffs = [1.0, quadraticTerms[i], linearTerms[i], constant]
	eqH3OConcs[i] = np.max(np.roots(polyCoeffs))

# Calculate Important Quantities
eqOHConcs = waterConstant / eqH3OConcs

eqAConcs = eqBaseCationConcs + eqAcidAnionConcs + eqH3OConcs - eqOHConcs

eqHAConcs = eqH3OConcs * eqAConcs / acidConstant

pHs = -np.log10(eqH3OConcs)
pOHs = -np.log10(eqOHConcs)
pHAs = -np.log10(eqHAConcs)
pAs = -np.log10(eqAConcs)

#print(np.argwhere(pHs > 9))
print("The pH at the equivalence point is " + str(pHs[(numPoints - 1) / 2]) +
						".")

# Plot Results
plt.figure(1)
plt.title('Power Transformations vs Base Added (mL)')
plt.plot(baseVols, pHs, 'r', label="pH")
plt.plot(baseVols, pOHs, 'c', label="pOH")
plt.plot(baseVols, pHAs, 'g', label="pHA")
plt.plot(baseVols, pAs, 'm', label="pA")
plt.grid('on')
plt.legend(loc=0)
plt.savefig('powerTransforms.png')
plt.show()

