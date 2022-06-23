from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
#prot = protein.load_pdb('m01_1gb1FH.pdb')
prot = protein.load_pdb('150_deg_rot_gb1-ntaco-leap_adj.pdb')

# Load the PCS data
rawData = dataparse.read_pcs('NTA-Cobalt(II)_19F_PCS.npc')

# Associate PCS data with atoms of the PDB
parsedData = prot.parse(rawData)

# Define an initial tensor
mStart = metal.Metal()

# Set the starting position to an atom close to the metal
mStart.position = prot[0]['A'][32]['NE2'].position

# Calculate an initial tensor from an SVD gridsearch
[mGuess], [data] = fit.svd_gridsearch_fit_metal_from_pcs(
	[mStart],[parsedData], radius=1, points=10)

# Refine the tensor using non-linear regression
[mFit], [data] = fit.nlr_fit_metal_from_pcs([mGuess], [parsedData], useracs=True)

# Calculate the Q-factor
qfac = fit.qfactor(data)

# Save the fitted tensor to file
mFit.save('test_PCS_tensor.txt')
print("mfit", mFit)

#### Plot the correlation ####
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

# Plot the data
ax.plot(data['exp'], data['cal'], marker='o', lw=0, ms=3, c='r',
	label="Q-factor = {:5.4f}".format(qfac))

# Plot a diagonal
l, h = ax.get_xlim()
ax.plot([l,h],[l,h],'-k',zorder=0)
ax.set_xlim(l,h)
ax.set_ylim(l,h)

# Make axis labels and save figure
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
#fig.savefig("pcs_fit.png")
#plt.show()
