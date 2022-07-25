"""
Manually input delta chi tensor
"""

from audioop import findmax
import matplotlib.pyplot as plt
from sympy import Q
from paramagpy import protein, fit, dataparse, metal
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def calc_racs(pdb, ax=-6.342E-32, rh=-1.412E-32, beta = 1.36718):
	"""
	Calculates racs 
	"""
	# Load the PDB file
	prot = protein.load_pdb(pdb)

	#fcsa = prot[0]['A'][43]['FZ3'].csa
	#print("csa", fcsa)

	f = prot[0]['A'][43]['FZ3']
	flib = f.csa_lib
	#print(flib['FZ3'])

	flib['FZ3'] = ([-82.7*1E-6, -140.2*1E-6, -151.8*1E-6]), beta
	#print(f.csa_lib)
	# change beta of csa lib
	fcsa = f.csa
    #print("csa", fcsa)

	# Define an initial tensor
	m = metal.Metal(position = (1.288E-10, 10.021E-10, 30.398E-10), 
					eulers = (60.420*(np.pi/180),75.172*(np.pi/180), 107.840*(np.pi/180)), 
					axrh = (ax, rh))
	racs = m.racs(fcsa)
	#print("racs", racs)
	return racs

"""
The following makes a list for ax, rh, and degree values
"""

ax = []
for i in range(-100, 100, 10):
	a = i*1E-32
	ax.append(a)

rh = []
for i in range (-100, 100, 10):
	r = i*1E-32
	rh.append(r)
	

deg = []
for i in range(0,360,30):
	n = str(i) + "_deg_rot_gb1-ntaco-leap_adj.pdb"
	deg.append(n)

bet = []
racsBeta = []
for i in np.linspace(0.1,np.pi*2,25):
	bet.append(i)

def MAXracs(deg, ax, rh, bet):
	"""
	The following gets the max racs value within a certain 
	degree range, ax range, and rh range. 
	"""
	results = []
	for i in range(0, 12, 1):
		for j in range(7, 13, 1):
			for k in range(10, 11, 1):
				for l in range(0, 6, 1):
					result = ("racs", calc_racs(deg[i], ax[j], rh[k], bet[l]),
							  "deg", deg[i],"ax", ax[j],"rh", rh[j], "beta", bet[l])
					results.append(result)

	maxRacs = max(results)
	print("max", max)
	index = results.index(maxRacs)
	return maxRacs


degVals = {i.replace('_deg_rot_gb1-ntaco-leap_adj.pdb', '') for i in deg}
b = [int(i) for i in degVals]
b.sort()

racsDeg = []
for i in range (0,12,1):
	c = [calc_racs(deg[i])]
	racsDeg.append(c)


racsAx = []
for i in range(0,20,1):
		d = [calc_racs('0_deg_rot_gb1-ntaco-leap_adj.pdb', ax[i])]
		racsAx.append(d)
"""
plt.plot(ax, racsAx)
plt.xlabel('ax values')
plt.ylabel('Racs Values')
plt.title('Racs vs. ax')
#plt.show()
"""

racsRh = []
for i in range(0,20,1):
		e = [calc_racs('0_deg_rot_gb1-ntaco-leap_adj.pdb', rh[i])]
		racsRh.append(e)
"""
plt.plot(rh, racsRh)
plt.xlabel('rh values')
plt.ylabel('Racs Values')
plt.title('Racs vs. rh')
#plt.show()
"""
racsBeta = []
for i in np.linspace(0.1,np.pi*2,25):
		x = [calc_racs('0_deg_rot_gb1-ntaco-leap_adj.pdb', beta = i)]
		racsBeta.append(x)
print(racsBeta)
"""
plt.plot(bet, racsBeta)
plt.xlabel('beta values')
plt.ylabel('Racs Values')
plt.title('Racs vs. Beta')
plt.show()
"""

def maxconstrh(deg, bet):
	results2 = []
	for i in range(0, 12, 1):
		for l in range(0, 7, 1):
			result2 = ("racs", calc_racs(deg[i], bet[l]),
					  "deg", deg[i], "beta", bet[l])
			results2.append(result2)

		maxRacs2 = max(results2)
		print("max", max)
		index = results2.index(maxRacs2)
		return maxRacs2

#print(MAXracs(deg,ax, rh,bet))

#calc_racs('0_deg_rot_gb1-ntaco-leap_adj.pdb')


fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(b, racsDeg)
axs[0, 0].set_title('Racs vs. Degree')
axs[0, 0].set_xlabel('Degree Values')
axs[0, 0].set_ylabel('Racs Values')
axs[0, 1].plot(ax, racsAx, 'tab:orange')
axs[0, 1].set_title('Racs vs. ax')
axs[0, 1].set_xlabel('ax Values')
axs[0, 1].set_ylabel('Racs Values')
axs[1, 0].plot(rh, racsRh, 'tab:green')
axs[1, 0].set_title('Racs vs. rh')
axs[1, 0].set_xlabel('rh Values')
axs[1, 0].set_ylabel('Racs Values')
axs[1, 1].plot(bet, racsBeta, 'tab:red')
axs[1, 1].set_title('Racs vs. Beta')
axs[1, 1].set_xlabel('Beta Values')
axs[1, 1].set_ylabel('Racs Values')
plt.show()

#Q = racsAx[0][0] - racsAx[1][0]
#print(Q)
#P = racsAx[4][0] - racsAx[5][0]
#print(P)

#plt.scatter(ax, racsAx)
#plt.show()

# make 2D empty array

z = np.zeros(shape=(len(ax), len(rh)))

for numAx, i in enumerate(ax):
    for numRh, j in enumerate(rh):
        ans = calc_racs('0_deg_rot_gb1-ntaco-leap_adj.pdb', ax[numAx], rh[numRh])
        z[numAx,numRh] = ans

#print(z)
plt.pcolormesh(ax, rh, z)
plt.colorbar()
plt.xlabel("ax values")
plt.ylabel("rh values")
plt.title("RACS Values for ax and rh combination")
plt.show()