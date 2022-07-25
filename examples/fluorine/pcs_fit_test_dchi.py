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

	fcsa = prot[0]['A'][43]['FZ3'].csa

	f = prot[0]['A'][43]['FZ3']
	flib = f.csa_lib

	#flib['FE3'] = ([-66.6*1E-6, -126.1*1E-6, -190.6*1E-6]), beta
	flib['FZ3'] = ([-82.7*1E-6, -140.2*1E-6, -151.8*1E-6]), beta
	#flib['FZ2'] = ([-73.2*1E-6, -126.1*1E-6, -201.1*1E-6]), beta
	#flib['FH2'] = ([-64.9*1E-6, -129*1E-6, -169.4*1E-6]), beta
	fcsa = f.csa

	# Define an initial tensor
	m = metal.Metal(position = (1.288E-10, 10.021E-10, 30.398E-10), 
					eulers = (60.420*(np.pi/180),75.172*(np.pi/180), 107.840*(np.pi/180)), 
					axrh = (ax, rh))
	racs = m.racs(fcsa)
	return racs

if __name__ == "__main__":

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

	"""
	deg = []
	for i in range(0,360,18):
		deg.append(i)	
	"""
	deg = []
	for i in range(0,360,18):
		n = str(i) + ".pdb"
		deg.append(n)

	degVals = []
	for i in range(0,360,18):
		n = i
		degVals.append(n)

	beta = []
	for i in np.linspace(0.1,np.pi*2,20):
		beta.append(i)

	def MAXracs(deg, ax, rh, beta):
		"""
		The following gets the max racs value within a certain 
		degree range, ax range, and rh range. 
		"""
		results = []
		for i in range(0, 12, 1):
			for j in range(7, 13, 1):
				for k in range(10, 11, 1):
					for l in range(0, 6, 1):
						result = ("racs", calc_racs(deg[i], ax[j], rh[k], beta[l]),
								"deg", deg[i],"ax", ax[j],"rh", rh[j], "beta", beta[l])
						results.append(result)

		maxRacs = max(results)
		print("max", max)
		index = results.index(maxRacs)
		return maxRacs
	"""
	allDegRacs = np.loadtxt('racs.txt')
	lines_to_read = []
	for i in range (0,360,18):
		lines_to_read.append(i)
	
	degRacs = []
	for position, line in enumerate(allDegRacs):
		if position in lines_to_read:
			degRacs.append(line)
	"""
	racsDeg = []
	for i in range (0,20,1):
		c = [calc_racs(deg[i])]
		racsDeg.append(c)

	racsAx = []
	for i in range(0,20,1):
			d = [calc_racs('allF_gb1-ntaco-leap_adj.pdb', ax = ax[i])]
			racsAx.append(d)

	racsRh = []
	for i in range(0,20,1):
			e = [calc_racs('allF_gb1-ntaco-leap_adj.pdb', rh = rh[i])]
			racsRh.append(e)

	racsBeta = []
	for i in range(0,20,1):
			x = [calc_racs('allF_gb1-ntaco-leap_adj.pdb', beta = beta[i])]
			racsBeta.append(x)

	def maxconstrh(deg, beta):
		results2 = []
		for i in range(0, 12, 1):
			for l in range(0, 7, 1):
				result2 = ("racs", calc_racs(deg[i], beta[l]),
						"deg", deg[i], "beta", beta[l])
				results2.append(result2)

			maxRacs2 = max(results2)
			print("max", max)
			index = results2.index(maxRacs2)
			return maxRacs2

	"""
	The following code makes graphs of 
	racs vs. beta, racs vs. deg, racs vs, ax,
	and racs vs rh.
	"""
	fig, axs = plt.subplots(2, 2)
	axs[0, 0].plot(deg, racsDeg)
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
	axs[1, 1].plot(beta, racsBeta, 'tab:red')
	axs[1, 1].set_title('Racs vs. Beta')
	axs[1, 1].set_xlabel('Beta Values')
	axs[1, 1].set_ylabel('Racs Values')
	#plt.savefig("Racs_subplots", dpi=300, transparent = True)
	plt.show()

	"""
	The following makes heatplots of the racs 
	values dependent on ax and rh, ax and beta,
	rh and beta, and degree and ax
	"""
	def heat_plot_AxRh(ax, rh): 

		z = np.zeros(shape=(len(ax), len(rh)))

		for numAx, i in enumerate(ax):
			for numRh, j in enumerate(rh):
				ans = calc_racs('allF_gb1-ntaco-leap_adj.pdb', ax[numAx], rh[numRh])
				z[numAx,numRh] = ans

		#print(z)
		plt.pcolormesh(ax, rh, z, cmap="seismic", vmin=-0.1, vmax=0.1)
		plt.colorbar()
		plt.xlabel("ax values")
		plt.ylabel("rh values")
		plt.title("RACS Values for ax and rh combination")
		plt.savefig("Racs_axRh_comb", dpi=300, transparent = True)
		plt.show()

	def heat_plot_AxBeta(ax, beta): 

		z = np.zeros(shape=(len(ax), len(beta)))

		for numAx, i in enumerate(ax):
			for numBet, j in enumerate(beta):
				ans = calc_racs('allF_gb1-ntaco-leap_adj.pdb', ax = ax[numAx], beta = beta[numBet])
				z[numAx,numBet] = ans

		#print(z)
		plt.pcolormesh(ax, beta, z, cmap="seismic", vmin=-0.1, vmax=0.1)
		plt.colorbar()
		plt.xlabel("ax values")
		plt.ylabel("beta values")
		plt.title("RACS Values for ax and beta combination")
		plt.savefig("Racs_axBeta_comb", dpi=300, transparent = True)
		plt.show()

	def heat_plot_RhBeta(rh, beta): 

		z = np.zeros(shape=(len(rh), len(beta)))

		for numRh, i in enumerate(rh):
			for numBet, j in enumerate(beta):
				ans = calc_racs('allF_gb1-ntaco-leap_adj.pdb', rh = rh[numRh], beta = beta[numBet])
				z[numRh,numBet] = ans

		#print(z)
		plt.pcolormesh(rh, beta, z, cmap="seismic", vmin=-0.1, vmax=0.1)
		plt.colorbar()
		plt.xlabel("rh values")
		plt.ylabel("beta values")
		plt.title("RACS Values for rh and beta combination")
		plt.savefig("Racs_rhBeta_comb", dpi=300, transparent = True)
		plt.show()
	
	
	def heat_plot_degAx(deg, ax): 

		z = np.zeros(shape=(len(deg), len(ax)))

		for numB, i in enumerate(deg):
			for numAx, j in enumerate(ax):
				ans = calc_racs(deg[numB], ax[numAx])
				z[numB,numAx] = ans

		#print(z)
		plt.pcolormesh(degVals, ax, z, cmap="seismic", vmin=-0.1, vmax=0.1)
		plt.colorbar()
		plt.xlabel("degree values")
		plt.ylabel("ax values")
		plt.title("RACS Values for degree and ax combination")
		plt.show()

	def heat_plot_degbeta(deg, beta): 

		z = np.zeros(shape=(len(deg), len(beta)))

		for numB, i in enumerate(deg):
			for numbet, j in enumerate(beta):
				ans = calc_racs(deg[numB], beta = beta[numbet])
				z[numB,numbet] = ans

		#print(z)
		plt.pcolormesh(degVals, beta, z, vmin=-0.1, vmax= 0.1, cmap='seismic')
		plt.colorbar()
		plt.xlabel("degree values")
		plt.ylabel("beta values")
		plt.title("RACS Values for degree and beta combination")
		plt.show()

	def heat_plot_degRh(deg, rh): 

		z = np.zeros(shape=(len(deg), len(rh)))

		for numB, i in enumerate(deg):
			for numRh, j in enumerate(rh):
				ans = calc_racs(deg[numB], rh = rh[numRh])
				z[numB,numRh] = ans

		print(z[np.where(deg=='0.pdb'),np.where(rh==0)])
		print(np.where(deg=='0.pdb'),np.where(rh==0))
		plt.pcolormesh(degVals, rh, z, cmap="seismic", vmin=-0.1, vmax=0.1)
		plt.colorbar()
		plt.xlabel("degree values")
		plt.ylabel("rh values")
		plt.title("RACS Values for degree and rh combination")
		plt.show()

#print(calc_racs('342.pdb', beta = 6.283))

heat_plot_degbeta(deg,beta)
#heat_plot_AxRh(ax, rh)
#heat_plot_AxBeta(ax,beta)
#heat_plot_degAx(deg,ax)
#heat_plot_RhBeta(rh,beta)
#heat_plot_degRh(deg,rh)

print(calc_racs('0.pdb'))