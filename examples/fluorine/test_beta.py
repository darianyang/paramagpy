from paramagpy import protein, fit, dataparse, metal
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def calc_racs(pdb, ax=-6.342E-32, rh=-1.412E-32, beta=1.36718):
	"""
	Calculates racs 
	"""
	# Load the PDB file
	prot = protein.load_pdb(pdb)

	f = prot[0]['A'][43]['FE3']
	flib = f.csa_lib
	print(flib['FE3'])

	flib['FE3'] = ([-66.6*1E-6, -126.1*1E-6 , -190.6*1E-6]), beta
	#print(f.csa_lib)
	# change beta of csa lib
	fcsa = f.csa
    #print("csa", fcsa)

	m = metal.Metal(position = (1.288E-10, 10.021E-10, 30.398E-10), 
					eulers = (60.420*(np.pi/180),75.172*(np.pi/180), 107.840*(np.pi/180)), 
					axrh = (ax, rh))
	racs = m.racs(fcsa)

	return racs

r = calc_racs('allF_gb1-ntaco-leap_adj.pdb', beta=1)
print(r)