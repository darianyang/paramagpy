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

	fcsa = prot[0]['A'][43]['F03'].csa
	print("csa", fcsa)

	f = prot[0]['A'][43]['F03']
	flib = f.csa_lib
	print(flib['F03'])

	flib['F03'] = ([-17.294  -7.082 -21.674]), beta
	print(f.csa_lib)
	# change beta of csa lib
	fcsa = f.csa
    #print("csa", fcsa)

	# Define an initial tensor
	m = metal.Metal(position = (1.288E-10, 10.021E-10, 30.398E-10), 
					eulers = (60.420*(np.pi/180),75.172*(np.pi/180), 107.840*(np.pi/180)), 
					axrh = (ax, rh))
	racs = m.racs(fcsa)
	print("racs", racs)
	return racs