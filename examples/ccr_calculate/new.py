from paramagpy import protein, fit, dataparse, metal
import numpy as np

data = """2 , 3.5±0.5, 3.5±0.3, 4.0±0.6
3 , 4.8±0.4, 5.3±0.4, 6.5±0.5
4 , 5.5±0.4, 6.3±0.4, 5.9±0.6
5 , 4.5±0.7, 4.4±0.3, 4.1±0.7
6 , 4.1±0.5, 3.9±0.3, 4.6±0.5
7 , 4.8±0.7, 4.6±0.3, 4.0±0.8
8 ,        , 4.1±0.4, 2.2±0.6 
9 , 3.4±0.3, 3.2±0.3, 4.4±0.6
10, 3.0±0.4, 3.1±0.2, 4.0±0.3
11, 4.8±0.4, 4.0±0.4, 5.2±0.5
12, 3.8±0.5, 4.0±0.3, 6.0±0.8
13, 2.4±0.4, 2.6±0.3, 4.1±0.6
14, 4.6±0.4, 4.8±0.3, 7.0±0.5
15, 3.8±0.5, 4.0±0.4, 3.5±0.5
16, 2.9±0.4, 3.2±0.3, 6.5±0.4
17, 3.0±0.4, 2.9±0.2, 5.5±0.4
18, 4.3±0.5, 5.2±0.2, 7.0±0.4
19, 2.8±0.4, 3.2±0.3, 4.4±0.5
20, 4.8±0.3, 4.5±0.3, 5.6±0.4
21, 3.6±0.5, 3.7±0.4, 4.2±0.7
22,        , 4.1±0.4, 0.9±0.4 
23, 5.0±0.4, 5.6±0.3, 3.6±0.3
24, 4.1±0.4,        ,
25, 5.7±0.4, 5.9±0.2, 13.6±0.5
26,        , 2.9±0.2, 0.3±0.3 
27, 3.5±0.5, 4.8±0.2, 8.8±0.4
28, 2.6±0.4, 4.2±0.3, 
29, 2.6±0.4, 5.3±0.2, 
30, 3.3±0.4, 4.3±0.3, 
31, 2.7±0.3,        ,
32, 3.4±0.3, 4.0±0.2, 9.3±0.4
33, 3.2±0.2, 4.9±0.2, 
34, 4.0±0.3, 4.3±0.2, 9.5±0.4
35, 3.2±0.4, 3.7±0.3, 1.7±0.3
36, 4.8±0.4, 4.7±0.2, 
38, 5.9±0.5, 6.1±0.3, 9.9±0.4
39, 4.8±0.5, 4.7±0.3, 7.1±0.4
40, 3.0±0.4, 1.9±0.2, -8.4±0.4
41, 3.4±0.4, 4.7±0.2, 10.3±0.3
42, 3.8±0.4, 5.1±0.3, 
43, 5.9±0.3, 5.4±0.1, 
44, 5.6±0.5, 5.8±0.3, 
45, 4.7±0.5, 1.9±0.3, 
46,        , 2.7±0.1,  
47, 2.7±0.4, 4.5±0.2, 12.8±0.4
48, 2.6±0.7, 3.3±0.5, 
49, 3.7±0.5, 2.9±0.2, -6.9±0.3
50, 2.6±0.3, 2.8±0.3, 3.1±0.5
51, 1.2±0.3, 1.5±0.2, 3.9±0.4
52, 7.0±0.6,        , 
53, 3.7±0.3,        ,  2.9±0.5
54, 4.0±0.2,        ,  6.9±0.3
55, 2.7±0.3, 3.5±0.2,
56,        , 4.0±0.2,  3.4±0.3
57, 3.2±0.4, 3.8±0.3,  5.0±0.4
58, 3.9±0.4, 4.4±0.2,  5.6±0.4
59, 5.3±0.5, 4.7±0.5,  4.6±0.5
60, 4.0±0.4, 3.2±0.3,  3.0±0.4
61, 4.9±0.4, 3.2±0.2, -5.7±0.4
62, 3.4±0.3,        , -6.9±0.3 
63, 4.6±0.5, 4.3±0.2, 
64,        , 3.4±0.3,   
65, 5.5±0.5, 3.8±0.2, 
66, 2.5±0.3, 3.4±0.2, 
67, 3.2±0.4, 4.8±0.3, 
68, 2.6±0.4,10.4±0.2, 
69, 2.5±0.2, 6.6±0.2, 
70,        , 4.3±0.3,   
71, 2.9±0.4,        ,
72, 3.7±0.3, 6.7±0.2, 
73, 5.9±0.4, 6.5±0.2, 
74, 4.0±0.5, 2.7±0.3, 
75, 2.0±0.3, 1.4±0.2, 
76, 4.7±0.4, 4.1±0.2, -2.1±0.4
77, 4.7±0.5, 3.5±0.3, -0.1±0.3
78,        , 3.4±0.3, -2.9±0.3 
79,        ,        ,
80,        , 3.2±0.7, -5.2±0.6 
81, 5.0±0.4, 5.5±0.2,  8.9±0.3
82,        ,-0.4±0.3,  1.6±0.3 
83, 4.2±0.3, 3.2±0.3, -2.2±0.5
84, 4.7±0.3, 5.0±0.5,  2.5±0.3
85, 3.7±0.4, 3.2±0.3,
86,        ,        ,
87,        ,        ,
89, 2.1±0.3, 1.2±0.1,
90, 5.2±0.4, 3.2±0.1,
91, 3.6±0.4, 7.4±0.2,
92, 3.3±0.4, 8.7±0.2,
93, 3.1±0.5,        ,
94, 4.6±0.5, 8.4±0.2,
95, 2.0±0.5, 5.1±0.2,
96, 2.4±0.4, 3.0±0.2,
97, 3.1±0.4, 9.1±0.2,
98, 2.5±0.4, 7.9±0.3,
99,        , 8.2±0.2, 
101,5.4±0.4, 4.0±0.2,
102,5.6±0.5,        ,  8.5±0.4
103,4.7±0.4, 6.7±0.2,
104,       , 4.3±0.3,
105,3.6±0.4, 2.9±0.3,
106,       , 3.7±0.3, 
107,4.0±0.5,        ,
108,       , 4.4±0.2, 
109,       ,        ,
110,3.6±0.5, 2.3±0.3, -3.6±0.5
111,       , 2.9±0.2,  
112,2.1±0.5,        ,  3.6±0.2 
113,4.2±0.6, 4.2±0.3, 
114,4.1±0.4,        ,
115,4.0±0.4, 3.3±0.3,  1.9±0.5
116,4.8±0.4,        ,  2.1±0.3  
117,3.9±0.4, 4.2±0.3,  0.9±0.6
118,4.4±0.5, 4.1±0.3,  0.9±0.5
119,       , 4.0±0.2,  3.2±0.3  
121,5.0±0.8,        ,
122,2.9±0.3,        ,
123,2.1±0.3, 2.4±0.2,  2.2±0.3
124,5.2±0.4, 5.4±0.4,  5.1±1.5
125,5.4±0.5,        ,  5.2±1.5  
126,2.3±0.6,        ,  1.8±0.8  
127,3.2±0.5, 2.9±0.4,  2.9±0.7
128,4.3±0.3, 3.9±0.3,  2.6±0.5
129,5.0±0.4, 4.8±0.3,  4.4±0.6
130,3.6±0.4, 4.2±0.3,  3.2±0.6
131,4.5±0.4, 4.9±0.3,  3.4±0.5
132,5.5±0.6,        ,  1.9±0.6  
133,3.2±0.3,        ,
134,4.9±0.4,        ,  6.2±0.4  
135,       , 3.5±0.3,  1.6±0.4  
136,3.6±0.5, 3.8±0.3,  0.9±0.4
137,4.5±0.5, 4.9±0.2,  8.4±0.5
138,3.4±0.4, 3.9±0.1,  
139,3.8±0.5, 5.9±0.2,  
140,4.5±0.8, 4.5±0.3,  8.9±0.5
141,       , 6.1±0.3, 12.7±0.4
142,3.8±0.3, 4.8±0.2, 14.0±0.3
143,4.9±0.5, 6.4±0.2, 14.4±0.3
144,3.4±0.3, 4.2±0.3, 10.2±0.4
145,3.5±0.4, 3.3±0.2,  4.0±0.3
146,5.9±0.4, 6.1±0.3, 11.2±0.4
147,       , 3.9±0.2,  6.5±0.3  
148,3.2±0.4, 3.2±0.2,  3.7±0.4
149,       , 3.0±0.2,  2.3±0.3  
150,5.3±0.4, 6.2±0.3,  9.1±0.3
151,4.3±0.3,        ,  6.7±0.3  
152,3.8±0.4, 4.5±0.4,  8.6±0.6
153,0.8±0.6, 1.6±0.4,  3.6±0.6"""

vals = {}
for line in data.split("\n"):
	seq, v1, v2, v3 = line.split(',')
	val1 = v1.split('±')[0].strip()
	val2 = v2.split('±')[0].strip()
	val3 = v3.split('±')[0].strip()

	try:
		vals[int(seq)] = float(val3) - float(val1)
	except ValueError:
		pass

	
# Load the PDB file
# prot = protein.load_pdb('../data_files/4icbH_mut.pdb')
prot = protein.load_pdb('1bzrH.pdb')



ironAtom = prot[0]['A'][("H_HEM",154," ")]['FE']
met = metal.Metal(position=ironAtom.position)
met.B0 = 18.79
met.T = 273.0 + 30.0

met.iso = 30.1E-32
# met.iso = 4.4E-32
met.taur = 5.7E-9
# met.axrh = (4.24E-32, 0.0)


compare = []
for i, val in vals.items():
# for a in prot.get_atoms():
	# if a.name=='H':
		# try:
			# H = a
			# N = a.parent['N']
		# except KeyError:
			# continue

	H = prot[0]['A'][i]['H']
	N = prot[0]['A'][i]['N']

		# vec = H.position - N.position
		# dd = N.dipole_shift_tensor(H.position)
		# metvec = H.position - met.position

		# theta = np.arccos(vec.dot(metvec)/(np.linalg.norm(vec)*np.linalg.norm(metvec)))

		# pf = -(1.)/(15.)
		# pf*= met.MU0/(np.pi*4)
		# pf*= (met.B0 * H.gamma**2 * N.gamma * met.HBAR)/(np.linalg.norm(vec)**3)
		# pf*= -(3*met.iso)/(4*np.pi*np.linalg.norm(metvec)**3)
		# pf*= (3.*np.cos(theta)**2-1.)/2.
		# pf*= 4*met.taur + (3*met.taur)/(1+(met.taur*met.B0*H.gamma)**2)


	dd = N.dipole_shift_tensor(H.position)
	r2 = met.ccr_r2(H.position, H.gamma, dd)
	compare.append((r2,val))
	# print(r2/pf)



from matplotlib import pyplot as plt
x, y = zip(*compare)
plt.scatter(x, y)
plt.show()


# for atom in prot.get_atoms():
# 	if atom.name == 'H':
# 		atomH = atom
# 		res = atom.parent
# 		if 'N' in res and res.resname!="PRO":
# 			atomN = res['N']
# 			dd = atomN.dipole_shift_tensor(atomH.position)
# 			dd2 = met.dipole_shift_tensor(atomH.position)
# 			r2 = met.ccr_r2(atomH.position, atomH.gamma, dd)
# 			# print(atomH, r2)
# 			# print(met.ccr_r2_old(atomH.position, atomH.gamma, atomN.gamma,
# 			 # np.linalg.norm(atomN.position-atomH.position)))

# print(dd)
# print(dd2)

# pf = met.MU0/(4*np.pi)
# r = 6.8E-10
# d = pf*((prot[0]['A'][2]['H'].gamma**2 * met.HBAR)/(r**3))
# print(d)
# print(4*600)

# print(met.MU0)
# print(met.HBAR)