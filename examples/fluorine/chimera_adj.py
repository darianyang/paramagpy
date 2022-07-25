import os
from examples.fluorine.pcs_fit_test_dchi import calc_racs
import subprocess
import numpy as np

# use chimera generate new pdb files with 43 rotation
fo = open("chimera.cmd", "w")
fo.write("open 0_deg_rot_gb1-ntaco-leap_adj.pdb\n")
fo.write("rotation 1 :43.A@CA,CB\n")
fo.write("display :43\n")
for i in range (0,360,18):
    fo.write(f"rotation 1 {i} \n")
    fo.write(f"write format pdb 0 {i}.pdb \n")
fo.close()

p = subprocess.Popen("powershell C:\Program` Files\Chimera` 1.16\\bin\chimera.exe chimera.cmd")
#p.wait()
#p.kill()
#subprocess.call(['kill', str(p.pid)])
#p.terminate()
#p.wait()
#########################

#########################
racs = []
for i in range (0,360,18):
    foo = open(f"{i}.pdb", "r")
    lines = foo.readlines()
    foo.close()
    fooo = open(f"{i}_adj.pdb", "w")
    
    for line in lines:
        if line.startswith("HETATM"):
            fooo.write(line.replace("HETATM", "ATOM  "))
        else:
            fooo.write(line)
    racsval = calc_racs(f"{i}_adj.pdb")
    racs.append(racsval)
    
    fooo.close()
    os.remove(f"{i}.pdb")
    os.remove(f"{i}_adj.pdb")

np.savetxt("racs.txt",racs)




# racs new pdb test

