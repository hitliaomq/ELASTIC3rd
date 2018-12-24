#!python
#

import energy.vasp
import energy.castep
import energy.glue
import esutils
import os

if os.path.exists("Al"):
    pass
else:
    os.mkdir("Al")
energy.vasp.copy_files("Al", "Al")

BaseName = "Al"
BaseVec = energy.vasp.get_base_vec(BaseName)
print BaseVec
BaseVec[0][2] = 3.0
energy.vasp.write_base_vec(BaseName, BaseVec)

E0 = energy.vasp.get_energy(BaseName)
print E0

RunStr = energy.vasp.run(4, BaseName)
#print RunStr

energy.glue.write_energyrun(RunStr)
energy.glue.run()

BaseName = "Si"
energy.castep.copy_files(BaseName, "Si")

ParaIn = esutils.read_input()
print ParaIn
