============
Get Started
============


Simplest way to run
===================

- Prepare the `INPUT` file and the input file for first principlec code
- Enter the folder
- Simply run ``elastic3rd run``


Inputs
======

INPUT
-----

===========  =======  =======================================================================================
Parameters   Default  Comments
===========  =======  =======================================================================================
BaseName     Si       The folder name of the result. In addition, 
                        For CASTEP, it must equals to the name of *cell* and *param* file
                        For VASP, no above required
CrystalType  Cubic1   The symmetry of the crystal. It should be one of the followings
                        triclinic or n, monoclinic or m, orthorhombic or o, tetragonal1 or t1, 
                        tetragonal2 or t2, rhombohedral1 or r1, rhombohedral2 or r2,
                        hexagonal1 or h1, hexagonal2 or h2, cubic1 or c1, cubic2 or c2
Ord          3        The order of elastic constants
EnergyCode   CASTEP   The first principles code, currently support CASTEP and VASP
EnergyRun    1        The way to run the first principles code
                        1 is the default way. For CASTEP: RunCASTEP -np NP BaseName
                                              For VASP: mpirun -np NP vasp_std
                        0 means user defined. It will read the *energyrun* file in elastic3rd/energy folder.
                          And it can be overwrite by elastic3rd.energy.glue.write_energyrun(RunStr) function
MaxStrain    5        The maximum strain, in unit of %. e.g. 5 means the strain range is -5% to 5%
STEPS        5        The number of points in each strain mode
NP           1        The number of cores used in the first principles codes
===========  =======  =======================================================================================

First principles inputs
-----------------------
- CASTEP
  For CASTEP, two files are required for calculating, *BaseName.cell* and *BaseName.param*. For details, ref. `CASTEP files`_ .

.. _`CASTEP files`: http://www.tcm.phy.cam.ac.uk/castep/documentation/WebHelp/content/modules/castep/expcastepfileformats.htm

- VASP
  To run VASP, four files are required: *INCAR*, *POSCAR*, *KPOINTS* and *POTCAR*. For details, ref. `VASP files`_ .

.. _`VASP files`: https://www.vasp.at/wiki/index.php/Category:Input_Files


Outputs
=======