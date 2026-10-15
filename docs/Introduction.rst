============
Introduction
============

Elastic3rd is an open-source tool to calculate the higher-order elastic constants (HOECs)
of materials from first-principles calculations or machine-learning force fields.
It generates symmetry-adapted strain modes for a crystal, deforms the structure over a
strain range, collects the energy or stress response from an external calculator, fits the
response, and solves the resulting symmetry equations for the independent elastic
constants of any requested order.

Features
--------
- Second-, third- and fourth-order elastic constants for all 11 elastic crystal classes
  (triclinic, monoclinic, orthorhombic, tetragonal 1/2, rhombohedral 1/2, hexagonal 1/2,
  cubic 1/2); fifth-order is partially supported (cubic, tetragonal,
  orthorhombic, monoclinic, triclinic)
- Both the strain-energy method and the strain-stress method (higher orders require
  substantially fewer deformation modes with the strain-stress method)
- Two coefficient-recovery schemes: polynomial fitting and numerical differentiation
  (arbitrary accuracy stencils, 2/4/6/8)
- Post-processing that takes higher-order effects into account when recovering the
  constants from the fitted coefficients (M. Liao et al., Revisiting the third-order
  elastic constants of diamond: the higher-order effect, Diamond Relat. Mater. 117
  (2021) 108490, doi:10.1016/j.diamond.2021.108490)
- First-principles providers: CASTEP (Materials Studio or standalone, Windows and Linux),
  VASP, a pymatgen-based VASP interface, ABACUS and Quantum ESPRESSO (energy only)
- Machine-learning providers executed in-process (no queue system needed): MACE, CHGNet,
  MatGL (M3GNet / TensorNet), SevenNet, MatterSim, ORB, GRACE, NequIP, DeepMD (DP),
  AlphaNet, TACE and Matris
- Command line interface (``elastic3rd run`` / ``elastic3rd post``) and a plain
  script-based API (``elastic3rd.elastic.elastic3``)
- Input files accepted as plain text, JSON or YAML; runs can be resumed with Continue mode
- Python 3.7+ on Windows and Linux

Requirements
------------
python, numpy, scipy, matplotlib, pyyaml.

Optional dependencies for specific providers: pymatgen (VASP-pymatgen interface, ABACUS,
and reading CIF/POSCAR structures for ML providers), ASE and the corresponding
framework package for each machine-learning provider (see :doc:`EnergyCodes`).
