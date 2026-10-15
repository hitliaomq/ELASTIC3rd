============
Installation
============

Requirements
============
- Python 3.5 or newer
- numpy, scipy, matplotlib, pyyaml (installed automatically by pip)
- Optional, depending on the calculator you use:

  - pymatgen: required by the ``vasp_pmg`` interface, the ABACUS interface and by all
    machine-learning providers (to read CIF/POSCAR structures)
  - ASE plus the provider package itself (``mace-torch``, ``chgnet``, ``matgl``,
    ``sevenn``, ``mattersim``, ``orb-models``, ``grace-atomistic-models``, ``nequip``,
    ``deepmd-kit``, ``alphanet``, ``tace``, ``matris``, ...) for machine-learning runs
  - abacuspy: required by the ``abacus`` interface

Install by pip
==============
::

    pip install elastic3rd

This also installs the ``elastic3rd`` command line entry point::

    elastic3rd -h

Install from source
===================
Clone or download the repository, enter the project root folder (the folder that
contains ``setup.py``), and run::

    python setup.py install

or, for a development install so that code edits take effect immediately::

    pip install -e .

After installation you can check that everything works by importing the package::

    python -c "import elastic3rd; print(elastic3rd.__version__)"

The machine-learning providers are imported lazily, so their heavy dependencies do not
need to be installed unless you use them.
