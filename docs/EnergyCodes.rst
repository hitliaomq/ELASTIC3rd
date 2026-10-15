===========================
Energy and Stress Providers
===========================

Elastic3rd does not compute energies or stresses itself. The value of ``EnergyCode`` in
``INPUT.e3rd`` selects a provider module ``elastic3rd.energy.<EnergyCode>``
(case-insensitive). Providers fall into two families:

- **First-principles codes**, launched as external shell commands per strain folder;
- **Machine-learning force fields**, executed in-process through ASE calculators
  (listed in ``elastic3rd.elastic.RUNDIRECTCODE``), which need no run command and no
  queue system.

First-principles providers
==========================

===============  =========  ==================================  ====================================================
EnergyCode       E/S        Files in the working folder         Notes
===============  =========  ==================================  ====================================================
``vasp``         E + S      ``INCAR``, ``POSCAR``,              Energy from ``OSZICAR``, stress
                            ``KPOINTS``, ``POTCAR``             from ``OUTCAR`` (requires
                                                                ``ISIF = 2``); runs
                                                                ``mpirun -np NP vasp_std``
``vasp_pmg``     E + S      ``POSCAR``                          pymatgen-based; energy from
                                                                ``OUTCAR``, stress from
                                                                ``vasprun.xml``
``castep``       E + S      ``BaseName.cell``,                  Results read from
                            ``BaseName.param``,                 ``BaseName.geom``, else
                            ``RunCASTEP.bat/sh``                ``BaseName.castep``; VCA supported
``castep_ac``    E + S      same as ``castep``                  older CASTEP variant
``abacus``       E + S      ``STRU``, ``KPT``, ``INPUT``        via ``abacuspy``; results from
                                                                ``abacus.json``; runs
                                                                ``mpirun -np NP abacus``
``qe``           E only     ``BaseName.scf``, reads             Energy from ``BaseName.out``;
                                                                strain-energy method only
===============  =========  ==================================  ====================================================

Machine-learning providers
==========================

For ML providers, place a structure file readable by pymatgen next to ``INPUT.e3rd``:
``BaseName.cif``, ``BaseName.poscar``, ``BaseName.json``, ``BaseName.yaml`` or a plain
``POSCAR``. The provider relaxes the structure (or runs a single point) with the ASE
calculator of the force field, then writes the ``ENERGY`` (eV) and ``STRESS`` (GPa) files
that Elastic3rd reads.

==============  ============================  ==========================================
EnergyCode      Package to install            Default model
==============  ============================  ==========================================
``mace``        ``mace-torch``                ``mace-matpes-pbe-0``
``chgnet``      ``chgnet``                    CHGNet pretrained
``matgl``       ``matgl``                     ``TensorNet-MatPES-PBE-v2025.1-PES``
                                              (also serves M3GNet models, e.g.
                                              ``M3GNet-MP-2021.2.8-PES``)
``sevenn``      ``sevenn``                    ``7net-mf-ompa``
``mattersim``   ``mattersim``                 ``MatterSim-v1.0.0-1M``
``orb``         ``orb-models``                ``orb-v3-conservative-inf-omat``
``grace``       ``grace-atomistic-models``    ``GRACE-2L-OMAT-large-ft-AM``
``nequip``      ``nequip``                    compiled model under ``~/.cache/nequip/``
``dpmd``        ``deepmd-kit``                frozen model under ``~/.cache/deepmd/``
``alphanet``    ``alphanet``                  checkpoint under ``~/.cache/alphanet/``
``tace``        ``tace``                      ``TACE-v1-OAM-M``
``matris``      ``matris``                    ``matris_10m_oam``
==============  ============================  ==========================================

The model, device and optimizer are configured in an optional ``settings.yaml`` placed
in the working folder, e.g.::

    model: M3GNet-MP-2021.2.8-PES
    device: cpu
    optimizer: FIRE
    ase_filter: ExpCellFilter
    steps: 400

Common keys: ``model``, ``model_path``/``weights_path`` (for locally trained models such
as NequIP/DP/AlphaNet), ``device`` (``cpu``/``cuda``), ``optimizer``
(BFGS/FIRE/LBFGS/FIRE2), ``ase_filter`` (ExpCellFilter/FrechetCellFilter/
StrainFilter/UnitCellFilter), ``steps``, ``fmax``, ``dtype``/``precision``
(sevenn/tace). Elastic3rd itself maintains the ``flageg`` and ``constrain_symmetry``
keys: ``flageg`` encodes what to do at each step (``EO`` single point, ``E`` relax ions
at fixed cell, ``G`` relax cell and ions, ``C`` cell only), matching the ``Geom`` and
``StableEst`` input parameters.

Energy and stress are taken from the final relaxed configuration
(``get_potential_energy``, ``get_stress``). Because the calculators run in-process, the
``EnergyRun``/``NP`` input values are ignored for this family.

Choosing a provider for validation work
=======================================

ML force fields make full HOEC pipelines (many strain points x several modes) very
cheap, which is convenient for testing an input setup before spending CPU-time on a
first-principles production run. Note, however, that the elastic constants inherit the
accuracy of the underlying model; report final values from the intended electronic
structure method.
