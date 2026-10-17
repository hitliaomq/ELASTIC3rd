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
                                                                strain-energy method only;
                                                                interface not yet tested
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

Provider references
===================

When a calculation uses one of the providers, the corresponding publications should be
cited.

First-principles codes:

- ``vasp``: G. Kresse, J. Furthmuller, Efficient iterative schemes for ab initio
  total-energy calculations using a plane-wave basis set, Phys. Rev. B 54 (1996)
  11169, doi:10.1103/PhysRevB.54.11169.
- ``castep``, ``castep_ac``: S. J. Clark et al., First principles methods using
  CASTEP, Z. Kristallogr. - Cryst. Mater. 220 (2005) 567,
  doi:10.1524/zkri.220.5.567.65075.
- ``abacus``: K. Zhou et al., ABACUS: An electronic structure analysis package for
  the AI era, J. Chem. Phys. 163 (2025) 192501, doi:10.1063/5.0297563.
- ``qe``: P. Giannozzi et al., QUANTUM ESPRESSO: a modular and open-source software
  project for quantum simulations of materials, J. Phys.: Condens. Matter 21 (2009)
  395502, doi:10.1088/0953-8984/21/39/395502.

Machine-learning models and libraries:

- ``mace``: I. Batatia, D. P. Kovacs, C. Simm, C. Ortner, G. Csanyi, MACE: Higher
  order equivariant message passing neural networks for fast and accurate force
  fields, Adv. Neural Inf. Process. Syst. 35 (2022) 11423, arXiv:2206.07697.
- ``chgnet``: B. Deng et al., CHGNet as a pretrained universal neural network
  potential for charge-informed atomistic modelling, Nat. Mach. Intell. 5 (2023)
  1031, doi:10.1038/s42256-023-00716-3.
- ``matgl``: C. Chen, S. P. Ong, A universal graph deep learning interatomic
  potential for the periodic table, Nat. Comput. Sci. 2 (2022) 718,
  doi:10.1038/s43588-022-00349-3 (M3GNet); H. Simeon, G. de Fabritiis, TensorNet:
  Cartesian tensor representations for efficient learning of molecular potentials,
  Adv. Neural Inf. Process. Syst. 36 (2023), arXiv:2306.06482 (TensorNet); Y.-C. Ko
  et al., Materials Graph Library (MatGL), an open-source graph deep learning
  library for materials science and chemistry, npj Comput. Mater. 11 (2025),
  doi:10.1038/s41524-025-01742-y (MatGL library).
- ``sevenn``: C. Park, K. Kim, Y. Hwang, I. S. Han, Scalable parallel algorithm for
  graph neural network interatomic potentials in molecular dynamics simulations,
  J. Chem. Theory Comput. 20 (2024) 4857, doi:10.1021/acs.jctc.4c00190.
- ``mattersim``: L. Yang et al., MatterSim: A deep learning atomistic model across
  elements, temperatures and pressures, arXiv:2405.04967 (2024) (preprint).
- ``orb``: W. Neumann et al. (Orbital Materials), Orb: A fast, scalable neural
  network potential, arXiv:2410.22570 (2024).
- ``grace``: P. Bochkarev, Y. Lysogorskiy, R. Drautz, Graph atomic cluster expansion
  for semilocal interactions beyond equivariant message passing, Phys. Rev. X 14
  (2024) 021036, doi:10.1103/PhysRevX.14.021036; Y. Lysogorskiy, P. Bochkarev,
  R. Drautz, Graph atomic cluster expansion for foundational machine learning
  interatomic potentials, arXiv:2508.17936 (2025) (preprint).
- ``nequip``: S. Batzner et al., E(3)-equivariant graph neural networks for
  data-efficient and accurate interatomic potentials, Nat. Commun. 13 (2022) 2453,
  doi:10.1038/s41467-022-29939-5.
- ``dpmd``: D. Zhang, et al., A graph neural network for the era of large atomistic 
  models, Npj Comput. Mater. 12 (2026) 276. doi: 10.1038/s41524-026-02146-2.
- ``alphanet``: D. Yin et al., AlphaNet: Scaling up local-frame-based neural network
  interatomic potentials, npj Comput. Mater. 11 (2025) 332, arXiv:2501.07155.
- ``tace``: Z. Xu, W. Xie, H. Xie, P. Hu, Towards universal property prediction in
  Cartesian space: TACE is all you need, arXiv:2509.14961 (2025) (preprint).
- ``matris``: Y. Zhou, et al., MatRIS: toward reliable and efficient pretrained machine 
  learning interatomic potentials, (2026). arXiv: 2603.02002.

