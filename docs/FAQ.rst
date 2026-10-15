============
FAQ
============

Which crystal type should I set for ``CrystalType``?
    Use the elastic (Laue) class of the crystal, not the space group: i, c1, c2, h1, h2,
    r1, r2, t1, t2, o, m, n (long names such as ``cubic1`` also work). The 32 crystallo-
    graphic point groups collapse into these classes; helpers in
    ``elastic3rd.symmetry.symdata`` (``group2crytyp``, ``space_group_number2CrystalType``)
    can map a point group or space group for you.

My VASP stress run returns nothing / wrong constants.
    For the strain-stress method (``FlagSE S``) the stress must be written to ``OUTCAR``:
    set ``ISIF = 2`` in the INCAR. The ``INCAR`` written by the VASP provider for ion
    relaxation already does this; a hand-written INCAR must.

What are the units?
    Strain amplitudes in ``MaxStrain``/``MinStrain``/``STRAINAMP`` are given in percent.
    Energies are eV per cell (converted also to Hartree, kJ/mol and kcal/mol in the
    outputs), stresses GPa, and all reported elastic constants GPa.

Can I use my own strain modes?
    Yes. Set ``CalMode S`` and provide a ``STRAINMODE`` file with one Voigt 6-vector
    (``xx yy zz yz xz xy``) per line. With ``CalMode W, GenSM 1`` the built-in modes are
    written to a ``STRAINMODE`` file, which you can inspect and then reuse.

Can I control the exact strain amplitudes?
    Yes. Put a file named ``STRAINAMP`` in the working folder containing the amplitudes
    (in percent); it overrides ``MinStrain``, ``MaxStrain``, ``STEPS`` and ``StrainType``.
    For ``CoefMethod num`` the list must be symmetric around 0 (e.g. ``-5 -2.5 0 2.5 5``)
    with ``STEPS = 5`` and the desired ``Acc``.

The job was interrupted. Must I start over?
    No. Set ``Continue 1`` in ``INPUT.e3rd`` and run again from the same folder: existing
    ``Energy_*.txt``/``Stress_*.txt`` files are read instead of recalculated. New strain
    points are still computed.

My relaxed structures at each strain drift / symmetry is not kept.
    Set ``StableEst 1`` (relax ions at fixed cell) if you want the lattice to stay exactly
    as deformed; relaxions with full cell freedom can drift. ``AddPert 1`` helps when ion
    relxations stall, and ``Geom`` handles the pre-deformation equilibrium structure.

How do I get 4th-order constants?
    Set ``Ord 4`` (built-in mode tables and coefficient relations exist for all 11 classes
    at 2nd/3rd/4th order). For 5th order, only c1/c2/t1/t2/o/m/n are tabulated. You can
    also always recover higher-order information afterwards from stored data with
    ``elastic3rd post -n N``.

How do I choose a different model for a machine-learning provider?
    Add a ``settings.yaml`` in the working folder with ``model``, ``device``, ``optimizer``,
    ``ase_filter``, ``steps`` keys (the accepted model names are the ones of the
    corresponding upstream package, e.g. ``M3GNet-MP-2021.2.8-PES`` for ``matgl``).
    Locally trained NequIP/DeepMD/AlphaNet models are read from their default cache
    locations or via ``model_path``/``weights_path``.

Can I run Elastic3rd inside a Jupyter notebook or a batch script?
    Yes, the CLI is a thin wrapper around ``elastic3rd.elastic.elastic3(INPUT="...")``.
    Copy the pattern of ``examples/runElastic3rd.py``.
