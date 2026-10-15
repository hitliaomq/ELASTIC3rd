==============
ELASTIC3rd API
==============

The public API is function-based (there is no ``Elastic3rd`` class). The entry points
most users need:

````elastic3rd.elastic````
    ``elastic3(INPUT="INPUT.e3rd")`` -- the full workflow described in Get Started
````elastic3rd.esutils````
    input parsing (``read_input``), strain list generation (``gen_strain_list``), strain-mode file IO (``read_strainmode``), logging, crystal type mapping
````elastic3rd.crystal.lattice````
    lattice utilities: ``print_lattice`` (lattice parameters and volume)
````elastic3rd.crystal.deform````
    strain to deformation gradient (``strain2deformgrad``)
````elastic3rd.symmetry````
    symmetry coefficient tables and mode generation: ``symmetry`` (energy method), ``symmetry_stress`` (stress method), ``symdata`` (independent-constant relations), ``gencoef`` (built-in strain-mode library), ``print_cijk`` (label results)
````elastic3rd.energy````
    one provider module per ``EnergyCode`` plus the shared helpers ``glue`` (units, run command) and ``ml_ff`` (base class for the ML providers)
````elastic3rd.post````
    ``post``/``post_mode``/``post_single`` (refit collected data, higher-order effects), ``esfit`` (polynomial fits and figures), ``post_nd`` (general numerical-differentiation post)
````elastic3rd.tools````
    auxiliary scripts: ``runElastic3rd.py`` (workflow launcher to copy into a work folder), ``runAll.py`` (pymatgen/dfttk-based VASP input generation and batch runs), ``nd_coef`` (finite-difference weights), ``shockwave`` (shock compression analysis), ``sep_nd`` (component-wise SOEC/TOEC by ND), ``post_asym`` (convergence vs strain-window study)
````elastic3rd.scripts````
    ``run_elastic3rd`` -- the ``elastic3rd`` console command

Typical programmatic use:

.. code::

    import elastic3rd.elastic as e3
    e3.elastic3(INPUT="INPUT.e3rd")

Post-processing only, with error estimates and R^2 from the weighted fit:

.. code::

    import elastic3rd.post.post2 as espost2
    C, Cerr = espost2.post_single(x, E, StrainIn, CoefMethod="fit",
                                  V0=V0, Flag_Ord=3, INPUT="INPUT.e3rd")

Provider API (to plug in a new code) is documented in :doc:`GetStart`.

A full autodoc reference is planned (TODO).
