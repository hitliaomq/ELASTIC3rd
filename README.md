# ELASTIC3rd

<p align="center">
  <img src="docs/Elastic3rd_LOGO.png" alt="ELASTIC3rd" width="320"/>
</p>

**ELASTIC3rd** calculates second-, third-, fourth- and (for selected crystal
classes) fifth-order elastic constants of crystals from first-principles
codes or machine-learning interatomic potentials, by fitting — or
numerically differentiating — the energy/stress response of a crystal to a
set of homogeneous strains.

[Documentation](https://elastic3rd.readthedocs.io/en/latest) ·
[PyPI](https://pypi.org/project/elastic3rd/) (may lag behind this repo)

## Highlights (v7.0.3)

- **Three extraction schemes**
  - *strain-energy*: fit the `(E − E₀)/V₀` expansion in the Lagrangian strain;
  - *strain-stress*: fit the second Piola-Kirchhoff stress response (needs one
    order less of strain data and far fewer modes at low symmetry);
  - *numerical differentiation*: component-wise derivatives of the stress–strain
    response (`elastic3rd nd`, ref. Phys. Rev. Lett. 121, 216001 (2018)).
- **All 11 elastic crystal classes**; built-in strain-mode sets for 2nd/3rd/4th
  order, 5th order for the cubic, tetragonal, orthorhombic, monoclinic and
  triclinic classes.
- **Higher-order effect post-processing**: `elastic3rd post -n 4` (up to order 9)
  recovers the stored data and refits with extra polynomial orders
  (ref. Diamond Relat. Mater. 117, 108490 (2021)).
- **First-principles providers**: `vasp`, `vasp_pmg`, `castep`, `castep_ac`,
  `abacus`, `qe` (strain-energy only; interface not yet tested).
- **Machine-learning providers** (in-process ASE calculators, no queue needed):
  `mace`, `chgnet`, `matgl`, `sevenn`, `mattersim`, `orb`, `grace`, `tace`,
  `matris`, `nequip`, `dpmd`, `alphanet`.
- **Stability checks** for the deformed structures and a continue mode to reuse
  completed calculations.

## Installation

```bash
pip install elastic3rd          # PyPI release (may lag behind this repo)
```

or from a clone of this repository:

```bash
pip install .
```

Python 3.7 or newer on Windows or Linux. Machine-learning providers require the
corresponding force-field package (e.g. `pip install mace-torch`).

## Quick start

Put an `INPUT.e3rd` file next to your calculation inputs, e.g. for cubic Si
with VASP and the strain-energy method:

```
BaseName        Si
CrystalType     c1
Ord             3
EnergyCode      vasp
FlagSE          E
CalMode         W
GenSM           1
CoefMethod      fit
MaxStrain       5.0
STEPS           5
NP              24
```

then run

```bash
elastic3rd run        # create all strained folders and collect energy/stress
elastic3rd post       # recover SOECs and TOECs (GPa) from the data
elastic3rd post -n 4  # include the higher-order effect
elastic3rd nd         # separate numerical-differentiation method
```

Worked examples for every provider are in the
[`examples/`](https://github.com/hitliaomq/ELASTIC3rd/tree/master/examples)
and [`test/`](https://github.com/hitliaomq/ELASTIC3rd/tree/master/test)
folders; the option list is in the
[documentation](https://elastic3rd.readthedocs.io/en/latest).

## Citation

If you use ELASTIC3rd, please cite:

- M. Liao, S. Shang, Y. Liu, Z. Chen, X. Liu, H. Song, F. Xu, *ELASTIC3rd: a tool for calculating third-order elastic constants from first principles*, Comput. Phys. Commun. 261 (2021) 107777.
- M. Liao et al., *A high-efficient strain-stress method for calculating higher-order elastic constants from first-principles*, Comput. Phys. Commun. 280 (2022) 108478.
- M. Liao et al., *Revisiting the third-order elastic constants of diamond: the higher-order effect*, Diamond Relat. Mater. 117 (2021) 108490.

and, when using the numerical-differentiation method:

- T. Cao, D. Cuffari, A. Bongiorno, *First-principles calculation of  third-order elastic constants via numerical differentiation of the second Piola-Kirchhoff stress tensor*, Phys. Rev. Lett. 121 (2018) 216001.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
