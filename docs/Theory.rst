=================
Theory Background
=================

Strain-energy method
====================

The strain-energy method to calculate TOECs is based on the continuum elasticity theory.
The base equation is as follows (note: the Einstein summation notation is used in the
following equations):

.. math::
    E=E_0+\frac{1}{2!} V_0 C_{ijkl} \eta_{ij} \eta_{kl} + \frac{1}{3!} V_0 C_{ijklmn} \eta_{ij} \eta_{kl} \eta_{mn} + ...

Where

- :math:`E` is the energy of the deformed structure

- :math:`E_0` is the energy of the initial structure

- :math:`V_0` is the volume of the initial structure

- :math:`C_{ijkl}` is SOECs

- :math:`C_{ijklmn}` is TOECs

- :math:`\eta_{ij}` is the Lagrangian strain

And the Lagrangian strain can be written as follows:

.. math::
    \eta_{ij} = \frac{1}{2} (F_{ki} F_{kj} - \delta_{ij})

The F is the deformation gradient, and can be expressed by the lattic vector of deformed structure(r') and intial structure (r).

.. math::
    r^' = Fr

when using symmetrical strain, the deformation gradient can be expressed by the Lagrangian strain.

.. math::
    F = Q \sqrt{\lambda} Q^T

    Y = 2\eta + I

    Q = (y_1, y_2, y_3)

    \lambda = diag(\lambda_1, \lambda_2, \lambda_3)

where :math:`y_i` and :math:`\lambda_i` are the eigenvector and eigenvalue of **Y**, I is the identity matrix.

For each strain mode, we only adjust the amplitude. Hence, we can express the strain as a function of strain amplitude :math:`\eta`, then we have

.. math::
    \frac{E-E_0}{V_0} = \frac{1}{2} A_2 \eta^2 + \frac{1}{6} A_3 \eta^3 + O(\eta^4)

Finally, we can choose different strain modes, then get different equations about TOECs and SOECs. By solving the equations, we can get the value of SOECs and TOECs.

Strain-stress method
====================

Differentiating the energy expansion above with respect to strain gives the conjugate
(second Piola-Kirchhoff) stress:

.. math::
    \sigma_{ij} = \frac{1}{V_0}\frac{\partial E}{\partial \eta_{ij}}
    = \frac{1}{2!} C_{ijkl} \eta_{kl} + \frac{1}{3!} C_{ijklmn} \eta_{kl} \eta_{mn} + ...

so fitting stress vs. strain yields constants of order :math:`N` from a polynomial of
degree :math:`N-1` instead of :math:`N`. Elastic3rd converts the Cauchy stress printed by
the first-principles code into the second Piola-Kirchhoff stress (``cauchy2pk``), which
means the strain-stress method (``FlagSE S``) needs one order less of strain data and
considerably fewer deformation modes for low-symmetry crystals (ref. Phys. Rev. Lett.
121, 216001 (2018)). The pre-factors differ accordingly: an :math:`N`-th order constant
is recovered from a fitted coefficient multiplied by :math:`N!` in the energy method and
:math:`(N-1)!` in the stress method.

Coefficient recovery
====================

Two schemes are available (``CoefMethod``):

- **fit**: the energy (or stress) of each strain mode is fitted with a polynomial whose
  degree equals the requested order (starting at :math:`\eta^2` for energy,
  :math:`\eta^1` for stress, since the unstrained structure is at equilibrium);
- **num**: the derivatives at zero strain are estimated directly by finite-difference
  stencils of accuracy 2/4/6/8 (``Acc``) on a symmetric strain list, following the
  separated-strain scheme of Phys. Rev. Lett. 121, 216001 (2018).

In both cases, the symmetry relations among the tensor components (Brugger for 2nd/3rd,
Brendel for 4th, and higher orders where tabulated) are encoded in
``elastic3rd/symmetry``: each strain mode contributes a linear equation that relates the
fitted coefficients to the independent constants, and the system is solved exactly when
the mode set is complete, or by least squares otherwise.

Higher-order effects
====================

Truncating the expansion at the order being fitted biases higher-order constants when
the strain window is large (ref. https://doi.org/10.1016/j.diamond.2021.108490). The ``elastic3rd post`` command refits the stored data with
an extra polynomial order and reports the constants after the higher-order terms are
taken into account (``Flag_Ord`` up to 9), together with the convergence of each
constant versus the strain window.

Supported orders
================

Fully supported independent-constant tables (and built-in strain modes) exist for 2nd,
3rd and 4th order over all 11 elastic crystal classes; 5th order is available for the
cubic, tetragonal, orthorhombic, monoclinic and triclinic classes.
