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
    r^{'} = Fr

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
the first-principles code into the second Piola-Kirchhoff stress (``cauchy2pk``). In the 
strain-stress method, one deformation can generate more independent equation about elastic 
consatnts (up to 6 if the strain mode is chosen carefully), which means the strain-stress 
method (``FlagSE S``) needs less strain modes compared to strain-energy method (ref. Comput. 
Phys. Commun. 
280 (2022) 108478).

Numerical differentiation method
================================

The separate numerical-differentiation scheme of Phys. Rev. Lett. 121, 216001 (2018)
(Cao *et al.*) treats every independent constant on its own, and is driven by the
``elastic3rd nd`` command (see :doc:`examples`). Instead of fitting the response of a
complete mode set with one polynomial, the constants are obtained directly as
derivatives of the second Piola-Kirchhoff stress at zero strain:

.. math::
    C_{ijkl}=\left.\frac{\partial\sigma_{ij}}{\partial\eta_{kl}}\right|_{\eta=0},
    \qquad
    C_{ijklmn}=\left.\frac{\partial^{2}\sigma_{ij}}
    {\partial\eta_{kl}\,\partial\eta_{mn}}\right|_{\eta=0}.

For each constant, the program generates a dedicated set of strain modes acting only on
the components that survive the permutation symmetry of that constant (Eq. 6-8 of the
reference, implemented in ``elastic3rd.tools.sep_nd``), so the requested derivative is
read off the conjugate stress component (converted from the Cauchy stress by
``cauchy2pk``) without solving any coupled system.

The derivatives are estimated with finite-difference stencils symmetric around
:math:`\eta = 0`. The accuracy parameter ``Acc`` sets the number of points:
:math:`\pm\varepsilon` for ``Acc 2``, :math:`\pm\varepsilon/2` and :math:`\pm\varepsilon`
for ``Acc 4``, and so on up to 8-point stencils. The stencil weights are generated for
the actual strain list by solving the linear system that makes the formula exact for
polynomials up to the stencil order (``elastic3rd.tools.nd_coef``); when energy rather
than stress is differentiated, the vanishing first derivative at equilibrium is imposed
as an extra condition, which reduces the number of points needed for the same accuracy.

Because no global polynomial model is assumed, each estimate uses the fewest possible
strain points and applies to every crystal class without a fitted mode set; the price
is that each constant needs its own calculations, and a bare central difference
inherits the numerical noise of the stress. Elastic3rd therefore evaluates every
constant at each requested strain amplitude, writes the individual estimates to
``EC.json`` and their average to ``EC_ave.json``, and prints the averaged constants in
GPa; the spread over amplitudes is a direct, assumption-free error estimate. The method
requires stress output, so it works with ``vasp``, ``castep`` and ``abacus``.

Coefficient relation and higher-order effects
=============================================

In both cases, the symmetry relations among the tensor components (Brugger for 2nd/3rd,
Brendel for 4th) are encoded in ``elastic3rd/symmetry``: each strain mode contributes a 
linear equation that relates the fitted coefficients to the independent constants, and 
the system is solved exactly when the mode set is complete, or by least squares otherwise.

Fully supported independent-constant tables (and built-in strain modes) exist for 2nd,
3rd and 4th order over all 11 elastic crystal classes.

Truncating the expansion at the order being fitted biases higher-order constants when
the strain window is large (ref. https://doi.org/10.1016/j.diamond.2021.108490). The 
``elastic3rd post`` command refits the stored data with an extra polynomial order and 
reports the constants after the higher-order terms are taken into account (``Flag_Ord`` 
up to 9), together with the convergence of each constant versus the strain window.


