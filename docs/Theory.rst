=================
Theory Background
=================

The strain-energy method to calculate TOECs is based on the conituumn elasticity theory. The base equation is as follows (note: the Einstein summation notation is used in the following equations):

.. math::
    $$E=E_{0}+\frac{1}{2 !} V_{0} C_{i j k} \eta_{i j} \eta_{H}+\frac{1}{3 !} V_{0} C_{i j k l m n} \eta_{i j} \eta_{k l} \eta_{m n}$$

Where :
- *E* is the energy of the deformed structure
- *E_0* is the energy of the initial structure
- *V_0* is the volume of the initial structure
- C_{ijkl} is SOECs
- C_{ijklmn} is TOECs
- eta_{ij} is the Lagrangian strain

And the Lagrangian strain can be written as follows:

.. math::
    \eta_{ij} = \frac{1}{2}(F_{ki} F_{kj} - \delta_{ij})

The F is the deformation gradient, and can be expressed by the lattic vector of deformed structure(r') and intial structure (r).

.. math::
    r' = Fr