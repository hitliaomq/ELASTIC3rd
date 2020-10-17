=================
Theory Background
=================

The strain-energy method to calculate TOECs is based on the conituumn elasticity theory. The base equation is as follows (note: the Einstein summation notation is used in the following equations):

|eq1|

Where

- :math:`E` is the energy of the deformed structure

- :math:`E_0` is the energy of the initial structure

- :math:`V_0` is the volume of the initial structure

- :math:`C_{ijkl}` is SOECs

- :math:`C_{ijklmn}` is TOECs

- :math:`\eta_{ij}` is the Lagrangian strain

And the Lagrangian strain can be written as follows:

|eq2|

The F is the deformation gradient, and can be expressed by the lattic vector of deformed structure(r') and intial structure (r).

|eq3|

when using symmetrical strain, the deformation gradient can be expressed by the Lagrangian strain.

|eq4|

where :math:`y_i` and :math:`\lamda_i` are the eigenvector and eigenvalue of **Y**.

For each strain mode, we only adjust the amplitude. Hence, we can express the strain as a function of strain amplitude :math:`\eta`, then we have

|eq5|

Finally, we can choose different strain modes, then get different equations about TOECs and SOECs. By solving the equations, we can get the value of SOECs and TOECs.

.. |eq1| image:: Eq-1.png
.. |eq2| image:: Eq-2.png
.. |eq3| image:: Eq-3.png
.. |eq4| image:: Eq-4.png
.. |eq5| image:: Eq-5.png