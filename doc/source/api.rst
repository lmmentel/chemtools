.. _api:

*************
API Reference
*************

Calculators
===========

.. _calculator-class:

Base Calculator
---------------

.. currentmodule:: chemtools.calculators.calculator

.. autoclass:: Calculator
   :members:

.. _dalton-class:

Dalton
------

.. currentmodule:: chemtools.calculators.dalton

.. autoclass:: Dalton
   :members:

.. _dmft-class:

DMFT
----

.. currentmodule:: chemtools.calculators.dmft

.. autoclass:: Dmft
   :members:

.. _gamessus-class:

GAMESS-US
---------

.. currentmodule:: chemtools.calculators.gamessus

.. autoclass:: GamessUS
   :members:

.. _molpro-class:

Molpro
------

.. currentmodule:: chemtools.calculators.molpro

.. autoclass:: Molpro
   :members:

.. _psi4-class:

PSI4
----

.. currentmodule:: chemtools.calculators.psi4

.. autoclass:: Psi4
   :members:

Basis Set Tools
===============

.. _basisset-class:

BasisSet module
---------------

.. currentmodule:: chemtools.basisset

.. autoclass:: BasisSet
   :members:


.. autofunction:: has_consecutive_indices

.. autofunction:: reorder_shell_to_consecutive

.. autofunction:: merge

.. autofunction:: primitive_overlap

.. autofunction:: nspherical

.. autofunction:: ncartesian

.. autofunction:: zetas2legendre

.. autofunction:: zetas2eventemp

.. autofunction:: eventemp

.. autofunction:: welltemp

.. autofunction:: legendre

.. autofunction:: xyzlist

.. autofunction:: zlmtoxyz

.. _basisparse-module:

basisparse module
-----------------

Parsing basis set from different formats.

.. automodule:: chemtools.basisparse
   :members:


basisopt module
---------------

Optimization of basis set exponents and contraction coefficients.

.. automodule:: chemtools.basisopt
   :members:

CBS module
----------

Complete Basis Set (CBS) extrapolation techniques.

.. automodule:: chemtools.cbs
   :members:

Molecule module
---------------

.. automodule:: chemtools.molecule
   :members:

parsetools module
-----------------

.. automodule:: chemtools.parsetools
   :members:

submitgamess module
-------------------

.. automodule:: chemtools.submitgamess
   :members:

submitmolpro module
-------------------

.. automodule:: chemtools.submitmolpro
   :members:

gamessorbitals module
---------------------

.. automodule:: chemtools.calculators.gamessorbitals
   :members:

gamessreader module
-------------------

.. automodule:: chemtools.calculators.gamessreader
   :members:
