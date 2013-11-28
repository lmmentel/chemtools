ChemTools: Python tools for Computational Chemistry
===================================================

.. toctree::
   :maxdepth: 2

Overview
========
Chemtools is a set of modules that is intended to help with more
advanced computations using common electronic structure methods/
programs. Currently the is some limited support for Gamess-US_ and
MolPro_ program packages.

.. _Gamess-US: http://www.msg.ameslab.gov/gamess
.. _MolPro: http://www.molpro.net/

Current Modules
---------------

* basisopt: module for optimizing one electron basis function exponents
* basisset: module for handling basis set in different formats and obtaining
            exponents from series like: even tempered, well tempered, Legendre
            expansion,
* dmft: utility functions for running DMFT calculations,
* gamessus: utility function for parsing input and log files from Gamess-US
  calculations, wrappers for running the code from python
* molpro: parser for the output file
* molecule: general purpose module intorducing moelcule class for handling
  molecules

License
-------
ChemTools are release under GPLv3.0.

Download and Installation
=========================

Prerequisites
-------------

Python_ version 2.7.3 or later
Periodic_ package
Scipy_ verison 0.11 or later

.. _Python: http://python.org/
.. _Periodic: https://pypi.python.org/pypi/periodic
.. _SciPy: http://www.scipy.org/
