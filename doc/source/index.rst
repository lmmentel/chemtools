===================================================
ChemTools: Python tools for Computational Chemistry
===================================================

Chemtools is a set of modules that is intended to help with more
advanced computations using common electronic structure methods/
programs. Currently the is some limited support for Gamess-US_ and
MolPro_ program packages but other codes can be easily interfaced.

.. _DALTON: http://daltonprogram.org/
.. _Gamess-US: http://www.msg.ameslab.gov/gamess
.. _MolPro: http://www.molpro.net/
.. _PSI4: http://www.psicode.org/

Current Modules
===============

* :mod:`basisopt`: module for optimizing one electron basis function exponents
* :mod:`basisset`: module for handling basis set in different formats and obtaining
  exponents from series like: even tempered, well tempered, Legendre expansion,
* :mod:`molecule`: general purpose module intorducing molecule class for handling
  molecules


Contents
========

.. toctree::
   :caption: Table of contents
   :maxdepth: 2

   Installation <install>
   Tutorials <tutorial>
   API Reference <_reference/modules>


Get involved
============

If you have some questions, remarks or requests email me at `<lmmentel@gmail.com> <mailto:lmmentel@gmail.com>`_.

Related projects
================

`cclib <http://cclib.github.io/>`_
    A python library for parsing 

`pygamess <https://github.com/kzfm/pygamess>`_
    Python wrapper for Gamess(US)


License
=======

.. include:: ../../LICENSE.rst


Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
