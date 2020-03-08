=====================================================
chemtools: Python toolbox for Computational Chemistry
=====================================================



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
   :maxdepth: 2

   Installation <install>
   Tutorials <tutorial>
   API Reference <api>


Get involved
============

If you have some questions, remarks or requests email me at `<lmmentel@gmail.com> <mailto:lmmentel@gmail.com>`_.

Citing
======

If you use *chemtools* in a scientific publication, please consider citing the software as 

  Łukasz Mentel, *chemtools* -- A Python toolbox for computational chemistry, 2014-- . Available at: `https://bitbucket.org/lukaszmentel/chemtools <https://bitbucket.org/lukaszmentel/chemtools>`_.


Here's the reference in the `BibLaTeX <https://www.ctan.org/pkg/biblatex?lang=en>`_ format

.. code-block:: latex

   @software{chemtools2014,
      author = {Mentel, Łukasz},
      title = {{chemtools} -- A Python toolbox for computational chemistry},
      url = {https://bitbucket.org/lukaszmentel/chemtools},
      version = {0.9.1},
      date = {2014--},
  }

or the older `BibTeX <http://www.bibtex.org/>`_ format

.. code-block:: latex

   @misc{chemtools2014,
      auhor = {Mentel, Łukasz},
      title = {chemtools} -- A Python toolbox for computational chemistry, ver. 0.9.1},
      howpublished = {\url{https://bitbucket.org/lukaszmentel/chemtools}},
      year  = {2014--},
   }


Funding
=======

This project was realized through the support from the National Science Center
(Poland) grant number UMO-2012/07/B/ST4/01347.

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
