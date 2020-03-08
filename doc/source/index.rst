=====================================================
chemtools: Python toolbox for Computational Chemistry
=====================================================

Chemtools is a set of modules that is intended to help with more
advanced computations using common electronic structure programs.

The main to goal was to enable convientnt `basis set <https://en.wikipedia.org/wiki/Basis_set_%28chemistry%29>`_ manipulations, including designing and optimizing exponents of basis sets. To acheive that there are several modules abstracting various functionalities:

Basis set object 
  :mod:`basisset <chemtools.basisset>` module that contains the :class:`chemtools.basisset.BasisSet`. See the :doc:`tutorial` page for a overview
  of the capabilities.

Basis set optimization
  :mod:`basisopt <chamtools.basisopt>` module that contains functionalities
  for building up basis set optimization tasks. See the :doc:`tutorial`
  page for a few illustrative examples.

Calculators
  :mod:`chemtools.calculators.calculator` contains wrappers for several packages
  handling the actual electronic structure calcualtions. Currently there is support for:

  * Dalton_
  * Gamess-US_
  * MolPro_
  * PSI4_

Molecule
  :mod:`chemtools.molecule` is general purpose module intorducing molecule class
  for handling molecules and defining atomic compositions and molecular geometries

Complete basis set extrapolation
  :mod:`chemtools.cbs` contains functions for performing complete basis set (CBS) extrapolation using different approaches.


Contents
========

.. toctree::
   :maxdepth: 2

   Installation <install>
   Tutorials <tutorial>
   API Reference <api>


Contributing
============

* `Source <https://github.com/lmmentel/chemtools>`_
* `Report a bug <https://github.com/lmmentel/chemtools/issues>`_
* `Request a feature <https://github.com/lmmentel/chemtools/issues>`_
* `Submit a pull request <https://github.com/lmmentel/chemtools/pulls>`_

Contact
=======

Łukasz Mentel 

*  github: `lmmentel <https://github.com/lmmentel>`_
*  email: lmmentel <at> gmail.com


Citing
======

If you use *chemtools* in a scientific publication, please consider citing the software as 

  Łukasz Mentel, *chemtools* -- A Python toolbox for computational chemistry, 2014-- . Available at: `https://bitbucket.org/lukaszmentel/chemtools <https://bitbucket.org/lukaszmentel/chemtools>`_.


Here's the reference in the `BibLaTeX <https://www.ctan.org/pkg/biblatex?lang=en>`_ format

.. code-block:: latex

   @software{chemtools2014,
      author = {Mentel, Łukasz},
      title = {{chemtools} -- A Python toolbox for computational chemistry},
      url = {https://github.com/lmmentel/chemtools},
      version = {0.9.1},
      date = {2014--},
  }

or the older `BibTeX <http://www.bibtex.org/>`_ format

.. code-block:: latex

   @misc{chemtools2014,
      auhor = {Mentel, Łukasz},
      title = {{chemtools} -- A Python toolbox for computational chemistry, ver. 0.9.1},
      howpublished = {\url{https://github.com/lmmentel/chemtools}},
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


.. _DALTON: http://daltonprogram.org/
.. _Gamess-US: http://www.msg.ameslab.gov/gamess
.. _MolPro: http://www.molpro.net/
.. _PSI4: http://www.psicode.org/
