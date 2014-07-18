ChemTools: Python tools for Computational Chemistry
===================================================


Chemtools is a set of modules that is intended to help with more
advanced computations using common electronic structure methods/
programs. Currently the is some limited support for Gamess-US_ and
MolPro_ program packages but other codes can be easily interfaced.

.. _Gamess-US: http://www.msg.ameslab.gov/gamess
.. _MolPro: http://www.molpro.net/



Current Modules
---------------

* :mod:`basisopt`: module for optimizing one electron basis function exponents
* :mod:`basisset`: module for handling basis set in different formats and obtaining
  exponents from series like: even tempered, well tempered, Legendre expansion,
* :mod:`dmft`: utility functions for running DMFT calculations,
* :mod:`gamessus`: utility function for parsing input and log files from Gamess-US
  calculations, wrappers for running the code from python
* :mod:`molpro`: parser for the output file
* :mod:`molecule`: general purpose module intorducing molecule class for handling
  molecules


.. toctree::
   :maxdepth: 2

   basisopt
   basisset
   code
   dmft
   gamessus
   molecule
   molpro

Download and Installation
=========================

Prerequisites
-------------

* Python_ version 2.7.3 or later
* Periodic_ package
* Scipy_ verison 0.11 or later
* setuptools_ if you want to install via ``python setup.py install`` or
  ``easy_install``


.. _Python: http://python.org/
.. _Periodic: https://pypi.python.org/pypi/periodic
.. _SciPy: http://www.scipy.org/
.. _setuptools: https://pypi.python.org/pypi/setuptools




Installation
============

Using pip
---------

If you don't have ``pip`` you should probably me your life easier and get it,
the installation instructions can be found `here <https://pip.pypa.io/en/latest/installing.html>`_.
chemtools is **NOT** hosted on pypi yet but in can be installed by pip from the
bibbucket repository with::

    $ [sudo] pip install https://bitbucket.org/lukaszmentel/chemtools/get/tip.tar.gz


Using setup.py
--------------

The code is hosted via public bitbucket repository named chemtools_
and a source distibution can be downloaded from there as a zip, gz or bz2
bundle.

.. _chemtools: https://bitbucket.org/lukaszmentel/chemtools/

Or from CLI::

    wget https://bitbucket.org/lukaszmentel/chemtools/get/tip.tar.gz

First we have to unpack the tarball by::

    $ tar xzvf tip.tar.gz

go to the the main directory::

    $ cd lukaszmentel-chemtools-xxxxxxxxxxxx

and install the package by::

    $ [sudo] python setup.py install

Using easy_install
------------------

To install chemtools using ``easy_install`` just type::

    $ [sudo] easy_install https://bitbucket.org/lukaszmentel/chemtools/get/tip.tar.gz

Help
====

If you have some questions, remarks or requests email me at `<lmmentel@gmail.com> <mailto:lmmentel@gmail.com>`_.

Similar projects
================

`cclib <http://cclib.github.io/>`_, `pygamess <https://github.com/kzfm/pygamess>`_

License
=======

.. include:: ../../LICENSE.rst


Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
