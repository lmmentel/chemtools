.. image:: https://readthedocs.org/projects/chemtools/badge/
   :target: https://chemtools.readthedocs.org
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/chemtools.svg?style=flat-square&label=PyPI%20version
   :target: https://pypi.python.org/pypi/chemtools
   :alt: Latest version released on PyPi

.. image:: https://www.travis-ci.org/lmmentel/chemtools.svg?branch=master
    :target: https://www.travis-ci.org/lmmentel/chemtools
    :alt: Build Status

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
    :target: https://opensource.org/licenses/MIT
    :alt: MIT license

.. image:: https://pepy.tech/badge/chemtools
    :target: https://pepy.tech/project/chemtools
    :alt: pepy


======================================================
chemtools_: Python toolbox for Computational Chemistry
======================================================

Chemtools is a set of modules that is intended to help with more
advanced computations using common electronic structure programs.
Currently there is support for:

* Dalton_
* Gamess-US_
* MolPro_
* PSI4_


Table of Contents
=================

* `Getting Started`_
  
  * Installation_
  * Documentation_

* Contributing_
* Contact_
* Citing_
* Funding_
* License_

Getting Started
===============

Installation
------------

Most convenient way to install the package is with `pip <https://pip.pypa.io/en/stable/>`_  

.. code-block:: bash

   pip install chemtools


Documentation
-------------

The documentation in hosted at `Read The Docs <http://chemtools.readthedocs.org/en/latest/>`_.


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

License
=======

The project is distributed under the MIT License. See `LICENSE <LICENSE.rst>`_ for more information.

.. _chemtools: http://chemtools.readthedocs.org
.. _Gamess-US: https://www.msg.chem.iastate.edu/gamess/gamess.html
.. _MolPro: http://www.molpro.net/
.. _Dalton: https://www.daltonprogram.org/
.. _PSI4: http://www.psicode.org/
