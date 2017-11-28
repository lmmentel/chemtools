============
Installation
============

Depenedencies
=============

- Python_ works with Python 2.7.x  and 3.x
- numba_
- numpy_
- mendeleev_
- scipy_
- setuptools_

Additional packages required to build the documentation:

- guzzle_sphinx_theme_
- nbsphinx_
- sphinx_


Using pip
=========

If you don't have ``pip`` you should probably me your life easier and get it,
the installation instructions can be found `here <https://pip.pypa.io/en/latest/installing.html>`_.
chemtools is **NOT** hosted on pypi yet but in can be installed by pip from the
bibbucket repository with::

    $ [sudo] pip install [-U] chemtools 


Using setup.py
==============

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

.. _guzzle_sphinx_theme: https://github.com/guzzle/guzzle_sphinx_theme
.. _mendeleev: http://mendeleev.readthedocs.io/en/stable/ 
.. _nbsphinx: http://nbsphinx.readthedocs.io/en/latest/
.. _numba: http://numba.pydata.org/
.. _numpy: http://www.numpy.org
.. _Python: http://python.org/
.. _scipy: http://www.scipy.org
.. _setuptools: https://pypi.python.org/pypi/setuptools
.. _sphinx: http://www.sphinx-doc.org/en/stable/
