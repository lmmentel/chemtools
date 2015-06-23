========================
basisset module tutorial
========================

The main purpose of the :mod:`basisset` module is to provide a common format
for manipulating basis sets in terms of the `BasisSet` class instances. The
`BasisSet` class can be instantiated from.

You can browse an ipython notebook with a short tutorial `here
<http://nbviewer.ipython.org/urls/bitbucket.org/lukaszmentel/chemtools/raw/tip/examples/ipython_notebooks/BasisSetTutorial.ipynb>`_

Initialization
==============

The `BasisSet` class can be instantiated through three basic mechanisms:
- intialize a basis set from a sequence
- parse a basis set from a file
- by directly specyfinig the `functions`

Initialize from a sequence
------------------------

There are three kinds of sequences that can be used to generate exponents:
- even tempered expansion,
- well tempered expansion,
- legendre polynomial exapansion.

First import the `BasisSet` object from the `basisset` module

.. code-block:: python

   >>> from chemtools.basisset import BasisSet

Class method `from_sequence` can be used to  generate a `BasisSet` object. The
method takes

.. code-block:: python

   >>> bas = BasisSet.from_sequence(formula='even tempered',
   ...                              name='eventemp basis',
   ...                              element='H',
   ...                              functs=[('s', 10, (0.5, 2.0)])
   >>> bas.name
   'eventemp basis'
   >>> bas.element
   'H'
   >>> bas.functions
