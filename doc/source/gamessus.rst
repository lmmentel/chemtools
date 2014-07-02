gamessus module
===============

``chemtools.gamessus`` module contains four classes for handling various tasks while
working with **GAMESS(US)**:

.. automodule:: gamessus
   :members:



*chemtools.gamessus.GamessReader*
    holding methods for reading gamess binary files:

    * dictionary file (\*.F10), contensts of this file are descibed in the **GAMESS(US)** documentations found `here <http://www.msg.ameslab.gov/gamess/GAMESS_Manual/prog.pdf>`_,
    * two-electron integrals in AO basis (\*.F08),
    * two-electron integrals in MO basis (\*.F09),
    * two-particle reduced density matrix 2RDM from CI calcualtions, (\*.F15).

*chemtools.gamessus.GamessDatParser*
    holding methods for parsig the PUNCH file (\*.F07).
