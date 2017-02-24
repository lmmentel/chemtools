
import os
from chemtools.calculators.dalton import Dalton
from chemtools.molecule import Molecule
from chemtools.basisopt import BSOptimizer
import pytest


@pytest.mark.skipif(os.getenv('DALTON_EXE') is None,
                    reason="<DALTON_EXE> undefined")
def test_optimize(tmpdir):

    tmpdir.chdir()

    moltemplate = '''INTGRL
basis set optimization
test dalton calculator
Atomtypes=1 Integrals=1.0D-20
%basis
'''

    hftemplate = '''**DALTON INPUT
.RUN WAVE FUNCTIONS
**WAVE FUNCTION
.HF
**END OF DALTON INPUT
'''

    daltemplate = '''**DALTON INPUT
.RUN WAVE FUNCTIONS
**WAVE FUNCTIONS
.HF
.GASCI
*LUCITA
.INIWFC
HF_SCF
.CITYPE
GASCI
.SYMMET
1
.MULTIP
1
.INACTIVE
%core
.GAS SHELLS
1
0 2 / 39   30   30   20   30   20   20   14
.NROOTS
1
**END OF DALTON INPUT
'''

    dalton = Dalton(exevar='DALTON_EXE',
                    runopts=['-nobackup', '-noarch', '-d'])

    fname = {'dal': 'hf.dal', 'mol': 'He_5s.mol'}
    template = {'mol': moltemplate, 'dal': hftemplate}
    fsopt = {'He': [('s', 'et', 5, (0.5, 2.0))]}
    he = Molecule(name='He', atoms=[('He', )])
    bso = BSOptimizer(objective='hf total energy', template=template,
                      verbose=False, code=dalton, mol=he, fsopt=fsopt,
                      fname=fname)

    bso.run()
    energy = -2.8586246608170001
    assert abs(bso.result.fun - energy) < 1.0e-8, 'wrong objective'
