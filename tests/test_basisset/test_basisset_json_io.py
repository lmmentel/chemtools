
import numpy as np
from chemtools.basisset import BasisSet


VTZGAUSSIAN = """
****
Be     0
S  11   1.00
     6863.0000000000     0.00023600
     1030.0000000000     0.00182600
      234.7000000000     0.00945200
       66.5600000000     0.03795700
       21.6900000000     0.11996500
        7.7340000000     0.28216200
        2.9160000000     0.42740400
        1.1300000000     0.26627800
        0.2577000000     0.01819300
        0.1101000000    -0.00727500
        0.0440900000     0.00190300
S  11   1.00
     6863.0000000000    -0.00004300
     1030.0000000000    -0.00033300
      234.7000000000    -0.00173600
       66.5600000000    -0.00701200
       21.6900000000    -0.02312600
        7.7340000000    -0.05813800
        2.9160000000    -0.11455600
        1.1300000000    -0.13590800
        0.2577000000     0.22802600
        0.1101000000     0.57744100
        0.0440900000     0.31787300
S   1   1.00
        0.2577000000     1.00000000
S   1   1.00
        0.0440900000     1.00000000
P   5   1.00
        7.4360000000     0.01073600
        1.5770000000     0.06285400
        0.4352000000     0.24818000
        0.1438000000     0.52369900
        0.0499400000     0.35342500
P   1   1.00
        0.1438000000     1.00000000
P   1   1.00
        0.0499400000     1.00000000
D   1   1.00
        0.3493000000     1.00000000
D   1   1.00
        0.1724000000     1.00000000
F   1   1.00
        0.3423000000     1.00000000
****
"""


def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


def test_to_json():

    bs = BasisSet.from_str(VTZGAUSSIAN, fmt='gaussian', name='VTZ')

    bsdumped = bs.to_json(indent=4)
    bsloaded = BasisSet.from_json(bsdumped)

    assert bs.name == bsloaded.name, 'inconsistent name'
    assert bs.element == bsloaded.element, 'inconsistent element'

    for shell, funs in bs.functions.items():

        assert shell in bsloaded.functions.keys(), 'missing shell {}'.format(shell)

        assert np.allclose(funs['e'], bsloaded.functions[shell]['e'])

        for f1, f2 in zip(funs['cf'], bsloaded.functions[shell]['cf']):
            assert np.allclose(f1['idx'], f2['idx']), 'inconsistent idx'
            assert np.allclose(f1['cc'], f2['cc']), 'inconsistent cc'
