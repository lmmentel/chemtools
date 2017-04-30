
import pytest
from chemtools.basisopt import BSOptimizer


def test_fsopt_validation():

    with pytest.raises(ValueError):

        fsopt = {'Be': [('z', 'exp', 1, (0.27173,))]}
        BSOptimizer(fsopt=fsopt)

    with pytest.raises(ValueError):

        fsopt = {'Be': [('s', 'other', 1, (0.27173,))]}
        BSOptimizer(fsopt=fsopt)

    with pytest.raises(ValueError):

        fsopt = {'Be': [('s', None, 1, (0.27173,))]}
        BSOptimizer(fsopt=fsopt)
