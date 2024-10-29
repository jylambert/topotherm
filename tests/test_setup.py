"""test_setup.py tests the correct setup of topotherm."""

from pytest import approx


def test_import():
    import topotherm


def test_functionality():
    from topotherm import single_timestep
    assert single_timestep.annuity(c_i=0.01, n=10) == approx(0.1055820766, rel=1e-2)
