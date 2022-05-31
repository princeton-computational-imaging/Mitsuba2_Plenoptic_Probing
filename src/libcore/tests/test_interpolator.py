import mitsuba
import pytest
import enoki as ek
import numpy as np

def test01_1d(variant_scalar_rgb):
    from mitsuba.core import LinearInterpolator1

    f = lambda x: x**2
    res = 50
    a = 0.0
    b = 3.0

    xx = np.linspace(a, b, res)
    yy = f(xx)

    interp = LinearInterpolator1(yy, [xx])

    assert np.allclose(interp.eval(0.0), 0.0)
    for x in [0.5, 1.0, 1.5, 2.5, 2.75]:
        assert np.allclose(f(x), interp.eval(x), rtol=1e-2)

def test02_1d_linear(variant_scalar_rgb):
    from mitsuba.core import LinearInterpolator1

    xx = [2.0, 5.0]
    yy = [1.0, -2.0]

    interp = LinearInterpolator1(yy, [xx])

    assert np.allclose(interp.eval(xx[0]), yy[0])
    assert np.allclose(interp.eval(xx[1]), yy[1])
    x_mid = 0.5*(xx[0] + xx[1])
    y_mid = 0.5*(yy[0] + yy[1])
    assert np.allclose(interp.eval(x_mid), y_mid)

def test03_2d(variant_scalar_rgb):
    from mitsuba.core import LinearInterpolator2

    f = lambda x, y: np.exp(-0.5*(x**2 + y**2))

    res_x = 50
    res_y = 80
    xx = np.linspace(-1, 1, res_x)
    yy = np.linspace(-1.5, 1.5, res_y)
    grid_x, grid_y = np.meshgrid(xx, yy)
    zz = f(grid_x.flatten(), grid_y.flatten()).reshape(res_y, res_x)

    interp = LinearInterpolator2(zz, [yy, xx])

    assert np.allclose(f(0, 0), interp.eval(0, 0), rtol=1e-2)
    for x, y in [(0.5, -0.5), (1.0, 0.3), (0.3, -0.8), (-0.77, 0.2)]:
        assert np.allclose(f(x, y), interp.eval([y, x]), rtol=1e-2)

def test04_2d_linear(variant_scalar_rgb):
    from mitsuba.core import LinearInterpolator2

    xx = [1.0, 3.0]
    yy = [10.0, 20.0]
    zz = np.array([[-5.0, 10.0], [8.0, 1.0]])

    interp = LinearInterpolator2(zz, [yy, xx])

    assert np.allclose(interp.eval([yy[0], xx[0]]), zz[0, 0])
    assert np.allclose(interp.eval([yy[0], xx[1]]), zz[0, 1])
    assert np.allclose(interp.eval([yy[1], xx[0]]), zz[1, 0])
    assert np.allclose(interp.eval([yy[1], xx[1]]), zz[1, 1])

    y0_mid = 0.5*(zz[0, 0] + zz[0, 1])
    y1_mid = 0.5*(zz[1, 0] + zz[1, 1])
    mid = 0.5*(y0_mid + y1_mid)

    x_mid = 0.5*(xx[0] + xx[1])
    y_mid = 0.5*(yy[0] + yy[1])

    assert np.allclose(interp.eval([y_mid, x_mid]), mid)
