import mitsuba
import pytest
import enoki as ek


def create_emitter(o, d, fov=34, fov_axis='x', s_open=1.5, s_close=5):
    from mitsuba.core.xml import load_string
    return load_string("""<emitter version='2.0.0' type='perspective_light'>
                              <float name='near_clip' value='1'/>
                              <float name='far_clip' value='35'/>
                              <float name='fov' value='{fov}'/>
                              <string name='fov_axis' value='{fov_axis}'/>
                              <integer name="width" value="512"/>
                              <integer name="height" value="256"/>
                              <transform name="to_world">
                                  <lookat origin="{ox}, {oy}, {oz}"
                                          target="{tx}, {ty}, {tz}"
                                          up    =" 0.0,  1.0,  0.0"/>
                              </transform>
                          </emitter> """.format(ox=o[0], oy=o[1], oz=o[2],
                                               tx=o[0]+d[0], ty=o[1]+d[1], tz=o[2]+d[2],
                                               fov=fov, fov_axis=fov_axis, so=s_open, sc=s_close))


origins    = [[1.0, 0.0, 1.5], [1.0, 4.0, 1.5]]
directions = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]


@pytest.mark.parametrize("origin", origins)
@pytest.mark.parametrize("direction", directions)
@pytest.mark.parametrize("s_open", [0.0, 1.5])
@pytest.mark.parametrize("s_time", [0.0, 3.0])
def test01_create(variant_scalar_rgb, origin, direction, s_open, s_time):
    from mitsuba.core import BoundingBox3f, Vector3f, Transform4f

    emitter = create_emitter(origin, direction, s_open=s_open, s_close=s_open + s_time)

    assert ek.allclose(emitter.near_clip(), 1)
    assert ek.allclose(emitter.far_clip(), 35)
    assert ek.allclose(emitter.shutter_open(), s_open)
    assert camera.bbox() == BoundingBox3f(origin, origin)
    assert ek.allclose(camera.world_transform().eval(0).matrix,
                       Transform4f.look_at(origin, Vector3f(origin) + direction, [0, 1, 0]).matrix)


@pytest.mark.parametrize("origin", origins)
@pytest.mark.parametrize("direction", directions)
def test02_sample_ray(variant_packet_spectral, origin, direction):
    # Check the correctness of the sample_ray() method
    from mitsuba.core import sample_shifted, sample_rgb_spectrum

    emitter = create_emitter(origin, direction)

    time = 0.5
    wav_sample = [0.5, 0.33, 0.1]
    pos_sample = [[0.2, 0.1, 0.2], [0.6, 0.9, 0.2]]

    ray, spec_weight = camera.sample_ray(time, wav_sample, pos_sample, aperture_sample)

    # Importance sample wavelength and weight
    wav, spec = sample_rgb_spectrum(sample_shifted(wav_sample))

    assert ek.allclose(ray.wavelengths, wav)
    assert ek.allclose(spec_weight, spec)
    assert ek.allclose(ray.time, time)
    assert ek.allclose(ray.o, origin)

    # Check that a [0.5, 0.5] position_sample generates a ray
    # that points in the camera direction
    ray, _ = camera.sample_ray(0, 0, [0.5, 0.5], 0)
    assert ek.allclose(ray.d, direction, atol=1e-7)



@pytest.mark.parametrize("origin", origins)
@pytest.mark.parametrize("direction", directions)
def test03_sample_ray_differential(variant_packet_spectral, origin, direction):
    # Check the correctness of the sample_ray_differential() method
    from mitsuba.core import sample_shifted, sample_rgb_spectrum

    camera = create_camera(origin, direction)

    time = 0.5
    wav_sample = [0.5, 0.33, 0.1]
    pos_sample = [[0.2, 0.1, 0.2], [0.6, 0.9, 0.2]]

    ray, spec_weight = camera.sample_ray_differential(time, wav_sample, pos_sample, 0)

    # Importance sample wavelength and weight
    wav, spec = sample_rgb_spectrum(sample_shifted(wav_sample))

    assert ek.allclose(ray.wavelengths, wav)
    assert ek.allclose(spec_weight, spec)
    assert ek.allclose(ray.time, time)
    assert ek.allclose(ray.o, origin)

    # Check that the derivatives are orthogonal
    assert ek.allclose(ek.dot(ray.d_x - ray.d, ray.d_y - ray.d), 0, atol=1e-7)

    # Check that a [0.5, 0.5] position_sample generates a ray
    # that points in the camera direction
    ray_center, _ = camera.sample_ray_differential(0, 0, [0.5, 0.5], 0)
    assert ek.allclose(ray_center.d, direction, atol=1e-7)

    # Check correctness of the ray derivatives

    # Deltas in screen space
    dx = 1.0 / camera.film().crop_size().x
    dy = 1.0 / camera.film().crop_size().y

    # Sample the rays by offsetting the position_sample with the deltas
    ray_dx, _ = camera.sample_ray_differential(0, 0, [0.5 + dx, 0.5], 0)
    ray_dy, _ = camera.sample_ray_differential(0, 0, [0.5, 0.5 + dy], 0)

    assert ek.allclose(ray_dx.d, ray_center.d_x)
    assert ek.allclose(ray_dy.d, ray_center.d_y)


@pytest.mark.parametrize("origin", [[1.0, 0.0, 1.5]])
@pytest.mark.parametrize("direction", [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
@pytest.mark.parametrize("fov", [34, 80])
def test04_fov_axis(variant_packet_spectral, origin, direction, fov):
    # Check that sampling position_sample at the extrimities of the unit square
    # along the fov_axis should generate a ray direction that make angle of fov/2
    # with the camera direction.

    from mitsuba.core import sample_shifted, sample_rgb_spectrum

    def check_fov(camera, sample):
        ray, _ = camera.sample_ray(0, 0, sample, 0)
        assert ek.allclose(ek.acos(ek.dot(ray.d, direction)) * 180 / ek.pi, fov / 2)

    # In the configuration, aspect==1.5, so 'larger' should give the 'x'-axis
    for fov_axis in ['x', 'larger']:
        camera = create_camera(origin, direction, fov=fov, fov_axis=fov_axis)
        for sample in [[0.0, 0.5], [1.0, 0.5]]:
            check_fov(camera, sample)

    # In the configuration, aspect==1.5, so 'smaller' should give the 'y'-axis
    for fov_axis in ['y', 'smaller']:
        camera = create_camera(origin, direction, fov=fov, fov_axis=fov_axis)
        for sample in [[0.5, 0.0], [0.5, 1.0]]:
                check_fov(camera, sample)

    # Check the 4 corners for the `diagonal` case
    camera = create_camera(origin, direction, fov=fov, fov_axis='diagonal')
    for sample in [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]]:
            check_fov(camera, sample)


