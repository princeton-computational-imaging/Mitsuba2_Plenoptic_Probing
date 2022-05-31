#include <iostream>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _emitter-perspective:

Distant directional emitter (:monosp:`perspective`)
Now perspective emitter
---------------------------------------------------

.. pluginparameters::

    * - irradiance
      - |spectrum|
      - Spectral irradiance, which corresponds to the amount of spectral power
        per unit area received by a hypothetical surface normal to the specified
        direction.
    * - to_world
      - |transform|
      - Emitter-to-world transformation matrix.
    * - direction
      - |vector|
      - Alternative (and exclusive) to `to_world`. Direction towards which the
        emitter is radiating in world coordinates.
    * - fov
      - |float|
      - denotes the emitters range of emission between 0 and 180, excluding
      - the extremes

This emitter plugin implements a perspective source

*/

MTS_VARIANT class PerspectiveEmitter final : public Emitter<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Emitter, m_flags, m_world_transform, m_needs_sample_3)
    MTS_IMPORT_TYPES(Scene, Texture)

    PerspectiveEmitter(const Properties &props) : Base(props) {
        m_width = props.float_("width");
        m_height = props.float_("height");
        // store so don't need to recompute aspect ratio
        m_aspect = m_width/m_height;
        m_x_fov = parse_fov(props, m_aspect);

        m_intensity_scale = props.float_("intensity_scale", 1.0f);
        // how far we're willing to let light travel.
        m_far_clip = props.float_("far_clip");
        m_near_clip = props.float_("near_clip");

        if (m_world_transform->has_scale())
          Throw("Scale factors are not allowed in local-to-world transform");

        m_flags = EmitterFlags::DeltaPosition |
                  EmitterFlags::SpatiallyVarying;

        m_irradiance = props.texture<Texture>("irradiance", Texture::D65(1.f));
        m_needs_sample_3 = false;
        update_local_transforms();
    }

    void update_local_transforms() {
        /**
         * These do the following (in reverse order):
         *
         * 1. Create transform from camera space to [-1,1]x[-1,1]x[0,1] clip
         *    coordinates (not taking account of the aspect ratio yet)
         *
         * 2+3. Translate and scale to shift the clip coordinates into the
         *    range from zero to one, and take the aspect ratio into account.
         *
         * Kept these in but replaced with width and height because those are the clipping
         * bounds locally. Not sure if this is correct. *Edit removed.
         * 4+5. Translate and scale the coordinates once more to account
         *     for a cropping window (if there is any)
         */
        m_local_to_sample =
            // [0, 1], [0, 1], [0, 1]
            ScalarTransform4f::scale(ScalarVector3f(-0.5f, -0.5f * m_aspect, 1.f)) *
            // [-2, 0], [-2, 0], [0, 1], shaves some stuff off with the aspect
            // ratio.
            ScalarTransform4f::translate(ScalarVector3f(-1.f, -1.f/m_aspect, 0.f)) *
            // [-1, 1], [-1, 1], [0, 1]
            ScalarTransform4f::perspective(m_x_fov, m_near_clip, m_far_clip);
            // camera space (camera at origin in positive z-direction)
        m_sample_to_local = m_local_to_sample.inverse();

        /* Precompute some data for importance(). Please
           look at that function for further details.

           This comment was in sensors/perspective.cpp but that file doesn't
           even have an importance function, it seems to be leftover from mitsuba 1. */
        ScalarPoint3f pmin(m_sample_to_local * ScalarPoint3f(0.f, 0.f, 0.f)),
                      pmax(m_sample_to_local * ScalarPoint3f(1.f, 1.f, 0.f));

        m_image_rect.reset();
        m_image_rect.expand(ScalarPoint2f(pmin.x(), pmin.y()) / pmin.z());
        m_image_rect.expand(ScalarPoint2f(pmax.x(), pmax.y()) / pmax.z());
        m_normalization = rcp(m_image_rect.volume());
    }

    std::pair<DirectionSample3f, Spectrum> sample_direction(const Interaction3f &it,
                                                            const Point2f & /* sample_ */,
                                                            Mask active) const override {
      MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);
      auto trafo = m_world_transform->eval(it.time, active);

      DirectionSample3f ds;
      // starts from the world position of the perspective emitter
      ds.p = trafo.translation();
      // doesn't particularly matter but translate fwd in z direction to world
      ds.n = trafo * Vector3f(0, 0, 1);

      // Whether this is drawn from a dirac delta or not, yes since we're only drawing from this
      // one point
      ds.delta = true;
      ds.pdf = 1.f;
      ds.time = it.time;
      ds.object = this;


      // this works because it converts the position to local space
      // then uses that as a direction since the emitter is at the origin
      // and converts that to a sample
      auto local_it = trafo.inverse() * it.p;
      auto inter_sample = m_local_to_sample * local_it;
      // convert the direction into sample coordinates
      if (any(inter_sample > 1.f) || any(inter_sample < 0.f))
        return std::make_pair(ds, Spectrum(0));
      ds.uv = Point2f(inter_sample.x(), inter_sample.y());

      auto si = zero<SurfaceInteraction3f>();
      // sample from the position on the emitter, not the intersection pt
      si.uv = ds.uv;
      si.wavelengths = it.wavelengths;
      Spectrum s = m_irradiance->eval(si, active);

      // world direction
      ds.d = ds.p - it.p;
      ds.dist = norm(ds.d);
      Float inv_dist = rcp(ds.dist);
      ds.d *= inv_dist;

      // apply damping factor to intensity
      return std::make_pair(ds, unpolarized<Spectrum>(s) * m_intensity_scale * inv_dist * inv_dist);
    }

    ScalarBoundingBox3f bbox() const override {
      return m_world_transform->translation_bounds();
    }

    void traverse(TraversalCallback *callback) override {
        Throw("test");
        callback->put_object("irradiance", m_irradiance.get());
        // callback->put_parameter("fov", m_x_fov);
    }
    Spectrum eval(const SurfaceInteraction3f &, Mask) const override {
      Throw("eval");
      return 1.f;
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "PerspectiveEmitter[" << std::endl
            << "  irradiance = " << string::indent(m_irradiance) << ","
            << "  fov = " << string::indent(m_x_fov) << ","
            << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    ref<Texture> m_irradiance;
    ScalarTransform4f m_local_to_sample;
    ScalarTransform4f m_sample_to_local;
    ScalarBoundingBox2f m_image_rect;
    ScalarFloat m_x_fov;
    ScalarFloat m_width;
    ScalarFloat m_intensity_scale;
    ScalarFloat m_height;
    ScalarFloat m_aspect;

    // TODO figure out what this was used for prior?
    ScalarFloat m_normalization;
    // These don't in particular make sense for an emitter
    // But not sure what they should be replaced with
    ScalarFloat m_near_clip;
    ScalarFloat m_far_clip;
};

MTS_IMPLEMENT_CLASS_VARIANT(PerspectiveEmitter, Emitter)
MTS_EXPORT_PLUGIN(PerspectiveEmitter, "Perspective emitter")
NAMESPACE_END(mitsuba)
