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

This emitter plugin implements an orthographic emitter.
*/

MTS_VARIANT class OrthographicEmitter final : public Emitter<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Emitter, m_flags, m_world_transform, m_needs_sample_3)
    MTS_IMPORT_TYPES(Scene, Texture)

    OrthographicEmitter(const Properties &props) : Base(props) {
        // these two specify resolution
        m_width = props.float_("width");
        m_height = props.float_("height");
        m_intensity_scale = props.float_("intensity_scale", 1.0f);
        // how far we're willing to let light travel altho it's kind of strange for this
        // context.
        m_near_clip = props.float_("near_clip");
        m_far_clip = props.float_("far_clip");

        m_delta_x = clamp(props.float_("delta_x", 1), 0, 1);
        m_delta_y = clamp(props.float_("delta_y", 1), 0, 1);

        if (m_world_transform->has_scale())
          Throw("Scale factors are not allowed in local-to-world transform");

        m_flags = EmitterFlags::DeltaPosition |
                  EmitterFlags::SpatiallyVarying;

        m_irradiance = props.texture<Texture>("irradiance", Texture::D65(1.f));
        m_needs_sample_3 = false;
        update_local_transforms();
    }

    void update_local_transforms() {
        m_sample_to_local =
            // scale [-.5, .5], " ", z to [-w/2, w/2], [-h/2, h/2], z
            ScalarTransform4f::scale(ScalarVector3f(-2.f * m_width, -2.f * m_height, 1.f)) *
            // move [0,1] to [-.5, .5]
            ScalarTransform4f::translate(ScalarVector3f(-0.5f, -0.5f, 0.f));

        m_local_to_sample = m_sample_to_local.inverse();
    }

    std::pair<DirectionSample3f, Spectrum> sample_direction(const Interaction3f &it,
                                                            // location on emitter plane
                                                            const Point2f &sample_,
                                                            Mask active) const override {
      MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);
      Point2f sample(sample_.x()/m_width, sample_.y()/m_height);
      auto trafo = m_world_transform->eval(it.time, active);

      DirectionSample3f ds;
      // starts from the world position of the perspective emitter
      ds.p = trafo.translation();
      // doesn't particularly matter but translate bkwd in z direction to world
      ds.n = trafo * Vector3f(0, 0, 1);

      ds.time = it.time;
      ds.object = this;

      // world direction
      ds.d = ds.p - it.p;
      ds.dist = norm(ds.p - it.p);
      Float inv_dist = rcp(ds.dist);
      ds.d *= inv_dist;

      auto si = zero<SurfaceInteraction3f>();
      // this works because it converts the position to local space
      // then uses that as a direction since the emitter is at the origin
      // and converts that to a sample
      auto local_it = trafo.inverse() * it.p;
      auto inter_sample = m_local_to_sample * local_it;
      // convert the direction into sample coordinates
      if (inter_sample.z() > m_far_clip || inter_sample.z() < m_near_clip)
        return std::make_pair(ds, Spectrum(0));

      si.uv = Point2f(inter_sample.x(), inter_sample.y());
      ds.uv = si.uv;
      // if out of range of image, then return no illumination
      if (any(si.uv < 0.f) || any(si.uv > 1.f)) return std::make_pair(ds, Spectrum(0));
      si.wavelengths = it.wavelengths;
      Spectrum s = m_irradiance->eval(si, active);


      ds.delta = true;
      ds.pdf = 1.0f;
      // apply damping factor to intensity
      return std::make_pair(ds, s * m_intensity_scale * inv_dist * inv_dist);
    }

    ScalarBoundingBox3f bbox() const override {
      return m_world_transform->translation_bounds();
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("irradiance", m_irradiance.get());
        // callback->put_parameter("fov", m_x_fov);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "OrthographicEmitter[" << std::endl
            << "  irradiance = " << string::indent(m_irradiance) << ","
            << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    ref<Texture> m_irradiance;
    ScalarTransform4f m_local_to_sample;
    ScalarTransform4f m_sample_to_local;
    ScalarFloat m_width;
    ScalarFloat m_height;
    ScalarFloat m_intensity_scale;

    ScalarFloat m_delta_x;
    ScalarFloat m_delta_y;
    // These don't in particular make sense for an emitter
    // But not sure what they should be replaced with
    ScalarFloat m_near_clip;
    ScalarFloat m_far_clip;
};

MTS_IMPLEMENT_CLASS_VARIANT(OrthographicEmitter, Emitter)
MTS_EXPORT_PLUGIN(OrthographicEmitter, "Orthographic emitter")
NAMESPACE_END(mitsuba)
