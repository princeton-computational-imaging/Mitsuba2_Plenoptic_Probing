#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/phase.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _phase-isotropic:

Isotropic phase function (:monosp:`isotropic`)
-----------------------------------------------

This phase function simulates completely uniform scattering,
where all directionality is lost after a single scattering
interaction. It does not have any parameters.

*/

// Derived from https://github.com/scottprahl/miepython/blob/master/miepython/miepython.py

template <typename Float, typename Spectrum>
class DepolarizePhaseFunction final : public PhaseFunction<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(PhaseFunction, m_flags)
    MTS_IMPORT_TYPES(PhaseFunctionContext)

    DepolarizePhaseFunction(const Properties & props) : Base(props) {
      m_flags = +PhaseFunctionFlags::Isotropic;
    }

    std::pair<Vector3f, Float> sample(const PhaseFunctionContext & /* ctx */,
                                      const MediumInteraction3f & /* mi */, const Point2f &sample,
                                      Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionSample, active);

        auto wo  = warp::square_to_uniform_sphere(sample);
        auto pdf = warp::square_to_uniform_sphere_pdf(wo);
        return std::make_pair(wo, pdf);
    }

    Float eval(const PhaseFunctionContext & /* ctx */, const MediumInteraction3f & /* mi */,
               const Vector3f &wo, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::PhaseFunctionEvaluate, active);
        return warp::square_to_uniform_sphere_pdf(wo);
    }

    Spectrum polarization(const PhaseFunctionContext &/*ctx*/, const Spectrum incoming,
               const MediumInteraction3f &/*mi*/, const Vector3f/*wo*/, Mask /*active*/) const override {
      return depolarize(incoming);
    }

    std::string to_string() const override { return "DepolarizePhaseFunction[]"; }


    MTS_DECLARE_CLASS()
private:

};

MTS_IMPLEMENT_CLASS_VARIANT(DepolarizePhaseFunction, PhaseFunction)
MTS_EXPORT_PLUGIN(DepolarizePhaseFunction, "Depolarizing phase function")
NAMESPACE_END(mitsuba)
