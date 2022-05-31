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
class MiePhaseFunction final : public PhaseFunction<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(PhaseFunction, m_flags)
    MTS_IMPORT_TYPES(PhaseFunctionContext)

    MiePhaseFunction(const Properties & props) : Base(props) {
      m_flags = +PhaseFunctionFlags::Isotropic;
      m_size = props.float_("size", 0.001);
      Assert(m_size <= 0.1);
      m_ior = props.float_("ior");
      init();
    }

    void init() {
      // precompute constants which will be used when calculating the mueller matrix
      auto m_2 = m_ior * m_ior;
      auto m_4 = m_2 * m_2;

      auto x_2 = m_size * m_size;
      x_3 = x_2 * m_size;
      auto x_4 = x_2 * x_2;
      auto ricatti_bessel = m_2 + 2.f + (1.f - 0.7f * m_2) * x_2 - (8.f * m_4 - 385.f * m_2
        + 350.f) * x_4 / 1400.f + Complex(0.f, 2.f) * (m_2 - 1.f)
          * x_3 * (1.f - 0.1f * x_2) / 3.f;

      a_hat_1 = Complex(0.f, 2.f) * (m_2 - 1.f) / 3.f * (1.f - 0.1f * x_2 + (4.f * m_2 + 5.f) *
      x_4 / 1400.f) / ricatti_bessel;

      b_hat_1 = Complex(0.f, 1.f) * x_2 * (m_2 - 1.f) / 45.f * (1.f + (2.f * m_2 - 5.f) /
                                       70.f * x_2) / (1.f - (2.f * m_2 - 5.f) / 30.f * x_2);

      a_hat_2 = Complex(0.f, 1.f) * x_2 * (m_2 - 1.f) / 15.f * (1.f - x_2 / 14.f) /
              (2.f * m_2 + 3.f - (2.f * m_2 - 7.f) / 14.f * x_2);

      norm = sqrt(Float(M_PI) * 6.f * x_3 * real(a_hat_1 + b_hat_1 + 5.f * a_hat_2 / 3.f));
    }

    std::pair<Complex<ScalarFloat>, Complex<ScalarFloat>> s1_s2(Float cos_theta) const {
      Float mu(cos_theta);
      Complex<ScalarFloat> s_1 = 1.5f * x_3 *
        (a_hat_1 + b_hat_1 * mu + 5.f / 3.f * a_hat_2 * mu);

      Complex<ScalarFloat> s_2 = 1.5f * x_3 *
        (b_hat_1 + a_hat_1 * mu + 5.f / 3.f * a_hat_2 * (2.f * mu*mu - 1.f));
      return { s_1/norm, s_2/norm };
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

    Spectrum polarization(const PhaseFunctionContext &ctx, const Spectrum incoming,
               const MediumInteraction3f & mi, const Vector3f wo, Mask /*active*/) const override {
      if constexpr (is_polarized_v<Spectrum>) {

        auto w_i = normalize(mi.sh_frame.to_world(mi.wi));
        auto w_o = normalize(wo);
        // Cosine theta of polar scattering angle
        auto cos_theta = dot(w_i, w_o);
        auto [s1, s2] = s1_s2(cos_theta);
        MuellerMatrix<Float> M = mueller::mie(s1, s2);

        if (ctx.mode != TransportMode::Radiance) {
          // Flip if doing importance sampling? This is probably never reached
          w_o = -w_o;
          w_i = -w_i;
        }
        if (w_o == w_i) {
          // same scattering plane with no rotation so can just apply.
          return M * incoming;
        } else {
          // project outgoing direction onto incoming to find normal.
          // This is to manually build a basis that is coplanar with the outgoing direction
          auto y_i = w_o - (w_i * dot(w_i, wo));
          auto x_i = normalize(cross(y_i, w_i));
          y_i = cross(x_i, w_i);
          M = mueller::rotate_mueller_basis_collinear(M, w_i,
            Vector3f(1.f, 0.f, 0.f), y_i);
          return M * incoming;
        }
      } else {
        // if non-polarized just do nothing
        return incoming;
      }
    }

    std::string to_string() const override { return "MiePhaseFunction[]"; }


    MTS_DECLARE_CLASS()
private:
  Complex<ScalarFloat> a_hat_1;
  Complex<ScalarFloat> a_hat_2;
  Complex<ScalarFloat> b_hat_1;
  ScalarFloat norm;
  ScalarFloat x_3;

  ScalarFloat m_size;
  ScalarFloat m_ior;
};

MTS_IMPLEMENT_CLASS_VARIANT(MiePhaseFunction, PhaseFunction)
MTS_EXPORT_PLUGIN(MiePhaseFunction, "Mie phase function")
NAMESPACE_END(mitsuba)
