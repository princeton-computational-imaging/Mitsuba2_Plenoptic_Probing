#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/tensor.h>
#include <mitsuba/core/interpolator.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/microfacet.h>

NAMESPACE_BEGIN(mitsuba)

// Normally, we sample with a mix of GGX and cosine hemisphere distributions
// With this flag we fall back to purely uniform sampling
// #define COSINE_SAMPLING

#define DIFFUSE_PDF 0.1f

template <typename Float, typename Spectrum>
class MeasuredPolarized final : public BSDF<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(BSDF, m_flags, m_components)
    MTS_IMPORT_TYPES(Texture, MicrofacetDistribution)

    using Interpolator6 = LinearInterpolator<Float, 6>;

    MeasuredPolarized(const Properties &props) : Base(props) {
        if constexpr (!is_polarized_v<Spectrum> || !is_spectral_v<Spectrum>)
            Throw("The measured polarized BSDF model requires that rendering takes place in spectral+polarized mode!");

        m_flags = BSDFFlags::GlossyReflection | BSDFFlags::FrontSide;
        m_components.push_back(m_flags);

        auto fs = Thread::thread()->file_resolver();
        fs::path file_path = fs->resolve(props.string("filename"));
        m_name = file_path.filename().string();

        // Optionally overwrite this to lookup only a specific wavelength
        m_lambda = props.float_("lambda", -1.f);

        ref<TensorFile> tf = new TensorFile(file_path);

        auto theta_h     = tf->field("theta_h");
        auto theta_d     = tf->field("theta_d");
        auto phi_d       = tf->field("phi_d");
        auto wavelengths = tf->field("wavelengths");
        auto pbrdf       = tf->field("pbrdf_f");

        if (!(theta_h.shape.size() == 1 &&
              theta_h.dtype == Struct::Type::Float32 &&

              theta_d.shape.size() == 1 &&
              theta_d.dtype == Struct::Type::Float32 &&

              phi_d.shape.size() == 1 &&
              phi_d.dtype == Struct::Type::Float32 &&

              wavelengths.shape.size() == 1 &&
              wavelengths.dtype == Struct::Type::Float32 &&

              pbrdf.dtype == Struct::Type::Float32 &&
              pbrdf.shape.size() == 6 &&
              pbrdf.shape[0] == 4 &&
              pbrdf.shape[1] == 4 &&
              pbrdf.shape[2] == phi_d.shape[0] &&
              pbrdf.shape[3] == theta_d.shape[0] &&
              pbrdf.shape[4] == theta_h.shape[0] &&
              pbrdf.shape[5] == wavelengths.shape[0])) {
            Throw("Invalid file structure: %s", tf->to_string());
        }

        ScalarFloat mm[4] = { 0.f, 1.f, 2.f, 3.f };

        m_interpolator = Interpolator6(
            (ScalarFloat *) pbrdf.data,
            {{ 4, 4,
                (uint32_t) phi_d.shape[0],
                (uint32_t) theta_d.shape[0],
                (uint32_t) theta_h.shape[0],
                (uint32_t) wavelengths.shape[0] }},
            {{ mm, mm,
               (const ScalarFloat *) phi_d.data,
               (const ScalarFloat *) theta_d.data,
               (const ScalarFloat *) theta_h.data,
               (const ScalarFloat *) wavelengths.data }}
        );

        m_alpha_sample = props.float_("alpha_sample", 0.5f);    // TODO: Make part of tensor
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float sample1,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs;
        if (unlikely(none_or<false>(active) || !ctx.is_enabled(BSDFFlags::GlossyReflection)))
            return { bs, 0.f };

#ifdef COSINE_SAMPLING
        bs.wo = warp::square_to_cosine_hemisphere(sample2);
        bs.pdf = warp::square_to_cosine_hemisphere_pdf(bs.wo);
#else
        MicrofacetDistribution distr(MicrofacetType::GGX,
                                     m_alpha_sample, m_alpha_sample, true);

        Float lobe_pdf_diffuse = DIFFUSE_PDF;
        Mask sample_diffuse    = active && sample1 < lobe_pdf_diffuse,
             sample_microfacet = active && !sample_diffuse;

        Vector3f wo_diffuse    = warp::square_to_cosine_hemisphere(sample2);
        auto [m, unused] = distr.sample(si.wi, sample2);
        Vector3f wo_microfacet = reflect(si.wi, m);

        bs.wo[sample_diffuse]    = wo_diffuse;
        bs.wo[sample_microfacet] = wo_microfacet;

        bs.pdf = pdf(ctx, si, bs.wo, active);
#endif

        bs.sampled_component = 0;
        bs.sampled_type = +BSDFFlags::GlossyReflection;
        bs.eta = 1.f;

        Spectrum value = eval(ctx, si, bs.wo, active);
        return { bs, select(active && bs.pdf > 0, value / bs.pdf, 0.f) };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);
        active &= (cos_theta_i > 0.f && cos_theta_o > 0.f);

        if (unlikely(none_or<false>(active) || !ctx.is_enabled(BSDFFlags::GlossyReflection)))
            return 0.f;

        /* Due to lack of reciprocity in polarization-aware pBRDFs, they are
           always evaluated w.r.t. the actual light propagation direction, no
           matter the transport mode. In the following, 'wi_hat' is toward the
           light source. */
        Vector3f wi_hat = ctx.mode == TransportMode::Radiance ? wo : si.wi,
                 wo_hat = ctx.mode == TransportMode::Radiance ? si.wi : wo;

        /* We now transform both directions to the standard frame defined in
           Figure (??). Here, one of the directions is aligned with the x-axis. */
        Float phi_std = phi(wo_hat);
        Vector3f wi_std = rotate_vector(wi_hat, Vector3f(0,0,1), -phi_std),
                 wo_std = rotate_vector(wo_hat, Vector3f(0,0,1), -phi_std);

        /* This representation can be turned into the (isotropic) Rusinkiewicz
           parameterization. */
        auto [phi_d, theta_h, theta_d] = directions_to_rusinkiewicz(wi_std, wo_std);

        Spectrum value;
        if constexpr (is_spectral_v<Spectrum>) {
            if constexpr (is_polarized_v<Spectrum>) {
                /* The Stokes reference frame vector of this matrix lies in the plane
                   of reflection. See Figure (??) */
                Vector3f zi_std = -wi_std,
                         ti_std = normalize(cross(wi_std - wo_std, zi_std)),
                         yi_std = normalize(cross(ti_std, zi_std)),
                         xi_std = cross(yi_std, zi_std),
                         zo_std = wo_std,
                         to_std = normalize(cross(wo_std - wi_std, zo_std)),
                         yo_std = normalize(cross(to_std, zo_std)),
                         xo_std = cross(yo_std, zo_std);

                // Reverse phi rotation from above on Stokes reference frames
                Vector3f xi_hat = rotate_vector(xi_std, Vector3f(0,0,1), +phi_std),
                         xo_hat = rotate_vector(xo_std, Vector3f(0,0,1), +phi_std);

                if (m_lambda == -1.f) {
                    for (int i = 0; i < 4; ++i) {
                        for (int j = 0; j < 4; ++j) {
                            UnpolarizedSpectrum tmp(0.f);
                            for (size_t k = 0; k < array_size_v<UnpolarizedSpectrum>; ++k) {
                                Float params[6] = {
                                    Float(j), Float(i),
                                    phi_d, theta_d, theta_h,
                                    si.wavelengths[k]
                                };
                                tmp[k] = m_interpolator.eval(params, active);
                            }
                            value(i, j) = tmp;
                        }
                    }
                } else {
                    for (int i = 0; i < 4; ++i) {
                        for (int j = 0; j < 4; ++j) {
                            Float params[6] = {
                                Float(j), Float(i),
                                phi_d, theta_d, theta_h,
                                Float(m_lambda)
                            };
                            value(i,j) = m_interpolator.eval(params, active);
                        }
                    }
                }

                if (any(isnan(value(0,0)))) {
                    value = 0.f;
                }

                // Make sure intensity is non-negative
                value(0,0) = max(0.f, value(0,0));

                /* Rotate in/out reference vector of value s.t. it aligns with the
                   implicit Stokes bases of -wi_hat & wo_hat. */
                value = mueller::rotate_mueller_basis(value,
                                                  -wi_hat, xi_hat, mueller::stokes_basis(-wi_hat),
                                                   wo_hat, xo_hat, mueller::stokes_basis(wo_hat));
            } else {
                if (m_lambda == -1.f) {
                    for (size_t k = 0; k < array_size_v<UnpolarizedSpectrum>; ++k) {
                        Float params[6] = {
                            0.f, 0.f,
                            phi_d, theta_d, theta_h,
                            si.wavelengths[k]
                        };
                        value[k] = m_interpolator.eval(params, active);
                    }
                } else {
                    Float params[6] = {
                        0.f, 0.f,
                        phi_d, theta_d, theta_h,
                        Float(m_lambda)
                    };
                    Float value_ = m_interpolator.eval(params, active);
                    value = Spectrum(value_);
                }

                // Make sure BRDF is non-negative
                value = max(0.f, value);
            }
        }

        return (value * cos_theta_o) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (unlikely(none_or<false>(active) || !ctx.is_enabled(BSDFFlags::GlossyReflection)))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

#ifdef COSINE_SAMPLING
        Float pdf = warp::square_to_cosine_hemisphere_pdf(wo);
#else
        MicrofacetDistribution distr(MicrofacetType::GGX,
                                     m_alpha_sample, m_alpha_sample, true);

        Vector3f H = normalize(wo + si.wi);

        Float pdf_diffuse = warp::square_to_cosine_hemisphere_pdf(wo);
        Float pdf_microfacet = distr.pdf(si.wi, H) / (4.f * dot(wo, H));

        Float pdf = 0.f;
        pdf += pdf_diffuse * DIFFUSE_PDF;
        pdf += pdf_microfacet * (1.f - DIFFUSE_PDF);

#endif
        return select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "MeasuredPolarized[" << std::endl
            << "  name = " << m_name << std::endl
            << "]";
        return oss.str();
    }

private:
    template <typename Vector3,
              typename Value = value_t<Vector3>>
    Value phi(const Vector3 &v) const {
        Value p = atan2(v.y(), v.x());
        masked(p, p < 0) += 2.f*math::Pi<Float>;
        return p;
    }

    template <typename Vector3, typename Value = value_t<Vector3>>
    MTS_INLINE
    Vector3 rotate_vector(const Vector3 &v, const Vector3 &axis_, Value angle) const {
        Vector3 axis = normalize(axis_);

        Value cos_angle = cos(angle),
              sin_angle = sin(angle);

        return v*cos_angle + axis*dot(v, axis)*(1.f - cos_angle) + sin_angle*cross(axis, v);
    }

    template <typename Vector3, typename Value = value_t<Vector3>>
    MTS_INLINE
    std::tuple<Value, Value, Value> directions_to_rusinkiewicz(const Vector3 &i, const Vector3 &o) const {
        Vector3 h = normalize(i + o);

        Vector3 n(0, 0, 1);
        Vector3 b = normalize(cross(n, h)),
                t = normalize(cross(b, h));

        Value td = safe_acos(dot(h, i)),
              th = safe_acos(dot(n, h));

        Vector3 i_prj = normalize(i - dot(i, h)*h);
        Value cos_phi_d = clamp(dot(t, i_prj), -1.f, 1.f),
              sin_phi_d = clamp(dot(b, i_prj), -1.f, 1.f);

        Value pd = atan2(sin_phi_d, cos_phi_d);

        return std::make_tuple(pd, th, td);
    }

    MTS_DECLARE_CLASS()
private:
    std::string m_name;
    ScalarFloat m_lambda;
    Interpolator6 m_interpolator;
    ScalarFloat m_alpha_sample;
};

MTS_IMPLEMENT_CLASS_VARIANT(MeasuredPolarized, BSDF)
MTS_EXPORT_PLUGIN(MeasuredPolarized, "Measured polarized material")
NAMESPACE_END(mitsuba)
