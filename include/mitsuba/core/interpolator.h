#pragma once

#include <mitsuba/core/fwd.h>
#include <mitsuba/core/math.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/core/util.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief Implements linear interpolation for data of arbitrary dimension.
 *
 * The input array should have dimensions <tt>N0 x N1 x ... x Nn</tt>
 * (where the last dimension is contiguous in memory). <tt>param_res</tt> should
 * be set to <tt>{ N0, N1, ..., Nn }</tt>, and <tt>param_values</tt> should
 * contain the parameter values where the distribution is discretized.
 * Linear interpolation is used when evaluating the function for in-between
 * values.
 *
 *
 * \remark The Python API exposes explicitly instantiated versions of this
 * class named LinearInterpolator1, LinearInterpolator2, etc.
 */
template <typename Float_, size_t Dimension_ = 1> class LinearInterpolator {
public:
    static constexpr size_t Dimension = Dimension_;
    using Float                       = Float_;
    using UInt32                      = uint32_array_t<Float>;
    using Mask                        = mask_t<Float>;
    using ScalarFloat                 = scalar_t<Float>;
    using FloatStorage                = DynamicBuffer<Float>;

    LinearInterpolator() = default;

    LinearInterpolator(const ScalarFloat *data,
                       const std::array<uint32_t, Dimension> &param_res = { },
                       const std::array<const ScalarFloat *, Dimension> &param_values = { }) {
        m_slices = 1;
        for (int i = (int) Dimension - 1; i >= 0; --i) {
            if (param_res[i] < 1)
                Throw("LinearInterpolator(): parameter resolution must be >= 1!");
            m_param_values[i] = FloatStorage::copy(param_values[i], param_res[i]);
            m_param_strides[i] = param_res[i] > 1 ? m_slices : 0;
            m_slices *= param_res[i];
        }

        m_data = empty<FloatStorage>(m_slices);
        m_data.managed();
        memcpy(m_data.data(), data, sizeof(ScalarFloat)*m_slices);
    }

    // Look up parameter-related indices and weights
    UInt32 interpolate_weights(const Float *param, Float *param_weight,
                               Mask active) const {
        ENOKI_MARK_USED(param);

        if constexpr (Dimension > 0) {
            MTS_MASK_ARGUMENT(active);

            UInt32 slice_offset = zero<UInt32>();
            for (size_t dim = 0; dim < Dimension; ++dim) {
                if (unlikely(m_param_values[dim].size() == 1)) {
                    param_weight[2 * dim] = 1.f;
                    param_weight[2 * dim + 1] = 0.f;
                    continue;
                }

                UInt32 param_index = math::find_interval(
                    (uint32_t) m_param_values[dim].size(),
                    [&](UInt32 idx) ENOKI_INLINE_LAMBDA {
                        return gather<Float>(m_param_values[dim], idx, active) <
                               param[dim];
                    });

                Float p0 = gather<Float>(m_param_values[dim], param_index, active),
                      p1 = gather<Float>(m_param_values[dim], param_index + 1, active);

                param_weight[2 * dim + 1] =
                    clamp((param[dim] - p0) / (p1 - p0), 0.f, 1.f);
                param_weight[2 * dim] = 1.f - param_weight[2 * dim + 1];
                slice_offset += m_param_strides[dim] * param_index;
            }

            return slice_offset;
        } else {
            ENOKI_MARK_USED(param);
            ENOKI_MARK_USED(param_weight);
            ENOKI_MARK_USED(active);
            return 0u;
        }
    }

    /**
     * \brief Evaluate the function at \c param.
     */
    Float eval(const Float *param = nullptr,
               Mask active = true) const {
        MTS_MASK_ARGUMENT(active);
        Float param_weight[2 * DimensionInt];
        UInt32 slice_offset = interpolate_weights(param, param_weight, active);
        return lookup<Dimension>(m_data.data(), slice_offset, param_weight, active);
    }

    std::string to_string() const {
        std::ostringstream oss;
        oss << "LinearInterpolator" << Dimension << "[" << std::endl;
        oss << "  param_size = [";
        for (size_t i = 0; i<Dimension; ++i) {
            if (i != 0)
                oss << ", ";
            oss << m_param_values[i].size();
        }
        oss << "]," << std::endl
            << "  param_strides = [";
        for (size_t i = 0; i<Dimension; ++i) {
            if (i != 0)
                oss << ", ";
            oss << m_param_strides[i];
        }
        oss << "]," << std::endl;
        oss << "  storage = { " << m_slices << " slice" << (m_slices > 1 ? "s" : "")
            << ", ";
        oss << util::mem_string(m_slices * sizeof(ScalarFloat)) << " }" << std::endl
            << "]";
        return oss.str();
    }

protected:
    template <size_t Dim = Dimension>
    MTS_INLINE Float lookup(const ScalarFloat *data,
                            UInt32 i0,
                            const Float *param_weight,
                            Mask active) const {
        if constexpr (Dim != 0) {
            UInt32 i1 = i0 + m_param_strides[Dim - 1];

            Float w0 = param_weight[2 * Dim - 2],
                  w1 = param_weight[2 * Dim - 1],
                  v0 = lookup<Dim - 1>(data, i0, param_weight, active),
                  v1 = lookup<Dim - 1>(data, i1, param_weight, active);

            return fmadd(v0, w0, v1 * w1);
        } else {
            ENOKI_MARK_USED(param_weight);
            return gather<Float>(data, i0, active);
        }
    }

protected:
#if !defined(_MSC_VER)
    static constexpr size_t DimensionInt = Dimension;
#else
    static constexpr size_t DimensionInt = (Dimension != 0) ? Dimension : 1;
#endif

    /// Stride per parameter in units of sizeof(ScalarFloat)
    uint32_t m_param_strides[DimensionInt];

    /// Discretization of each parameter domain
    FloatStorage m_param_values[DimensionInt];

    /// Total number of slices
    uint32_t m_slices;

    /// Function values
    FloatStorage m_data;
};

NAMESPACE_END(mitsuba)
