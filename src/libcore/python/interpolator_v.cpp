#include <mitsuba/python/python.h>
#include <mitsuba/core/interpolator.h>
#include <pybind11/numpy.h>

template <typename Interp> auto bind_interpolator(py::module &m,
                                                  const char *name,
                                                  const char *doc,
                                                  const char *doc_constructor,
                                                  const char *doc_eval) {
    using Float                = typename Interp::Float;
    using ScalarFloat          = scalar_t<Float>;
    using NumPyArray           = py::array_t<ScalarFloat, py::array::c_style | py::array::forcecast>;
    using Vector2f             = enoki::Array<Float, 2>;
    using ScalarVector2u       = enoki::Array<uint32_t, 2>;
    using Mask                 = mask_t<Float>;

    // py::object zero = py::cast(enoki::zero<Array<ScalarFloat, Interp::Dimension>>());

    auto constructor =
        py::init([](const NumPyArray &data,
                    const std::array<std::vector<ScalarFloat>, Interp::Dimension> &param_values_in) {
            if (data.ndim() != Interp::Dimension)
                throw std::domain_error("'data' array has incorrect dimension");

            std::array<uint32_t, Interp::Dimension> param_res;
            std::array<const ScalarFloat *, Interp::Dimension> param_values;

            for (size_t i = 0; i < Interp::Dimension; ++i) {
                if (param_values_in[i].size() != (size_t) data.shape(i))
                    throw std::domain_error(
                        "'param_values' array has incorrect dimension");
                param_values[i] = param_values_in[i].data();
                param_res[i]    = param_values_in[i].size();
            }

            return Interp(data.data(),
                          param_res, param_values);
        });

    auto interpolator = py::class_<Interp>(m, name, py::module_local(), doc);

    interpolator.def(std::move(constructor),
         "data"_a, "param_values"_a, doc_constructor);

    interpolator.def("eval",
         vectorize([](const Interp *i,
                      const Array<Float, Interp::Dimension> &param,
                      Mask active) {
            return i->eval(param.data(), active);
         }),
         "param"_a, "active"_a = true, doc_eval);
}


MTS_PY_EXPORT(Interpolator) {
    MTS_PY_IMPORT_TYPES()

    bind_interpolator<LinearInterpolator<Float, 1>>(m, "LinearInterpolator1",
                                                    "", "", "");
    bind_interpolator<LinearInterpolator<Float, 2>>(m, "LinearInterpolator2",
                                                    "", "", "");
    bind_interpolator<LinearInterpolator<Float, 3>>(m, "LinearInterpolator3",
                                                    "", "", "");
    bind_interpolator<LinearInterpolator<Float, 6>>(m, "LinearInterpolator6",
                                                    "", "", "");
}

// TODO: make docstrings doesn't work anymore at the moment, replace doc strings above:
// D(LinearInterpolator),
// D(LinearInterpolator, LinearInterpolator),
// D(LinearInterpolator, eval)
