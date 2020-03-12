#include "MultiComplex/MultiComplex.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>

namespace py = pybind11;

void init_MultiComplex(py::module &m){
    typedef MultiComplex<double> MCD;

    using TN = double;
    m.def("diff_mcx1", py::overload_cast<const std::function<MultiComplex<TN>(const MultiComplex<TN>&)>&, TN, int, bool>(&diff_mcx1<TN>), 
        py::arg("f"), py::arg("x"), py::arg("numderiv"), py::arg("and_val") = false);
    using tuplefunction = std::function<std::tuple<MultiComplex<TN>, MultiComplex<TN>>(const MultiComplex<TN>&)>;
    m.def("diff_mcx1", py::overload_cast<const tuplefunction&, TN, int, bool>(&diff_mcx1<TN>), 
        py::arg("f"), py::arg("x"), py::arg("numderiv"), py::arg("and_val") = false);
    m.def("diff_mcxN", &diff_mcxN<double>);

    py::class_<MCD>(m, "MultiComplex")
        .def(py::init<const std::complex<double> &>())
        .def(py::init<const std::valarray<double>&>())
        // Operators for pairs of MultiComplex instances
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(py::self / py::self)
        // Bidirectional operators working with double variables
        .def(py::self + double())
        .def(double() +  py::self)
        .def(py::self - double())
        .def(double() - py::self)
        .def(py::self / double())
        .def(double() / py::self)
        .def(py::self * double())
        .def(double() * py::self)
        // Unary operators
        .def(-py::self)
        // Comparison operators
        .def("__lt__", [](const MCD& a, const double b) { return a.real() < b; }, py::is_operator())
        .def("__lt__", [](const MCD& a, const MCD& b) { return a.real() < b.real(); }, py::is_operator())
        .def("__gt__", [](const MCD& a, const double b) { return a.real() > b; }, py::is_operator())
        .def("__gt__", [](const MCD& a, const MCD& b) { return a.real() > b.real(); }, py::is_operator())

        .def("get_coef", &MCD::get_coef)
        .def("dim", &MCD::dim)
        .def("real", &MCD::real)
        .def("complex", &MCD::complex)
        .def("__pow__", [](MCD& mc, int exponent) { return mc.pow(exponent); })
        .def("__pow__", [](MCD& mc, double exponent) { return mc.pow(exponent); })
        .def("__pow__", [](MCD& mc, const std::vector<double> & exponents) { 
            std::vector<MCD> o;
            for (auto exponent: exponents){
                o.push_back(exp(exponent*log(mc))); 
            }
            return o;
        })
        .def("__getitem__", [](MCD& mc, int index) { return mc[index]; })
        .def("__setitem__", [](MCD& mc, int index, double value) { mc.set_coef(index,value); })
        .def("sin", [](MCD&mc) { return sin(mc); })
        .def("cos", [](MCD& mc) { return cos(mc); })
        .def("sinh", [](MCD& mc) { return sinh(mc); })
        .def("cosh", [](MCD& mc) { return cosh(mc); })
        .def("exp", [](MCD& mc) { return exp(mc); })
        .def("log", [](MCD& mc) { return log(mc); })
        ;
}

PYBIND11_MODULE(pymcx, m) {
    m.doc() = " A Python wrapper of C++ library for working with multicomplex mathematics";
    init_MultiComplex(m);
}