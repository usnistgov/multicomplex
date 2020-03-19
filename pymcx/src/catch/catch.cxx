#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"

#include "MultiComplex/MultiComplex.hpp"

/// Return true if all elements of argument are finite and representable in double precision
template<typename VEC>
inline bool all_finite(const VEC& vec) {
    bool all_ok = true;
    for (auto v : vec) {
        if (!std::isfinite(v)) {
            all_ok = false;
        }
    }
    return all_ok;
}

TEST_CASE("Myslice", "[myslice]") {
    std::slice even_increments = myslice(0, 10, 2);
    std::slice uneven_increments = myslice(0, 10, 3);
    REQUIRE(even_increments.size() == 5);
    REQUIRE(uneven_increments.size() == 3);
}

TEST_CASE("log2i", "[log2i]") {
    CHECK(log2i(1) == 0); 
    CHECK(log2i(2) == 1);
    CHECK(log2i(4) == 2);
    CHECK_THROWS(log2i(7));
}

TEST_CASE("real(MCX)", "[1D]") {
    MultiComplex<double> z{ std::complex<double>{3.0, 1e-100} };
    REQUIRE(z.real() == 3.0);
}

TEST_CASE("MCX < double", "[1D]") {
    MultiComplex<double> z{ std::complex<double>{3.0, 1e-100} };
    REQUIRE(z < 4.0);
}

TEST_CASE("exp(-big)", "[1D]") {
    MultiComplex<double> n{{-100000, 1e-50, 1e-50, 1e-150}}, 
                         nexp = exp(n);
    REQUIRE(all_finite(nexp.get_coef()) == true);
}

TEST_CASE("1/n derivs","[1D]") {
    typedef std::function<MultiComplex<double>(const MultiComplex<double>&)> fcn_t;
    fcn_t ff = [](const MultiComplex<double>& z) {
        return 1.0 / z;
    };
    double x = 0.1234; 
    int numderiv = 6;
    auto fo = diff_mcx1(ff, x, numderiv);
    std::vector<double> exacts;
    auto factorial = [](int n) {
        return std::tgamma(n + 1);
    };
    for (int n = 1; n <= numderiv; ++n) {
        exacts.push_back(pow(-1,n)*factorial(n)/pow(x,n+1));
    }
    double abs_rel_errs = 0;
    for (auto i = 0; i < numderiv; ++i) {
        abs_rel_errs += std::abs((fo[i] - exacts[i])/exacts[i]);
    }
    REQUIRE(abs_rel_errs < 1e-12);
}

TEST_CASE("x^4 derivs", "[1D]") {
    typedef std::function<MultiComplex<double>(const MultiComplex<double>&)> fcn_t;
    fcn_t ff = [](const MultiComplex<double>& z) {
        return z.pow(4);
    };
    double x = 0.1234;
    int numderiv = 6;
    auto fo = diff_mcx1(ff, x, numderiv);
    std::vector<double> exacts(numderiv);
    exacts[0] = 4*std::pow(x,3);
    exacts[1] = 12*std::pow(x, 2);
    exacts[2] = 24*std::pow(x, 1);
    exacts[3] = 24;
    exacts[4] = 0;
    exacts[5] = 0;
    double abs_errs = 0;
    for (auto i = 0; i < numderiv; ++i) {
        abs_errs += std::abs((fo[i] - exacts[i]));
    }
    REQUIRE(abs_errs < 1e-12);
}

TEST_CASE("x^4 derivs; 4 double", "[1D]") {
    typedef std::function<MultiComplex<double>(const MultiComplex<double>&)> fcn_t;
    fcn_t ff = [](const MultiComplex<double>& z) {
        return z.pow(4.0);
    };
    double x = 0.1234;
    int numderiv = 6;
    auto fo = diff_mcx1(ff, x, numderiv);
    std::vector<double> exacts(numderiv);
    exacts[0] = 4.0*std::pow(x, 3.0);
    exacts[1] = 12.0*std::pow(x, 2.0);
    exacts[2] = 24.0*std::pow(x, 1.0);
    exacts[3] = 24.0;
    exacts[4] = 0;
    exacts[5] = 0;
    double abs_errs = 0;
    for (auto i = 0; i < numderiv; ++i) {
        abs_errs += std::abs((fo[i] - exacts[i]));
    }
    REQUIRE(abs_errs < 1e-12);
}

TEST_CASE("x^4 derivs, returned as tuple", "[1D]") {
    using fcn_t = std::function<std::tuple<MultiComplex<double>, MultiComplex<double>>(const MultiComplex<double>&)>;
    fcn_t ff = [](const MultiComplex<double>& z) {
        return std::make_tuple(z.pow(4.0), z.pow(4.0));
    };
    double x = 0.1234;
    int numderiv = 6;
    auto [fo, err] = diff_mcx1(ff, x, numderiv); // A dummy error output for testing
    std::vector<double> exacts(numderiv);
    exacts[0] = 4.0 * std::pow(x, 3.0);
    exacts[1] = 12.0 * std::pow(x, 2.0);
    exacts[2] = 24.0 * std::pow(x, 1.0);
    exacts[3] = 24.0;
    exacts[4] = 0;
    exacts[5] = 0;
    double abs_errs = 0;
    for (auto i = 0; i < numderiv; ++i) {
        abs_errs += std::abs((fo[i] - exacts[i]));
    }
    REQUIRE(abs_errs < 1e-12);
}

TEST_CASE("Higher derivatives","[ND]") {
    auto func = [](const std::vector<MultiComplex<double>>& zs) {
        return cos(zs[0]) * sin(zs[1]) * exp(zs[2]);
    };
    SECTION("110"){
        std::vector<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 1, 1, 0 };
        auto exact = -sin(xs[0]) * cos(xs[1]) * exp(xs[2]);
        auto num = diff_mcxN<double>(func, xs, order);
        auto abs_err = std::abs(exact-num);
        REQUIRE(abs_err < 1e-15);
    }
    SECTION("114") {
        std::vector<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 1, 1, 4 };
        auto exact = -sin(xs[0]) * cos(xs[1]) * exp(xs[2]);
        auto num = diff_mcxN<double>(func, xs, order);
        auto abs_err = std::abs(exact - num);
        REQUIRE(abs_err < 1e-15);
    }
    SECTION("414") {
        std::vector<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 4, 1, 4 };
        auto exact = cos(xs[0]) * cos(xs[1]) * exp(xs[2]);
        auto num = diff_mcxN<double>(func, xs, order);
        auto abs_err = std::abs(exact - num);
        REQUIRE(abs_err < 1e-15);
    }
    SECTION("Bad") {
        std::vector<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 4 };
        CHECK_THROWS(diff_mcxN<double>(func, xs, order));
    }
}