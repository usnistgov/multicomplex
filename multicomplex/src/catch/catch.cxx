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

template<typename T>
inline bool almost_equal_complex(const std::complex<T>& v1, const std::complex<T>& v2, double re_tol, double im_tol) {
    auto re_diff = std::abs(v1.real() - v2.real());
    auto im_diff = std::abs(v1.imag() - v2.imag());
    return re_diff < re_tol && im_diff < im_tol;
}


TEST_CASE("Myslice", "[myslice]") {
    std::slice even_increments = mcx::myslice(0, 10, 2);
    std::slice uneven_increments = mcx::myslice(0, 10, 3);
    REQUIRE(even_increments.size() == 5);
    REQUIRE(uneven_increments.size() == 3);
}

TEST_CASE("log2i", "[log2i]") {
    CHECK(mcx::log2i(1) == 0);
    CHECK(mcx::log2i(2) == 1);
    CHECK(mcx::log2i(4) == 2);
    CHECK(mcx::log2i(8) == 3);
    CHECK(mcx::log2i(16) == 4);
    CHECK_THROWS(mcx::log2i(7));
}

TEST_CASE("real(MCX)", "[1D]") {
    mcx::MultiComplex<double> z{ std::complex<double>{3.0, 1e-100} };
    REQUIRE(z.real() == 3.0);
}

TEST_CASE("inplace operations", "[1D]") {
    mcx::MultiComplex<double> z{ std::complex<double>{3.0, 1e-100} };
    z += 3.7;
    z -= 3.7;
    z *= 1.8;
    z /= 1.8;
    REQUIRE(z.real() == 3.0);
}


TEST_CASE("MCX < double", "[1D]") {
    mcx::MultiComplex<double> z{ std::complex<double>{3.0, 1e-100} };
    REQUIRE(z < 4.0);
}

TEST_CASE("exp(-big)", "[1D]") {
    mcx::MultiComplex<double> n{{-100000, 1e-50, 1e-50, 1e-150}},
                         nexp = exp(n);
    REQUIRE(all_finite(nexp.get_coef()) == true);
}

TEST_CASE("square of negative number", "[ops]") {
    std::complex<double> ee = std::complex<double>(-0.1, 1e-100);
    mcx::MultiComplex<double> mce(ee);
    CHECK(almost_equal_complex(ee*ee, (mce*mce).complex(), 1e-14, 1e-16));
    CHECK(almost_equal_complex(ee*ee, (mce.pow(2)).complex(), 1e-14, 1e-16));
    CHECK(almost_equal_complex(ee*ee, (mce.pow(2.0)).complex(), 1e-14, 1e-16));
    CHECK(almost_equal_complex(ee*ee, (pow(mce, 2.0)).complex(), 1e-14, 1e-16));
    CHECK_THROWS(mce.pow(2.1));
}

TEST_CASE("Integer and non-integer powers", "[ops]") {
    std::complex<double> ee = std::complex<double>(0.1, 1e-100);
    mcx::MultiComplex<double> mce(ee);
    CHECK(almost_equal_complex(ee * ee, (mce * mce).complex(), 1e-14, 1e-16));
    CHECK(almost_equal_complex(ee * ee, (mce.pow(2)).complex(), 1e-14, 1e-16));
    CHECK(almost_equal_complex(ee * ee, (mce.pow(2.0)).complex(), 1e-14, 1e-16));
    CHECK(almost_equal_complex(ee * ee, (pow(mce, 2.0)).complex(), 1e-14, 1e-16));
    CHECK(almost_equal_complex(pow(ee, 2.1), (mce.pow(2.1)).complex(), 1e-14, 1e-16));
}

TEST_CASE("1/x derivs","[1D]") {
    typedef std::function<mcx::MultiComplex<double>(const mcx::MultiComplex<double>&)> fcn_t;
    fcn_t ff = [](const mcx::MultiComplex<double>& z) {
        return 1.0 / z;
    };
    double x = 0.1234; 
    int numderiv = 4;
    auto fo = mcx::diff_mcx1(ff, x, numderiv, false /* and_val */);
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
    typedef std::function<mcx::MultiComplex<double>(const mcx::MultiComplex<double>&)> fcn_t;
    fcn_t ff = [](const mcx::MultiComplex<double>& z) {
        return z.pow(4);
    };
    double x = 0.1234;
    int numderiv = 6;
    auto fo = mcx::diff_mcx1(ff, x, numderiv);
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
    REQUIRE(abs_errs < 1e-10);
}

TEST_CASE("x^4 derivs; 4 double", "[1D]") {
    typedef std::function<mcx::MultiComplex<double>(const mcx::MultiComplex<double>&)> fcn_t;
    fcn_t ff = [](const mcx::MultiComplex<double>& z) {
        return z.pow(4.0);
    };
    double x = 0.1234;
    int numderiv = 6;
    auto fo = mcx::diff_mcx1(ff, x, numderiv);
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
    REQUIRE(abs_errs < 1e-10);
}

TEST_CASE("x^4 derivs, returned as tuple", "[1D]") {
    using fcn_t = std::function<std::tuple<mcx::MultiComplex<double>, mcx::MultiComplex<double>>(const mcx::MultiComplex<double>&)>;
    fcn_t ff = [](const mcx::MultiComplex<double>& z) {
        return std::make_tuple(z.pow(4.0), z.pow(4.0));
    };
    double x = 0.1234;
    int numderiv = 6;
    auto [fo, err] = mcx::diff_mcx1(ff, x, numderiv); // A dummy error output for testing
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
    REQUIRE(abs_errs < 1e-10);
}

TEST_CASE("Higher derivatives","[ND]") {
    using fcn_t = std::function<mcx::MultiComplex<double>(const std::vector<mcx::MultiComplex<double>>&)>;
    fcn_t func = [](const std::vector<mcx::MultiComplex<double>>& zs){
        return cos(zs[0]) * sin(zs[1]) * exp(zs[2]);
    };
    SECTION("110"){
        std::vector<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 1, 1, 0 };
        auto exact = -sin(xs[0]) * cos(xs[1]) * exp(xs[2]);
        auto num = diff_mcxN(func, xs, order);
        auto abs_err = std::abs(exact-num);
        REQUIRE(abs_err < 1e-15);
    }
    SECTION("114") {
        std::vector<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 1, 1, 4 };
        auto exact = -sin(xs[0]) * cos(xs[1]) * exp(xs[2]);
        auto num = diff_mcxN(func, xs, order);
        auto abs_err = std::abs(exact - num);
        REQUIRE(abs_err < 1e-15);
    }
    SECTION("414") {
        std::vector<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 4, 1, 4 };
        auto exact = cos(xs[0]) * cos(xs[1]) * exp(xs[2]);
        auto num = diff_mcxN(func, xs, order);
        auto abs_err = std::abs(exact - num);
        REQUIRE(abs_err < 1e-15);
    }
    SECTION("Bad") {
        std::vector<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 4 };
        CHECK_THROWS(diff_mcxN(func, xs, order));
    }
}

TEST_CASE("Higher derivatives w/ std::valarray", "[ND]") {
    using fcn_t = std::function<mcx::MultiComplex<double>(const std::valarray<mcx::MultiComplex<double>>&)>;
    fcn_t func = [](const std::valarray<mcx::MultiComplex<double>>& zs) {
        return cos(zs[0]) * sin(zs[1]) * exp(zs[2]);
    };
    SECTION("110") {
        std::valarray<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 1, 1, 0 };
        auto exact = -sin(xs[0]) * cos(xs[1]) * exp(xs[2]);
        auto num = diff_mcxN(func, xs, order);
        auto abs_err = std::abs(exact - num);
        REQUIRE(abs_err < 1e-15);
    }
    SECTION("114") {
        std::valarray<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 1, 1, 4 };
        auto exact = -sin(xs[0]) * cos(xs[1]) * exp(xs[2]);
        auto num = diff_mcxN(func, xs, order);
        auto abs_err = std::abs(exact - num);
        REQUIRE(abs_err < 1e-15);
    }
    SECTION("414") {
        std::valarray<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 4, 1, 4 };
        auto exact = cos(xs[0]) * cos(xs[1]) * exp(xs[2]);
        auto num = diff_mcxN(func, xs, order);
        auto abs_err = std::abs(exact - num);
        REQUIRE(abs_err < 1e-15);
    }
    SECTION("Bad") {
        std::valarray<double> xs = { 0.1234, 20.1234, -4.1234 };
        std::vector<int> order = { 4 };
        CHECK_THROWS(diff_mcxN(func, xs, order));
    }
}

TEST_CASE("Hessian", "[ND]") {
    using fcn_t = std::function< mcx::MultiComplex<double>(const std::valarray<mcx::MultiComplex<double>>&)>;
    fcn_t func = [](const std::valarray<mcx::MultiComplex<double>>& zs) -> mcx::MultiComplex<double> {
        return cos(zs[0]) * sin(zs[1]);
    };
    double x =0.1234, y=20.1234, z=-4.1234;
    std::valarray<std::valarray<double>> Hessianexact = {{-sin(y) * cos(x), -sin(x) * cos(y)}, {-sin(x) * cos(y), -sin(y) * cos(x)}};
    
    SECTION("Hessian multiple") {
        std::valarray<double> pt = {x, y};
        using mattype = std::valarray<std::valarray<double>>;
        using functype = decltype(func);
        auto H = mcx::get_Hessian<mattype, functype, std::valarray<double>, mcx::HessianMethods::Multiple>(func, pt);
        for (auto i = 0; i < 2; ++i) {
            CHECK(std::abs((H[i]-Hessianexact[i]).max()) < 1e-14);
        }
        int rr =0;
    }
}
