#include <string>
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#define CATCH_CONFIG_MAIN
#include <catch/catch.hpp>

#include "MultiComplex/MultiComplex.hpp"

using mcxf = std::function<mcx::MultiComplex<double>(const mcx::MultiComplex<double>&)>;

TEST_CASE("Time some bit things", "[bits]")
{
	BENCHMARK("exp2i") {
		return mcx::exp2i(3); // 2^3
	};
	BENCHMARK("ln2i") {
		return mcx::log2i(8UL); // ln2(8), for integer inputs
	};
}

TEST_CASE("time x^3", "[benchmark]")
{
	double x = 8.1;
	auto f = [](const auto& x) { return x*x*x; };

	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("call with mcx"){
		auto xmcx = mcx::MultiComplex<double>(x);
		return static_cast<mcxf>(f)(xmcx);
	};
	BENCHMARK("first deriv"){
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 1, false);
	};
}

TEST_CASE("time x^3, but with pow(x,int)", "[benchmark]")
{
	double x = 8.1;
	auto f = [](const auto& x) { return pow(x, 3); };

	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("call with mcx") {
		auto xmcx = mcx::MultiComplex<double>(x);
		return static_cast<mcxf>(f)(xmcx);
	};
	BENCHMARK("first deriv") {
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 1, false);
	};
}

TEST_CASE("time x^3, but with pow(x,double)", "[benchmark]")
{
	double x = 8.1;
	auto f = [](const auto& x) { return pow(x, 3.0); };

	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("call with mcx") {
		auto xmcx = mcx::MultiComplex<double>(x);
		return static_cast<mcxf>(f)(xmcx);
	};
	BENCHMARK("first deriv") {
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 1, false);
	};
}

TEST_CASE("time 1/x", "[benchmark]")
{
	double x = 8.1;
	auto f = [](const auto& x) { return 1.0/x; };

	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("call with mcx") {
		auto xmcx = mcx::MultiComplex<double>(x);
		return static_cast<mcxf>(f)(xmcx);
	};
	BENCHMARK("first deriv") {
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 1, false);
	};
}

TEST_CASE("time sin(x)", "[benchmark],[sin]")
{
	double x = 8.1;
	auto f = [](const auto& x) { return sin(x); };

	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("first deriv"){	
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 1, false);
	};
}

TEST_CASE("time cosh(x)", "[benchmark],[cosh]")
{
	double x = 8.1;
	auto f = [](const auto& x) { return cosh(x); };

	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("first deriv") {
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 1, false);
	};
}

TEST_CASE("time exp(x)", "[benchmark]")
{
	double x = 8.1;
	auto f = [](const auto &x) { return exp(x); };
	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("first deriv"){
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 1, false);
	};
}

TEST_CASE("time ln(x)", "[benchmark]")
{
	double x = 8.1;
	auto f = [](const auto& x) { return log(x); };
	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("first deriv"){
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 1, false);
	};
}

TEST_CASE("time cos(x)*sin(x)", "[benchmark]")
{
	double x = 8.1;
	auto f = [](const auto& x) { return cos(x) * sin(x); };
	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("first deriv") {
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 1, false);
	};
}

TEST_CASE("time alphar", "[benchmark]")
{
	double x = 8.1;
	
	// Critical point values for argon
	double R = 8.31446261815324,
		   Tc = 150.687, // K
		   pc = 4863000.0, // Pa
		   a = (27.0/64.0)*pow(R*Tc, 2)/pc, 
		   b = (1.0/8.0)*(R*Tc)/pc,
		   T = 298.15;
	auto f = [&](const auto& rho) { return -log(1.0-b*rho)-a*rho/(R*T); };

	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("first deriv") {
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 1, false);
	};
	BENCHMARK("second deriv") {
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 2, false);
	};
	BENCHMARK("third deriv") {
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 3, false);
	};
	BENCHMARK("fourth deriv") {
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 4, false);
	};
	BENCHMARK("fifth deriv") {
		return mcx::diff_mcx1(static_cast<mcxf>(f), x, 5, false);
	};
}