#include <string>
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#define CATCH_CONFIG_MAIN
#include <catch/catch.hpp>

#include "MultiComplex/MultiComplex.hpp"

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
	auto f = [](double x) { return x*x*x; };
	typedef std::function<mcx::MultiComplex<double>(const mcx::MultiComplex<double>&)> fcn_t;
	fcn_t ff = [](const mcx::MultiComplex<double>& x) {
		return x*x*x;
	};
	auto xmcx = mcx::MultiComplex<double>(x);
	auto xcx = std::complex<double>(x,x);

	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("call with mcx"){
		return ff(xmcx);
	};
	BENCHMARK("first deriv"){
		return mcx::diff_mcx1(ff, x, 1, false);
	};
}

TEST_CASE("time cos*sin", "[benchmark]")
{
	double x = 8.1;
	auto f = [](double x) { return cos(x) * sin(x); };
	typedef std::function<mcx::MultiComplex<double>(const mcx::MultiComplex<double>&)> fcn_t;
	fcn_t ff = [](const mcx::MultiComplex<double>& x) {
		return cos(x) * sin(x);
	};
	auto xmcx = mcx::MultiComplex<double>(x);

	BENCHMARK("function"){ 		
		return f(x);
	};
	BENCHMARK("call with mcx"){
		return ff(xmcx);
	};
	BENCHMARK("first deriv"){	
		return mcx::diff_mcx1(ff, x, 1, false);
	};
}

TEST_CASE("time sin", "[benchmark],[sin]")
{
	double x = 8.1;
	auto f = [](double x) { return sin(x); };
	typedef std::function<mcx::MultiComplex<double>(const mcx::MultiComplex<double>&)> fcn_t;
	fcn_t ff = [](const mcx::MultiComplex<double>& x) {
		return sin(x);
	};
	auto fgen = [](const std::complex<double>& x) {
		return sin(x);
	};
	auto xcx = std::complex<double>(x, 1e-100);

	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("first deriv"){	
		return mcx::diff_mcx1(ff, x, 1, false);
	};
	BENCHMARK("increment") {
		return mcx::increment(4);
	};
	BENCHMARK("allocate 2 element valarray<double>") {
		return std::valarray<double>(0.0, 2);
	};
	BENCHMARK("get the complex") {
		auto xcx = mcx::MultiComplex<double>(x);
		return xcx.complex();
	}; 
	BENCHMARK("call w/ MultiComplex<double>") {
		auto xcx = mcx::MultiComplex<double>(x);
		return ff(xcx);
	};
}

TEST_CASE("time ln(x)", "[benchmark]")
{
	double x = 8.1;
	auto f = [](double x) { return log(x); };
	BENCHMARK("function") {
		return f(x);
	};
	BENCHMARK("first deriv")
	{
		typedef std::function<mcx::MultiComplex<double>(const mcx::MultiComplex<double>&)> fcn_t;
		fcn_t ff = [](const mcx::MultiComplex<double>& x) {
			return log(x);
		};
		return mcx::diff_mcx1(ff, x, 1, false);
	};
}