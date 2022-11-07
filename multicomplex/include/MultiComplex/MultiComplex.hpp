#ifndef MULTICOMPLEX_HPP
#define MULTICOMPLEX_HPP

#include <chrono>
#include <vector>
#include <iostream>
#include <complex>
#include <atomic>
#include <valarray>
#include <array>
#include <tuple>
#include <functional>
#include <numeric>
#include <cmath>

#if !defined(MULTICOMPLEX_NO_MULTIPRECISION)
#ifdef __has_include                                           // Check if __has_include is present
#if __has_include(<boost/multiprecision/cpp_bin_float.hpp>)  // Check for presence of boost/multiprecision
#define  BOOST_MULTIPRECISION_FOUND
#include "boost/multiprecision/cpp_bin_float.hpp"
#include "boost/multiprecision/cpp_complex.hpp"
#endif
#endif
#endif

namespace mcx {

// A custom slice, more like Python, where the arguments to the slice are start, stop, increment
// rather than start, number of times, increment
inline auto myslice(std::size_t start, std::size_t stop, std::size_t increment) {
    return std::slice(start, (stop - start) / increment, increment);
}

/// Raise 2 to the power i, where i is an integer
// see https://stackoverflow.com/a/5345369 for explanation
inline int exp2i(int i) {
    return 1 << i;
}

/// Log2i, where i is an integer
/// find the exponent j that yields 2^j = i if possible
/// This call could be further optimized in a non-portable way: https://stackoverflow.com/questions/671815/what-is-the-fastest-most-efficient-way-to-find-the-highest-set-bit-msb-in-an-i
inline int log2i(unsigned int i) {
    unsigned int counter = 0; // How many times you need to multiply by 2...
    for (unsigned int j = 1; j <= i; j *= 2) {
        if (j == i) {
            return counter;
        }
        counter++;
    }
    throw std::invalid_argument("The argument to log2i " + std::to_string(i) + " is not a power of 2");
}

inline double increment(std::size_t l) {
    return ldexp(1.0, -664 / static_cast<int>(l));
}

// Scalar versions of math functions that help 
// the compiler with overload resolution
template<typename T> inline auto scalar_cosh(const T& x) { if constexpr (std::is_arithmetic_v<T>) { return ::cosh(x); } else { return cosh(x); } };
template<typename T> inline auto scalar_sinh(const T& x) { if constexpr (std::is_arithmetic_v<T>) { return ::sinh(x); } else { return sinh(x); } };
template<typename T> inline auto scalar_cos(const T& x) { if constexpr (std::is_arithmetic_v<T>) { return ::cos(x); } else { return cos(x); } };
template<typename T> inline auto scalar_sin(const T& x) { if constexpr (std::is_arithmetic_v<T>) { return ::sin(x); } else { return sin(x); } };
template<typename T> inline auto scalar_exp(const T& x) { if constexpr (std::is_arithmetic_v<T>) { return ::exp(x); } else { return exp(x); } };
template<typename T> inline auto scalar_log(const T& x) { if constexpr (std::is_arithmetic_v<T>) { return ::log(x); } else { return log(x); } };
template<typename T> inline auto scalar_abs(const T& x) { if constexpr (std::is_arithmetic_v<T>) { return std::abs(x); } else { return abs(x); } };
template<typename TN>
TN scalar_pow(const TN& x, int e) {
    if constexpr (std::is_floating_point<TN>::value) {
        // base is a number (float, double, long double)
        return std::pow(x, e);
    }
    else {
        // This is something more interesting
        return pow(x, e);
    }
}

template<typename T >
struct MultiComplex
{
    int m_d = 0;
    std::valarray<T> coef;

    /// Default constructor
    MultiComplex() : m_d(0) { coef = {}; };

    /// This is a "normal" floating point number (re+i*0), stored as a complex number
    /// All floating point values *other* than the numerical type used in the class are enabled as candidate types 
    template<typename TI, typename std::enable_if<!std::is_same_v<TI, T> && std::is_arithmetic_v<TI>>::type* = nullptr>
    MultiComplex(const TI& re) : m_d(1) { coef.resize(2); coef[0] = re, coef[1] = 0; };

    /// This is a "normal" floating point number (re+i*0), stored as a complex number
    MultiComplex(const T& re) : m_d(1) { coef.resize(2); coef[0] = re, coef[1] = 0; };

    // This is a "normal" complex number, stored as such 
    MultiComplex(const std::complex<T>& c) : m_d(1) { coef = { c.real(), c.imag() }; };

    /// Constructor in which a vector of values are provided
    /// These are the values of the terminals, and the length of the values will be 2^d, 
    /// where d is the dimension (or level) of the MultiComplex object
    MultiComplex(const std::valarray<T>& v) {
        coef = v;
        m_d = log2i(static_cast<int>(v.size()));
    }
    MultiComplex(std::valarray<T>&& v) {
        m_d = log2i(static_cast<int>(v.size())); 
        coef = std::move(v);
    }

    /// The dimension of the multicomplex number
    int dim() const { return m_d; }

    /// The "most real" component; infinity if not initialized yet
    T real() const {
        return (m_d > 0) ? coef[0] : std::numeric_limits<T>::infinity();
    }

    auto complex() const {
#if defined(BOOST_MULTIPRECISION_FOUND)
        if constexpr (boost::multiprecision::is_number<T>::value) { // Is something from boost::multiprecision
            return boost::multiprecision::cpp_complex<std::numeric_limits<T>::max_digits10>(coef[0], coef[1]);
        }
        else {
            return std::complex<T>(coef[0], coef[1]);
        }
#else
        return std::complex<T>(coef[0], coef[1]);
#endif
    }

    /// Get a const reference to the coefficients
    const std::valarray<T>& get_coef() const { return coef; }

    /// Unary negation parameter (-this)
    const MultiComplex operator-() const {
        decltype(coef) new_coef = -coef;
        return new_coef;
    }
    /// Right addition by a scalar (this+r)
    template<typename TI, typename std::enable_if<!std::is_same_v<TI, T> && std::is_arithmetic_v<TI>>::type* = nullptr>
    const MultiComplex operator+(TI r) const {
        decltype(coef) new_coef = coef; new_coef[0] += r;
        return new_coef;
    }
    /// Inplace addition of a scalar (this += r)
    template<typename TI, typename std::enable_if<!std::is_same_v<TI, T> && std::is_arithmetic_v<TI>>::type* = nullptr>
    auto& operator+=(const TI& r) {
        coef[0] += r;
        return *this;
    }
    /// Inplace addition of a scalar (this += r)
    auto& operator+=(const T& r) {
        coef[0] += r;
        return *this;
    }
    auto& operator+=(const MultiComplex& w) {
        auto n = *this + w;
        std::swap(*this, n);
        return *this;
    }
    /// Right subtraction by a scalar (this-r)
    template<typename TI, typename std::enable_if<!std::is_same_v<TI, T> && std::is_arithmetic_v<TI>>::type* = nullptr>
    const MultiComplex operator-(TI r) const {
        decltype(coef) new_coef = coef; new_coef[0] -= r;
        return new_coef;
    }
    /// Inplace subtraction of a scalar (this -= r)
    template<typename TI, typename std::enable_if<!std::is_same_v<TI, T> && std::is_arithmetic_v<TI>>::type* = nullptr>
    auto& operator-=(const TI& r) {
        coef[0] -= r;
        return *this;
    }
    /// Inplace subtraction of a scalar (this -= r)
    auto& operator-=(const T& r) {
        coef[0] -= r;
        return *this;
    }
    auto& operator-=(const MultiComplex& w) {
        auto n = *this - w;
        std::swap(*this, n);
        return *this;
    }
    /// Right multiplication by a scalar (this*r)
    template<typename TI, typename std::enable_if<!std::is_same_v<TI, T> && std::is_arithmetic_v<TI>>::type* = nullptr>
    const MultiComplex operator*(TI r) const {
        decltype(coef) new_coef = coef*r;
        return new_coef;
    }
    /// Inplace multiplication by a scalar (this *= r)
    template<typename TI, typename std::enable_if<!std::is_same_v<TI, T> && std::is_arithmetic_v<TI>>::type* = nullptr>
    auto& operator*=(const TI& r) {
        coef *= r;
        return *this;
    }
    /// Inplace multiplication by a scalar (this *= r)
    auto& operator*=(const T& r) {
        coef *= r;
        return *this;
    }
    auto& operator*=(const MultiComplex& w) {
        auto n = *this * w;
        std::swap(*this, n);
        return *this;
    }
    /// Right division by a scalar (this/r)
    template<typename TI, typename std::enable_if<!std::is_same_v<TI, T> && std::is_arithmetic_v<TI>>::type* = nullptr>
    const MultiComplex operator/(TI r) const {
        decltype(coef) new_coef = coef/r;
        return new_coef;
    }
    /// Inplace division by a scalar (this /= r)
    template<typename TI, typename std::enable_if<!std::is_same_v<TI, T> && std::is_arithmetic_v<TI>>::type* = nullptr>
    auto& operator/=(const TI& r) {
        coef /= r;
        return *this;
    }
    /// Inplace division by a scalar (this /= r)
    auto& operator/=(const T& r) {
        coef /= r;
        return *this;
    }
    auto& operator/=(const MultiComplex& w) {
        auto n = *this / w;
        std::swap(*this, n);
        return *this;
    }
    /// Right addition by another MultiComplex
    const MultiComplex operator+(const MultiComplex& w) const {
        const MultiComplex& z = *this;
        auto Nw = w.dim(), Nz = z.dim();
        if (Nw == Nz) {
            return MultiComplex<T>(z.get_coef() + w.get_coef());
        }
        else {
            const MultiComplex& big = (Nw > Nz) ? w : z;
            const MultiComplex& small = (Nw < Nz) ? w : z;
            auto slice = std::slice(0, small.get_coef().size(), 1);
            std::valarray<T> new_coef = big.get_coef(); // copy
            new_coef[slice] += small.get_coef()[slice];
            return new_coef;
        }
    }
    /// Right subtraction by another MultiComplex
    const MultiComplex operator-(const MultiComplex& w) const {
        const MultiComplex& z = *this;
        auto Nw = w.dim(), Nz = z.dim();
        const MultiComplex& big = (Nw > Nz) ? w : z;
        const MultiComplex& small = (Nw < Nz) ? w : z;
        auto slice = std::slice(0, small.get_coef().size(), 1);
        if (Nw == Nz) {
            return MultiComplex<T>(z.get_coef() - w.get_coef());
        }
        else if (Nw > Nz) {
            std::valarray<T> new_coef = -big.get_coef(); // copy
            new_coef[slice] += small.get_coef()[slice];
            return new_coef;
        }
        else {
            std::valarray<T> new_coef = big.get_coef(); // copy
            new_coef[slice] -= small.get_coef()[slice];
            return new_coef;
        }
    }

    template<typename VEC>
    void _times(const VEC& z, std::size_t iz, const VEC& w, std::size_t iw, VEC& p, std::size_t ip, std::size_t L, int sign) const{
        if (L == 2) {
            p[ip] += sign * (z[iz] * w[iw] - z[iz + 1] * w[iw + 1]);
            p[ip + 1] += sign * (z[iz] * w[iw + 1] + z[iz + 1] * w[iw]);
            return;
        }
        auto L2 = L/2;
        _times(z, iz,    w, iw,    p, ip,    L2,  sign);
        _times(z, iz+L2, w, iw+L2, p, ip,    L2, -sign);
        _times(z, iz,    w, iw+L2, p, ip+L2, L2,  sign);
        _times(z, iz+L2, w, iw,    p, ip+L2, L2,  sign);
    }

    /// Right multiplication by another MultiComplex
    const MultiComplex operator*(const MultiComplex& w) const {
        const MultiComplex& z = *this;
        if (w.dim() == 1 && z.dim() == 1) {
            auto zw = z.complex() * w.complex();
            std::valarray<T> coef(2); coef[0] = static_cast<T>(zw.real()); coef[1] = static_cast<T>(zw.imag());
            return coef;
        }
        int Lz = static_cast<int>(z.get_coef().size()),
            Lw = static_cast<int>(w.get_coef().size());
        // Allocate a work array, filled with zeros
        std::valarray<T> p(static_cast<T>(0), std::max(Lz, Lw));

        if (Lz < Lw) {
            for (int i = 0; i < Lw/Lz; ++i) {
                _times(z.get_coef(), 0, w.get_coef(), 0 + i*Lz, p, 0 + i*Lz, Lz, 1);
            }
        }
        else if (Lz > Lw) {
            for (int i = 0; i < Lz/Lw; ++i) {
                _times(w.get_coef(), 0, z.get_coef(), 0 + i*Lw, p, 0 + i*Lw, Lw, 1);
            }
        }
        else {
            _times(z.get_coef(), 0, w.get_coef(), 0, p, 0, Lz, 1);
        }
        return p;
    }

    template<typename TN>
    static auto walloc(std::size_t L) {
        return std::make_tuple(
            std::valarray<TN>(0.0, L),
            std::valarray<TN>(0.0, L),
            std::valarray<TN>(0.0, 2 * L));
    }

    template<typename VEC>
    void _cossin(
        const VEC& d, std::size_t id, VEC& c, std::size_t ic, VEC& s, std::size_t is,
        VEC& w1, std::size_t i1, VEC& w2, std::size_t i2,
        VEC& w3, std::size_t i3, VEC& w4, std::size_t i4,
        std::size_t L) const
    {
        if (L == 2) {
            // In order to avoid duplicate calculations (maybe compiler would 
            // optimize this away, but not clear), pre-calculate the values
            T coshdp1 = scalar_cosh(d[id + 1]), sinhdp1 = scalar_sinh(d[id + 1]);
            T cosd = scalar_cos(d[id]), sind = scalar_sin(d[id]);
            c[ic] = cosd * coshdp1;
            c[ic + 1] = -sind * sinhdp1;
            s[is] = sind * coshdp1;
            s[is + 1] = cosd * sinhdp1;
            return;
        }
        auto L2 = L / 2;
        _cossin(d, id, w1, i1, w2, i2, c, ic, c, ic + L2, s, is, s, is + L2, L2);
        _coshsinh(d, id + L2, w3, i3, w4, i4, c, ic, c, ic + L2, s, is, s, is + L2, L2);
        c[myslice(ic, ic + L, 1)] = 0; s[myslice(is, is + L, 1)] = 0;
        _times(w1, i1, w3, i3, c, ic, L2, 1);
        _times(w2, i2, w4, i4, c, ic + L2, L2, -1);
        _times(w2, i2, w3, i3, s, is, L2, 1);
        _times(w1, i1, w4, i4, s, is + L2, L2, 1);
    }

    template<typename VEC>
    void _coshsinh(
        const VEC& d, std::size_t id,
        VEC& c, std::size_t ic, VEC& s, std::size_t is,
        VEC& w1, std::size_t i1, VEC& w2, std::size_t i2,
        VEC& w3, std::size_t i3, VEC& w4, std::size_t i4,
        std::size_t L) const
    {
        if (L == 2) {
            // In order to avoid duplicate calculations (maybe compiler would 
            // optimize this away, but not clear), pre-calculate the values
            T cosdp1 = scalar_cos(d[id + 1]), sindp1 = scalar_sin(d[id + 1]);
            T coshd = scalar_cosh(d[id]), sinhd = scalar_sinh(d[id]);
            c[ic] = coshd * cosdp1;
            c[ic + 1] = sinhd * sindp1;
            s[is] = sinhd * cosdp1;
            s[is + 1] = coshd * sindp1;
            return;
        }
        std::size_t L2 = L / 2;
        _cossin(d, id + L2, w1, i1, w2, i2, c, ic, c, ic + L2, s, is, s, is + L2, L2);
        _coshsinh(d, id, w3, i3, w4, i4, c, ic, c, ic + L2, s, is, s, is + L2, L2);
        c[myslice(ic, ic + L, 1U)] = 0; s[myslice(is, is + L, 1U)] = 0;
        _times(w1, i1, w3, i3, c, ic, L2, 1);
        _times(w2, i2, w4, i4, c, ic + L2, L2, 1);
        _times(w1, i1, w4, i4, s, is, L2, 1);
        _times(w2, i2, w3, i3, s, is + L2, L2, 1);
    }

    template<typename VEC>
    auto cossin(const VEC& d) const {
        VEC w1, w2, w3;
        std::size_t L = d.size(), L2 = L / 2;
        std::tie(w1, w2, w3) = walloc<T>(L);
        _cossin(d, 0, w1, 0, w2, 0, w3, 0, w3, 0 + L2, w3, 0 + 2 * L2, w3, 0 + 3 * L2, L);
        return std::make_tuple(w1, w2);
    }

    template<typename VEC>
    auto coshsinh(const VEC& d) const {
        VEC w1, w2, w3;
        std::size_t L = d.size(), L2 = L / 2;
        std::tie(w1, w2, w3) = walloc<T>(L);
        _coshsinh(d, 0, w1, 0, w2, 0, w3, 0, w3, 0 + L2, w3, 0 + 2 * L2, w3, 0 + 3 * L2, L);
        return std::make_tuple(w1, w2);
    }

    template<typename VEC>
    void _exp(
        const VEC& d, std::size_t id,
        VEC& e, std::size_t ie, 
        VEC& w1, std::size_t i1,
        std::size_t L) const
    {
        if (L == 2) {
            T expd = scalar_exp(d[id]);
            e[ie] = expd * scalar_cos(d[id + 1]);
            e[ie + 1] = expd * scalar_sin(d[id + 1]);
            return;
        }
        std::size_t L2 = L / 2;
        _cossin(d, id + L2, w1, i1, w1, i1 + L2, w1, i1 + 2 * L2, w1, i1 + 3 * L2, e, ie, e, ie + L2, L2);
        _exp(d, id, w1, i1 + 2 * L2, e, ie, L2);
        e[myslice(ie, ie + L, 1U)] = 0; 
        _times(w1, i1 + 2* L2, w1, i1, e, ie, L2, 1);
        _times(w1, i1 + 2* L2, w1, i1 + L2, e, ie + L2, L2, 1);
    }

    MultiComplex exp() const {
        if (dim() == 1){
            auto z = scalar_exp(complex());
            std::valarray<T> coef(2); coef[0] = static_cast<T>(z.real()); coef[1] = static_cast<T>(z.imag());
            return coef;
        }
        std::size_t L = coef.size(), L2 = L / 2; 
        std::valarray<T> w1 = std::valarray<T>(0.0, L), w2=std::valarray<T>(0.0, 2*L);
        _exp(coef, 0, w1, 0, w2, 0, L);
        return w1;
    }
    T exp(T x) const { return scalar_exp(x); }

    MultiComplex cos() const {
        if (dim() == 1) {
            return std::cos(complex());
        }
        else {
            return std::get<0>(cossin(get_coef()));
        }
    }
    T cos(T x) const { return scalar_cos(x); }
    
    
    MultiComplex sin() const {
        if (dim() == 1) {
            return std::sin(complex());
        }
        else {
            return std::get<1>(cossin(get_coef()));
        }
    }
    T sin(T x) const { return scalar_sin(x); }
    

    MultiComplex cosh() const {
        if (dim() == 1) {
            return std::cosh(complex());
        }
        else {
            return std::get<0>(coshsinh(get_coef()));
        }
    }
    T cosh(T x) const { return scalar_cosh(x); }

    MultiComplex sinh() const {
        if (dim() == 1) {
            return std::sinh(complex());
        }
        else {
            return std::get<1>(coshsinh(get_coef()));
        }
    }
    T sinh(T x) const { return scalar_sinh(x); }

    MultiComplex conj() const {
        std::valarray<T> new_coef = coef;
        for (auto i = coef.size() / 2; i < coef.size(); ++i) {
            new_coef[i] *= -1;
        }
        return new_coef;
    }

    MultiComplex squared() const {
        return (*this)*(*this);
    }

    MultiComplex pow(int n) const {
        if (dim()==1){
            auto z = scalar_pow(complex(), n);
            std::valarray<T> coef(2); coef[0] = static_cast<T>(z.real()); coef[1] = static_cast<T>(z.imag());
            return coef;
        }
        int absn = std::abs(n);
        if (absn == 0) {
            return 1.0;
        }
        else if (absn == 1) {
            return (absn == n) ? *this : 1.0/(*this);
        }
        else if (absn == 2) {
            return (absn == n) ? squared() : 1.0/squared();
        }
        else {
            MultiComplex z = *this;
            for (auto i = 0; i < absn-1; ++i) {
                z = z*(*this);
            }
            return (absn == n) ? z : 1.0/z;
        }
    }

    /** Raise to an exponent
    * If real-most part of number is negative, power must be an integer
    */
    MultiComplex pow(double exponent) const {
        int integer_exponent = static_cast<int>(exponent);
        if (dim()==1){
            auto c = complex();
            if (c.real() < 0){
                // Special case 'negative' arguments
                if (exponent == integer_exponent){
                    // exponent (as double) is exact integer representation
                    return pow(integer_exponent);
                }
                else {
                    throw std::range_error("Cannot use 'negative' complex numbers with non-integer exponents");
                }
            }
            else {
                return (exponent * log(*this)).exp();
            }
        }
        else{
            auto c0 = coef[0];
            if (c0 >= 0){
                // a^b = exp(ln(a^b)) = exp(b*ln(a))
                return (exponent*log(*this)).exp();
            }
            else {
                // Special case 'negative' arguments
                if (exponent == integer_exponent) {
                    // exponent (as double) is exact integer representation
                    return pow(integer_exponent);
                }
                else {
                    throw std::range_error("Cannot use 'negative' complex numbers with non-integer exponents");
                }
            }
        }
    }

    /// Multiplicative inverse: multinv(z) = 1/z
    MultiComplex multinv() const {
        int L = static_cast<int>(coef.size()), L2 = L/2;
        if (L2 == 1) {
            return 1/(coef[0]*coef[0] + coef[1]*coef[1])*this->conj();
        }
        else {
            auto left = MultiComplex(coef[myslice(0, L2, 1)]);
            auto right = MultiComplex(coef[myslice(L2, L, 1)]);
            return (left.squared() + right.squared()).multinv()*this->conj();
        }
    }

    /// Right division by another MultiComplex
    const MultiComplex operator/(const MultiComplex& w) const {
        return (*this)*w.multinv();
    }

    /// Convert a set of indices, where each is a 0 or 1, to an index into coef
    /// 
    /// We start off with this object, and walk the pointers to the left
    /// and right depending on the index.  A 0 takes us to the left, a 
    /// 1 to the right
    /// 
    template<typename VEC>
    int branches_to_index(const VEC& inds) const {
        if (inds.size() != m_d) {
            throw std::invalid_argument("Length of index vector of "
                + std::to_string(inds.size())
                + " is not equal to the dimension of"
                + std::to_string(m_d));
        }
        int o = 0, m = 1;
        for (int i = static_cast<int>(inds.size()) - 1; i >= 0; --i) {
            o += inds[i] * m;
            m *= 2;
        }
        return o;
    }

    /// Index into the MultiComplex object by an array of indices, in which each
    /// index corresponds to 0 (take the real branch) or 1 (take the complex branch)
    /// The most real value corresponds to index 0
    template<typename VEC>
    T operator[](const VEC& inds) const {
        return coef[branches_to_index(inds)];
    }

    /// Index into the MultiComplex object by single index
    T operator[](int ind) const {
        return coef[ind];
    }

    /// Set a coefficient by numerical index.  This setter function was written
    /// because the operator[] overloads are const
    void set_coef(int ind, T val) {
        coef[ind] = val;
    }

    friend bool operator< (const MultiComplex<T>& lhs, const MultiComplex<T>& rhs) { 
        return lhs.real() < rhs.real(); 
    };
};

// A type trait to define acceptable numerical types in the prefix operator functions
template<typename T> struct is_acceptable_number : public std::false_type {};
template<> struct is_acceptable_number<int> : public std::true_type {};
template<> struct is_acceptable_number<double> : public std::true_type {};
template<> struct is_acceptable_number<std::complex<double>> : public std::true_type {};

#if defined(BOOST_MULTIPRECISION_FOUND)
template<typename T> struct is_acceptable_number<boost::multiprecision::number<T>> : public std::true_type {};
#endif

/// See https://stackoverflow.com/a/14294277/1360263 for description of type limiting
/// Helper function that allows for pre-addition by calling the postfix function for numerical types
template <
    typename TN, typename T, typename = typename std::enable_if<is_acceptable_number<T>::value, T>::type>
MultiComplex<TN> operator+(T value, const MultiComplex<TN>& mc) {
    return mc + value;
};
/// Helper function that allows for pre-multiplication by calling the postfix function for numerical types
template <
    typename TN, typename T, typename = typename std::enable_if<is_acceptable_number<T>::value, T>::type>
MultiComplex<TN> operator*(T value, const MultiComplex<TN>& mc) {
    return mc * value;
};
/// Helper function that allows for pre-subtraction by calling the postfix function for numerical types
template <
    typename TN, typename T, typename = typename std::enable_if<is_acceptable_number<T>::value, T>::type>
const MultiComplex<TN> operator-(const T value, const MultiComplex<TN>& mc) {
    return -(mc - value);
};
/// Helper function that allows for pre-division by calling the postfix function for numerical types
template <
    typename TN, typename T,typename = typename std::enable_if<is_acceptable_number<T>::value, T>::type>
const MultiComplex<TN> operator/(const T value, const MultiComplex<TN>& mc) {
    return value * mc.multinv();
};

template<typename TN>
MultiComplex<TN> cos(const MultiComplex<TN>& z) {
    return z.cos();
}
template<typename TN>
MultiComplex<TN> sin(const MultiComplex<TN>& z) {
    return z.sin();
}
template<typename TN>
MultiComplex<TN> cosh(const MultiComplex<TN>& z) {
    return z.cosh();
}
template<typename TN>
MultiComplex<TN> sinh(const MultiComplex<TN>& z) {
    return z.sinh();
}
template<typename TN>
MultiComplex<TN> exp(const MultiComplex<TN>& z) {
    return z.exp();
}

template<typename TN>
MultiComplex<TN> log(const MultiComplex<TN>& z) {
    // Normally we would like to short-circuit for normal complex number, but can't here
    // because of the possible branch-cuts
    TN re = z[0], re_old = z[0];
    MultiComplex<TN> y = scalar_log(re); // Start off with the most-real component
    MultiComplex<TN> expny = exp(-y);
    for (auto counter = 0; counter <= 6; ++counter) {
        y = y - 2.0 * (1.0 - z * expny) / (1.0 + z * expny);
        expny = exp(-y);
        re = y[0];
        TN diff = re - re_old; re_old = re;
        if (scalar_abs(diff) < 1e-15 && counter > 0){
            break;
        }
    }
    return y;
}

template<typename TN>
MultiComplex<TN> pow(const MultiComplex<TN>& z, double e) {
    return z.pow(e);
}

template<typename TN>
MultiComplex<TN> pow(const MultiComplex<TN>& z, int e) {
    return z.pow(e);
}

template<typename TN>
MultiComplex<TN> sqrt(const MultiComplex<TN>& z) {
    return z.pow(0.5);
}

/**
 The derivatives of order 1 to n (inclusive) of a function that takes a single variable

 @param f The function to be called. Takes a MultiComplex, returns a MultiComplex.  Consider writing a lambda function that does the calculation
 @param x The real value at which the multicomplex should be instantiated
 @param numderiv The maximum number of derivatives to take.  The larger this number the more computational effort required
 @parma and_val Also return the function value in slot 0
 @return out The vector of numerical derivatives that were obtained.  The first entry is the first derivative, the second is the second derivative, and so on
 */
template<typename TN>
std::vector<TN> diff_mcx1(const std::function<MultiComplex<TN>(const MultiComplex<TN>&)>& f,
    TN x, int numderiv, bool and_val = false)
{
    // The tiny step
    TN DELTA = increment(numderiv);
    // Coeffs of the multicomplex number, filled by default
    // with zeros.  The array of coefficients is of length 2^(numderiv)
    std::valarray<TN> c(0.0, exp2i(numderiv));
    // The real component as passed to function
    c[0] = x;
    // The very, very small offset number goes in all the indices
    // 2^k for 0 <= k < numderiv. For a "normal" complex number,
    // this would be the second entry (index=1)
    for (auto k = 0; k < numderiv; ++k) {
        c[exp2i(k)] = DELTA;
    }
    // Call the function with our multicomplex number
    // (implicit instantiation of MC argument to f)
    auto o = f(c);
    //for(auto i:o.get_coef()){ std::cout << i << std::endl;}
    // Store all the derivatives that were calculated
    std::vector<TN> ders;
    if (and_val) {
        ders.push_back(o[0]);
    }
    for (auto L = 1; L <= numderiv; ++L) {
        // The calculated value is sitting in the 2^L-1 index,
        // and then need to divide by DELTA^L
        ders.push_back(o[int(exp2i(L) - 1)] / scalar_pow(DELTA, L));
    }
    return ders;
}

/**
 The derivatives of order 1 to n (inclusive) of a function that takes a single variable

 @param f The function to be called. Takes a MultiComplex, returns a tuple of MultiComplex values.  Consider writing a lambda function that does the calculation
 @param x The real value at which the multicomplex should be instantiated
 @param numderiv The maximum number of derivatives to take.  The larger this number the more computational effort required
 @parma and_val Also return the function value in slot 0
 @return out The vector of numerical derivatives that were obtained.  The first entry is the first derivative, the second is the second derivative, and so on

 This function was written in particular to take derivatives of integrands for which the value is associated with an error estimate
 */
template<typename TN>
std::tuple<std::vector<TN>, std::vector<TN>> diff_mcx1(
    const std::function<std::tuple<MultiComplex<TN>, MultiComplex<TN>>(const MultiComplex<TN>&)>& f,
    TN x, int numderiv, bool and_val = false)
{
    // The tiny step
    TN DELTA = increment(numderiv);
    // Coeffs of the multicomplex number, filled by default
    // with zeros.  The array of coefficients is of length 2^(numderiv)
    std::valarray<TN> c(0.0, exp2i(numderiv));
    // The real component as passed to function
    c[0] = x;
    // The very, very small offset number goes in all the indices
    // 2^k for 0 <= k < numderiv. For a "normal" complex number,
    // this would be the second entry (index=1)
    for (auto k = 0; k < numderiv; ++k) {
        c[exp2i(k)] = DELTA;
    }
    // Call the function with our multicomplex number
    // (implicit instantiation of MC argument to f)
    auto [o,e] = f(c);
    //for(auto i:o.get_coef()){ std::cout << i << std::endl;}
    // Store all the derivatives that were calculated
    std::vector<TN> ders,errs;
    if (and_val) {
        ders.push_back(o[0]);
        errs.push_back(e[0]);
    }
    for (auto L = 1; L <= numderiv; ++L) {
        // The calculated value is sitting in the 2^L-1 index,
        // and then need to divide by DELTA^L
        ders.push_back(o[int(exp2i(L) - 1)] / scalar_pow(DELTA, L));
        errs.push_back(e[int(exp2i(L) - 1)] / scalar_pow(DELTA, L));
    }
    return std::make_tuple(ders,errs);
}

// Base template for function traits
template<typename T> struct function_traits;

// Partial specialization for std::function argument
template<typename R, typename ...Args>
struct function_traits<std::function<R(Args...)>>
{
    using argtype = typename std::decay<typename std::tuple_element<0, std::tuple<Args...>>::type>::type;
};

/**
 The specified derivatives of a function that takes multiple variables

 @param f The function to be called. Takes a vector of MultiComplex, returns a MultiComplex.  Consider writing a lambda function that does the calculation.
 @param x The vector of real values at which the multicomplex instances should be instantiated
 @param order The order of derivatives w.r.t. the independent variables. [1,0,1] would a first derivative with respect to the first variable and a first partial derivative with respect to the third variable
 @return out The numerical derivative that was obtained.  
 */
template<typename FuncType, typename PointType>
typename PointType::value_type diff_mcxN(
    const FuncType &f,
    const PointType &x,
    const std::vector<int> &orders)
{
    if (x.size() != orders.size()) {
        throw std::invalid_argument(
            "Length of x: "
            + std::to_string(x.size()) 
            + " does not equal that of the order: " 
            + std::to_string(orders.size())
        );
    }
    // The total number of derivatives to take
    int numderiv = std::accumulate(orders.begin(), orders.end(), 0);

    using TN = typename PointType::value_type;
    
    // The tiny step
    TN DELTA = increment(numderiv);
    
    // Get the type of the container of arguments to function, to allow for std::vector, std::array, etc.
    constexpr std::size_t zero = 0;
    using MCVecType = typename function_traits<FuncType>::argtype;
    MCVecType zs(x.size());

    int k_counter = 0;
    for (auto i = 0; i < x.size(); ++i){
        // Coeffs of the multicomplex number, filled by default
        // with zeros.  The array of coefficients is of length 2^(numderiv)
        std::valarray<TN> c(0.0, exp2i(numderiv));
        // The real component as passed to function
        c[0] = x[i];
        // The very, very small offset number goes in all the indices
        // 2^k for 0 <= k < order[i]. 
        for (int k = 0; k < orders[i]; ++k) {
            c[exp2i(k_counter)] = DELTA;
            k_counter++;
        }
        zs[i] = c;
    }
    // Call the function with our multicomplex arguments
    auto o = f(zs);
    // Return the desired derivative
    return o[exp2i(numderiv)-1]/scalar_pow(DELTA, numderiv);
}

namespace detail{
    // Generic setting functions to handle Eigen types and STL types with the same interface
    template<typename MatrixLike, typename Integer, typename ValType>
    void setval(MatrixLike& m, Integer i, Integer j, const ValType val) {
        m(i, j) = val;
    }

    // Partial specialization for valarray "matrix"
    template <> inline void setval<std::valarray<std::valarray<double>>, std::size_t, double>(std::valarray<std::valarray<double>>& m, std::size_t i, std::size_t j, const double val) {
        m[i][j] = val;
    }

    // Generic resizing functions to handle Eigen types and STL types with the same interface
    template<typename MatrixLike, typename Integer>
    void resizemat(MatrixLike& m, Integer i, Integer j) {
        m.resize(i, j);
    }

    // Partial specialization for valarray "matrix"
    template <> inline void resizemat<std::valarray<std::valarray<double>>, std::size_t>(std::valarray<std::valarray<double>>& m, std::size_t M, std::size_t N) {
        m.resize(M);
        for (auto i = 0; i < M; ++i){
            m[i] = std::valarray<double>(0.0, N);
        }
    }
}

enum class HessianMethods {OneBig, Multiple};
template<typename MatType, typename FuncType, typename ArgType, HessianMethods method = HessianMethods::OneBig>
auto get_Hessian(const FuncType &f, const ArgType& x)
{
    if (x.size() != 2) {
        throw std::invalid_argument(
            "Length of x: "
            + std::to_string(x.size())
            + " does not equal 2"
        );
    }
    MatType H;
    std::size_t N = 2;
    detail::resizemat(H, N, N);

    if constexpr (method == HessianMethods::OneBig) {
        throw std::invalid_argument("Not yet implemented; also probably a bad idea anyway");
        //// Take a single derivative w.r.t. independent variables

        //// The total number of derivatives to take
        //std::vector<int> orders = {3, 3};
        //int numderiv = 6; // 3 in each variable so we get all the second partial derivatives

        //using TN = ArgType::value_type;

        //// The tiny step
        //auto DELTA = increment(numderiv);

        //using MCVecType = typename function_traits<FuncType>::arg<0>::type;
        //MCVecType zs(x.size());
        //int k_counter = 0;
        //for (auto i = 0; i < x.size(); ++i) {
        //    // Coeffs of the multicomplex number, filled by default
        //    // with zeros. The array of coefficients is of length 2^(numderiv)
        //    std::valarray<TN> c(0.0, exp2i(numderiv));
        //    // The real component as passed to function
        //    c[0] = x[i];
        //    // The very, very small offset number goes in all the indices
        //    // 2^k for 0 <= k < order[i]. 
        //    for (int k = 0; k < orders[i]; ++k) {
        //        c[exp2i(k_counter)] = DELTA;
        //        k_counter++;
        //    }
        //    zs[i] = c;
        //}
        //// Call the function with our multicomplex arguments
        //auto o = f(zs);
        //// Extract the desired derivative
        //for (auto i = 0; i < 2; ++i) {
        //    for (auto j = 0; j < 2; ++j) {
        //        numderiv = i+j;
        //        H[i][j] = o[exp2i(numderiv) - 1] / pow(DELTA, numderiv);
        //    }
        //}
        //return H;
    }
    else if constexpr (method == HessianMethods::Multiple) {
        // Take multiple second derivatives, for each commbination
        for (std::size_t i = 0; i < x.size(); ++i) {
            for (std::size_t j = i; j < x.size(); ++j) {
                std::vector<int> order = { 0, 0 };
                order[i] += 1;
                order[j] += 1;
                auto val = diff_mcxN(f, x, order);
                detail::setval(H, i, j, val);
                detail::setval(H, j, i, val);
            }
        }
        return H;
    }
    else {
        std::invalid_argument("incomplete options");
    }
}

}; // namespace mcx
#endif
