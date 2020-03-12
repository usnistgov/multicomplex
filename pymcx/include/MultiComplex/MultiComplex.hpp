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

// A custom slice, more like Python, where the arguments to the slice are start, stop, increment
// rather than start, number of times, increment
auto myslice(std::size_t start, std::size_t stop, std::size_t increment) {
    return std::slice(start, (stop - start) / increment, increment);
}

/// Raise 2 to the power i, where i is an integer
// see https://stackoverflow.com/a/5345369 for explanation
int exp2i(int i) {
    return 1 << i;
}

/// Log2i, where i is an integer
/// Could still be more optimized if needed
int log2i(int i) {
    int val = int(log2(i));
    if (exp2i(val) != i) {
        throw std::invalid_argument("The argument to log2i " + std::to_string(i) + " is not a power of 2");
    }
    return val;
}

inline double increment(std::size_t l) {
    return exp2(-664 / static_cast<int>(l));
}

template<typename T >
struct MultiComplex
{
    int m_d = 0;
    std::valarray<T> coef;

    /// Default constructor
    MultiComplex() : m_d(0) { coef = {}; };

    /// This is a "normal" real number (re+i*0), stored as a complex number
    MultiComplex(const T re) : m_d(1) { coef = { re, 0 }; };

    // This is a "normal" complex number, stored as such 
    MultiComplex(const std::complex<T>& c) : m_d(1) { coef = { c.real(), c.imag() }; };

    /// This is a multicomplex complex number (re+i1*im), with two children,
    /// the left child is the "real" part (itself a multicomplex number), 
    /// and the right child is the "imaginary" part (itself a multicomplex number)
    //MultiComplex(const MultiComplex<T>& re, const MultiComplex<T>& im, int d = 2) : m_d(d) { children.push_back(re); children.push_back(im);  };

    //MultiComplex(const MultiComplex<T>&& re, const MultiComplex<T>&& im, int d = 2) : m_d(d) { children.emplace_back(std::move(re)); children.emplace_back(std::move(im)); move++;  };

    /// Constructor in which a vector of values are provided
    /// These are the values of the terminals, and the length of the values will be 2^d, 
    /// where d is the dimension (or level) of the MultiComplex object
    MultiComplex(const std::valarray<T>& v) {
        coef = v;
        m_d = log2i(static_cast<int>(v.size()));
    }
    MultiComplex(const std::valarray<T>&& v) {
        coef = std::move(v);
        m_d = log2i(static_cast<int>(v.size()));
    }

    /// The dimension of the multicomplex number
    int dim() const { return m_d; }

    /// The "most real" component; infinity if not initialized yet
    T real() const {
        return (m_d > 0) ? coef[0] : std::numeric_limits<double>::infinity();
    }

    std::complex<T> complex() const {
        return { coef[0], coef[1] };
    }

    /// Get a const reference to the coefficients
    const std::valarray<T>& get_coef() const { return coef; }

    /// Unary negation parameter (-this)
    const MultiComplex operator-() const {
        return MultiComplex<T>(-coef);
    }
    /// Right addition by a float (this+r)
    const MultiComplex operator+(T r) const {
        std::valarray<double> new_coef = coef; new_coef[0] += r;
        return new_coef;
    }
    /// Right subtraction by a float (this-r)
    const MultiComplex operator-(T r) const {
        std::valarray<double> new_coef = coef; new_coef[0] -= r;
        return new_coef;
    }
    /// Right multiplication by a float (this*r)
    const MultiComplex operator*(T r) const {
        return MultiComplex<T>(coef*r);
    }
    /// Right division by a float (this*r)
    const MultiComplex operator/(T r) const {
        return MultiComplex<T>(coef/r);
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
    static void _times(const VEC& z, std::size_t iz, const VEC& w, std::size_t iw, VEC& p, std::size_t ip, std::size_t L, int sign) {
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
            return z.complex() * w.complex();
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
    static void _cossin(
        const VEC& d, std::size_t id, VEC& c, std::size_t ic, VEC& s, std::size_t is,
        VEC& w1, std::size_t i1, VEC& w2, std::size_t i2,
        VEC& w3, std::size_t i3, VEC& w4, std::size_t i4,
        std::size_t L)
    {
        if (L == 2) {
            // In order to avoid duplicate calculations (maybe compiler would 
            // optimize this away, but not clear), pre-calculate the values
            T coshdp1 = cosh(d[id + 1]), sinhdp1 = sinh(d[id + 1]);
            T cosd = cos(d[id]), sind = sin(d[id]);
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
    static void _coshsinh(
        const VEC& d, std::size_t id,
        VEC& c, std::size_t ic, VEC& s, std::size_t is,
        VEC& w1, std::size_t i1, VEC& w2, std::size_t i2,
        VEC& w3, std::size_t i3, VEC& w4, std::size_t i4,
        std::size_t L)
    {
        if (L == 2) {
            // In order to avoid duplicate calculations (maybe compiler would 
            // optimize this away, but not clear), pre-calculate the values
            T cosdp1 = cos(d[id + 1]), sindp1 = sin(d[id + 1]);
            T coshd = cosh(d[id]), sinhd = sinh(d[id]);
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
    static auto cossin(const VEC& d) {
        VEC w1, w2, w3;
        std::size_t L = d.size(), L2 = L / 2;
        std::tie(w1, w2, w3) = walloc<double>(L);
        _cossin(d, 0, w1, 0, w2, 0, w3, 0, w3, 0 + L2, w3, 0 + 2 * L2, w3, 0 + 3 * L2, L);
        return std::make_tuple(w1, w2);
    }

    template<typename VEC>
    static auto coshsinh(const VEC& d) {
        VEC w1, w2, w3;
        std::size_t L = d.size(), L2 = L / 2;
        std::tie(w1, w2, w3) = walloc<double>(L);
        _coshsinh(d, 0, w1, 0, w2, 0, w3, 0, w3, 0 + L2, w3, 0 + 2 * L2, w3, 0 + 3 * L2, L);
        return std::make_tuple(w1, w2);
    }

    template<typename VEC>
    static void _exp(
        const VEC& d, std::size_t id,
        VEC& e, std::size_t ie, 
        VEC& w1, std::size_t i1,
        std::size_t L)
    {
        if (L == 2) {
            T expd = std::exp(d[id]);
            e[ie] = expd * cos(d[id + 1]);
            e[ie + 1] = expd * sin(d[id + 1]);
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
            return std::exp(complex());
        }
        std::size_t L = coef.size(), L2 = L / 2; 
        std::valarray<T> w1 = std::valarray<T>(0.0, L), w2=std::valarray<T>(0.0, 2*L);
        _exp(coef, 0, w1, 0, w2, 0, L);
        return w1;
    }

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
            return std::pow(complex(), n);
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

    MultiComplex pow(double exponent) const {
        if (dim()==1){
            return std::pow(complex(), exponent);
        }
        else{
            return (exponent*log(*this)).exp();
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

/// Helper function that allows for pre-addition by calling the postfix function
template <typename TN>
MultiComplex<TN> operator+(TN value, const MultiComplex<TN>& mc) {
    return mc + value;
};
/// Helper function that allows for pre-multiplication by calling the postfix function
template <typename TN>
MultiComplex<TN> operator*(TN value, const MultiComplex<TN>& mc) {
    return mc * value;
};
/// Helper function that allows for pre-subtraction by calling the postfix function
template <typename TN>
const MultiComplex<TN> operator-(const TN value, const MultiComplex<TN>& mc) {
    return -(mc - value);
};
/// Helper function that allows for pre-division by calling the postfix function
template <typename TN>
const MultiComplex<TN> operator/(const TN value, const MultiComplex<TN>& mc) {
    return value * mc.multinv();
};

template<typename TN>
MultiComplex<TN> cos(const MultiComplex<TN>& z) {
    return (z.dim() == 1) ? std::cos(z.complex()) : MultiComplex<TN>(std::get<0>(MultiComplex<TN>::cossin(z.get_coef())));
}
template<typename TN>
MultiComplex<TN> sin(const MultiComplex<TN>& z) {
    return (z.dim() == 1) ? std::sin(z.complex()) : MultiComplex<TN>(std::get<1>(MultiComplex<TN>::cossin(z.get_coef())));
}
template<typename TN>
MultiComplex<TN> cosh(const MultiComplex<TN>& z) {
    return (z.dim() == 1) ? std::cosh(z.complex()) : MultiComplex<TN>(std::get<0>(MultiComplex<TN>::coshsinh(z.get_coef())));
}
template<typename TN>
MultiComplex<TN> sinh(const MultiComplex<TN>& z) {
    return (z.dim() == 1) ? std::sinh(z.complex()) : MultiComplex<TN>(std::get<1>(MultiComplex<TN>::coshsinh(z.get_coef())));
}
template<typename TN>
MultiComplex<TN> exp(const MultiComplex<TN>& z) {
    if (z.dim() == 1) {
        return std::exp(z.complex());
    }
    else {
        return z.exp();
    }
}

template<typename TN>
MultiComplex<TN> log(const MultiComplex<TN>& z) {
    // Normally we would like to short-circuit for normal complex number, but can't here
    // because of the possible branch-cuts
    TN re = z[0], re_old = z[0];
    MultiComplex<TN> y = log(z[0]); // Start off with the most-real component
    MultiComplex<TN> expny = exp(-y);
    for (auto counter = 0; counter <= 6; ++counter) {
        y = y - 2.0 * (1.0 - z * expny) / (1.0 + z * expny);
        expny = exp(-y);
        re = y[0];
        TN diff = re - re_old; re_old = re;
        if (std::abs(diff) < 1e-15 && counter > 0){
            break;
        }
    }
    return y;
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
        ders.push_back(o[int(exp2i(L) - 1)] / pow(DELTA, L));
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
        ders.push_back(o[int(exp2i(L) - 1)] / pow(DELTA, L));
        errs.push_back(e[int(exp2i(L) - 1)] / pow(DELTA, L));
    }
    return std::make_tuple(ders,errs);
}

/**
 The specified derivatives of a function that takes multiple variables

 @param f The function to be called. Takes a vector of MultiComplex, returns a MultiComplex.  Consider writing a lambda function that does the calculation.
 @param x The vector of real values at which the multicomplex instances should be instantiated
 @param order The order of derivatives w.r.t. the independent variables. [1,0,1] would a first derivative with respect to the first variable and a first partial derivative with respect to the third variable
 @return out The numerical derivative that was obtained.  
 */
template<typename TN>
TN diff_mcxN(
    std::function<MultiComplex<TN>(const std::vector<MultiComplex<TN>> &)> f,
    const std::vector<TN> &x, 
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
    
    // The tiny step
    TN DELTA = increment(numderiv);
    
    std::vector<MultiComplex<TN>> zs;
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
        zs.emplace_back(c);
    }
    // Call the function with our vector of multicomplex arguments
    auto o = f(zs);
    // Return the desired derivative
    return o[exp2i(numderiv)-1]/pow(DELTA, numderiv);
}

#endif
