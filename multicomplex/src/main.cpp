#include <chrono>
#include <vector>
#include <iostream>
#include <complex>
#include <atomic>

template<typename T >
struct MultiComplex
{
    int m_d = 0;
    std::complex<T> m_complex;
    std::vector<MultiComplex> children;

    /// This is a "normal" complex number (re+i*im), stored as such
    MultiComplex(T re, T im) : m_d(1) { m_complex = { re, im }; };

    /// This is a "normal" complex number, stored as such 
    MultiComplex(const std::complex<T>& c) : m_d(1) { m_complex = c; };

    /// This is a multicomplex complex number (re+i1*im), with two children,
    /// the left child is the "real" part (itself a multicomplex number), 
    /// and the right child is the "imaginary" part (itself a multicomplex number)
    MultiComplex(const MultiComplex<T>& re, const MultiComplex<T>& im, int d = 2) : m_d(d) { children.push_back(re); children.push_back(im);  };

    //MultiComplex(const MultiComplex<T>&& re, const MultiComplex<T>&& im, int d = 2) : m_d(d) { children.emplace_back(std::move(re)); children.emplace_back(std::move(im)); move++;  };

    /// Constructor in which a vector of values are provided
    /// These are the values of the terminals, and the length of the values will be 2^d, 
    /// where d is the dimension (or level) of the MultiComplex object
    MultiComplex(const std::vector<double>& v) {
        if (v.size() < 2) {
            throw std::invalid_argument("MultiComplex: length(v)<2 not meaningful");
        }
        else if (v.size() == 2) {
            // Complex number
            m_complex = { v[0], v[1] }; m_d = 1;
        }
        else if (v.size() % 2 != 0) {
            // Though this check is just for length being divisible by 2, it is enough because any recursive 
            // call into this function will also potentially cause this test to fail
            throw std::invalid_argument("MultiComplex: length(v) not power of 2");
        }
        else {
            children.emplace_back(std::vector<double>(v.begin(), v.begin() + v.size()/2));
            children.emplace_back(std::vector<double>(v.begin() + v.size()/2, v.end()));
            m_d = children[0].dim() + 1;
        }
    }

    /// The dimension of the multicomplex number
    int dim() const { return m_d; }

    // Return the "real" leaf in the tree
    // See https://www.fluentcpp.com/2018/07/13/the-incredible-const-reference-that-isnt-const/
    std::remove_reference_t<MultiComplex> const& re() const { return children[0]; }

    // Return the "imaginary" leaf in the tree
    std::remove_reference_t<MultiComplex> const& im() const { return children[1]; }

    /// Unary negation parameter (-this)
    const MultiComplex operator-() const {
        if (m_d == 1) {
            return -m_complex;
        }
        else {
            return MultiComplex(-re(), -im(), m_d);
        }
    }

    /// Right addition by a float (this+r)
    const MultiComplex operator+(T r) const {
        if (m_d == 1) {
            return m_complex + r;
        }
        else {
            return MultiComplex(re()+r, im(), m_d);
        }
    }
    
    /// Right addition by another MultiComplex
    const MultiComplex operator+(const MultiComplex& w) const{
        const MultiComplex& z = *this;
        if (z.dim() == 1 && w.dim() == 1) {
            return z.complex() + w.complex();
        }
        else if (dim() < w.dim()) {
            return MultiComplex(z+w.re(), w.im(), w.dim());
        }
        else if (dim() > w.dim()) {
            return MultiComplex(z.re()+w, z.im(), z.dim());
        }
        else {
            return MultiComplex(z.re()+w.re(), z.im()+w.im(), z.dim());
        }
    }
    
    /// Right multiplication by a float (this*r)
    const MultiComplex operator*(T r) const {
        if (m_d == 1) {
            return m_complex * r;
        }
        else {
            return MultiComplex(re() * r, im() * r, m_d);
        }
    }
    
    /// Right multiplication by another MultiComplex
    const MultiComplex operator*(const MultiComplex &w) const {
        const MultiComplex& z = *this;
        if (z.dim() == 1 && w.dim() == 1) {
            return z.complex() * w.complex();
        }
        else if (z.dim() < w.dim()) {
            return MultiComplex(z*w.re(), z*w.im(), w.dim());
        }
        else if (dim() > w.dim()) {
            return MultiComplex(z.re()*w, z.im()*w, z.dim());
        }
        else {
            return MultiComplex(z.re()*w.re() - (z.im()*w.im()), z.re()*w.im() + z.im()*w.re(), dim());
        }
    }
    
    /// Right subtraction by a float (this-r)
    const MultiComplex operator-(T r) const {
        if (m_d) {
            return m_complex - r;
        }
        else {
            return MultiComplex(re() - r, im(), m_d);
        }
    }
    
    /// Right subtraction by another MultiComplex
    const MultiComplex operator-(const MultiComplex& w) const {
        const MultiComplex& z = *this;
        if (z.dim() == 1 && w.dim() == 1) {
            return z.complex() - w.complex();
        }
        else if (z.dim() < w.dim()) {
            return MultiComplex(z - w.re(), -w.im(), w.dim());
        }
        else if (z.dim() > w.dim()) {
            return MultiComplex(z.re() - w, z.im() - w.im(), z.dim());
        }
        else {
            return MultiComplex(z.re() - w.re(), z.im()-w.im(), dim());
        }
    }
    
    /// Return a std::complex copy of this class
    const std::complex<T> & complex() const {
        switch (m_d) {
        case 1:
            return m_complex;
        default:
            throw -1;
        }
    }

    /// Index into the MultiComplex object
    T operator[](const std::vector<int>& inds) const{
        if (inds.size() != m_d) {
            throw std::invalid_argument("Length of index vector of "
                +std::to_string(inds.size())
                +" is not equal to the dimension of"
                + std::to_string(m_d));
        }
        
        // We start off with this object, and walk the pointers to the left
        // and right depending on the index.  A 0 takes us to the left, a 
        // 1 to the right
        // 
        // 1. There are no copies anywhere, only pointer manipulations
        // 2. The const decasting is needed to satisfy the fact that
        //    this function is const

        // Pointer to const MultiComplex
        MultiComplex const * pMC = const_cast<MultiComplex const *>(this);

        // Iterate over the indices
        for (auto&& ind : inds) {
            // Dereference the pointer for concision (not strictly needed)
            MultiComplex const& mc = *pMC;
            if (ind < 0 || ind > 1) {
                throw std::invalid_argument("An invalid index of "
                    + std::to_string(inds.size())
                    + " was found. Options are 0 or 1");
            }
            if (mc.dim() == 1) {
                // Return the value of the real or complex branch
                return (ind == 0) ? mc.complex().real() : mc.complex().imag();
            }
            else {
                // Take the real branch if 0, or the imaginary branch if 1.
                // Need to cast away the constness of the reference
                pMC = const_cast<MultiComplex const *> (&((ind == 0) ? mc.re() : mc.im()));
            }
        }
        throw std::invalid_argument("Cannot get here, something bad happened");
    }
};

/// Helper function that allow for pre-addition,subtraction,multiplication (e.g., value+mc) by callling
/// the post-addition function
template <typename TN>
MultiComplex<TN> operator+(TN value, const MultiComplex<TN>& mc) {
    return mc + value;
};
template <typename TN>
MultiComplex<TN> operator*(TN value, const MultiComplex<TN>& mc) {
    return mc * value;
};
template <typename TN>
const MultiComplex<TN> operator-(const TN value, const MultiComplex<TN>& mc) {
    return -(mc - value);
};

template<typename TN>
MultiComplex<TN> cos(const MultiComplex<TN>& z) {
    if (z.dim() == 1) {
        return std::cos(z.complex());
    }
    else {
        // A general multicomplex number
        auto r = z.re(), i = z.im();
        return MultiComplex<TN>(cos(r)*cosh(i), -sin(r)*sinh(i), z.dim());
    }
}
template<typename TN>
MultiComplex<TN> sin(const MultiComplex<TN>& z) {
    if (z.dim() == 1) {
        return std::sin(z.complex());
    }
    else {
        // A general multicomplex number
        auto r = z.re(), i = z.im();
        return MultiComplex<TN>(sin(r)*cosh(i), cos(r)*sinh(i), z.dim());
    }
}
template<typename TN>
MultiComplex<TN> cosh(const MultiComplex<TN>& z) {
    if (z.dim() == 1) {
        return std::cosh(z.complex());
    }
    else {
        // A general multicomplex number
        auto r = z.re(), i = z.im();
        return MultiComplex<TN>(cosh(r) * cos(i), sinh(r) * sin(i), z.dim());
    }
}
template<typename TN>
MultiComplex<TN> sinh(const MultiComplex<TN>& z) {
    if (z.dim() == 1) {
        return std::sinh(z.complex());
    }
    else {
        // A general multicomplex number
        auto r = z.re(), i = z.im();
        return MultiComplex<TN>(sinh(r) * cos(i), cosh(r) * sin(i), z.dim());
    }
}
template<typename TN>
MultiComplex<TN> exp(const MultiComplex<TN>& z) {
    if (z.dim() == 1) {
        return std::exp(z.complex());
    }
    else {
        return sinh(z) + cosh(z);
    }
}

void time_one(int d) {
    double DELTA = 1e-100, rx = 0.1234;

    std::vector<double> r(int(exp2(d)), 0.0);
    r[0] = rx;
    r[1] = DELTA;

    std::vector<int> indices(d, 0); indices.back() = 1;
    
    double exact = rx * cos(rx) + sin(rx);
    MultiComplex<double> x1(r);
    double errs = 0;
    auto startTime = std::chrono::system_clock::now();
    for (auto ii = 0; ii < 1e6; ++ii) {

        // A real number and its complex step
        errs += (x1 * exp(x1))[indices] / DELTA - exact;
    }

    auto endTime = std::chrono::system_clock::now();
    double elap = std::chrono::duration<double>(endTime - startTime).count();
    std::cout << "run:" << elap << " us/evaluation\n";
    std::cout << "err/call:" << errs / 1e6 << " \n";
}

int main() {
    
    for (int d = 1; d < 5; ++d) {
        std::cout << "************ " << d << " ******************\n";
        time_one(d);
    }
}