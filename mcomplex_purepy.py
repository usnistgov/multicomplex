import numpy as np
import numbers

class MultiComplex:
    """

    Special methods: https://docs.python.org/3/reference/datamodel.html#emulating-numeric-types
    """
    def __init__(self, *, a=None, b=None, d=1, v=None):

        if v is not None:
            L = len(v)
            if L < 1:
                raise ValueError("Doesn't make sense")
            elif L == 1:
                self.a = v[0]
                self.b = 0
                self.d = 1
            elif L == 2:
                self.a = v[0]
                self.b = v[1]
                self.d = 1
            elif L % 2 != 0:
                raise ValueError("Length of coeffs [{0:d}] must be power of 2".format(L))
            else:
                r = MultiComplex(v=v[0:L//2])
                s = MultiComplex(v=v[L//2::])
                self.a = r
                self.b = s
                self.d = r.d+1
        else:
            self.a = a
            self.b = b
            self.d = d

    def __array_ufunc__(self, ufunc, method, *args, **kwargs):
        """
        Overload some methods from numpy to handle multicomplex numbers

        N.B.: Not all numpy methods are implemented as ufunc (e.g., np.real), and therefore, 
        they will just return the object back again.

        The expansion for sin goes something like this in sympy:
        >>> a, b = symbols('a,b', real=True)
        >>> sin(a + I*b).expand(complex=True)
        """
        if method == '__call__':
            z = self
            if ufunc is np.cos:
                a, b = reim(z); return MultiComplex(a=np.cos(a)*np.cosh(b),b=-np.sin(a)*np.sinh(b), d=dim(z))
            elif ufunc is np.sin:
                a, b = reim(z); return MultiComplex(a=np.sin(a)*np.cosh(b), b=np.cos(a)*np.sinh(b), d=dim(z))
            elif ufunc is np.sinh:
                a, b = reim(z); return MultiComplex(a=np.sinh(a)*np.cos(b), b=np.cosh(a)*np.sin(b), d=dim(z))
            elif ufunc is np.cosh:
                a, b = reim(z); return MultiComplex(a=np.cosh(a)*np.cos(b), b=np.sinh(a)*np.sin(b), d=dim(z))
            elif ufunc is np.exp:
                a, b = reim(z); return np.exp(a) * MultiComplex(a=np.cos(b), b=np.sin(b), d=dim(z))

    def __neg__(self):
        z = self
        return MultiComplex(a=-real(z), b=-imag(z), d=dim(z))

    def __mul__(self, w):
        z = self
        if dim(z) < dim(w):
            return MultiComplex(a=z*real(w), b=z*imag(w), d=dim(w)) 
        elif dim(z) > dim(w):
            return MultiComplex(a=real(z)*w, b=imag(z)*w, d=dim(z))
        else:
            return MultiComplex(a=real(z)*real(w)-imag(z)*imag(w),
                                b=real(z)*imag(w)+imag(z)*real(w), d=dim(z))

    def __sub__(self, w):
        z = self
        if dim(z) < dim(w):
            return MultiComplex(a=z-real(w),b=-imag(w),d=dim(w)) 
        elif dim(z) > dim(w):
            return MultiComplex(a=real(z)-w,b=imag(z),d=dim(z))
        else:
            return MultiComplex(a=real(z)-real(w),b=imag(z)-imag(w),d=dim(z))

    def __add__(self, w):
        z = self
        if isinstance(w, numbers.Number):
            return MultiComplex(a=w+real(z),b=w+imag(z),d=dim(z))
        elif dim(z) < dim(w):
            return MultiComplex(a=z+real(w),b=imag(w),d=dim(w))
        elif dim(z) > dim(w):
            return MultiComplex(a=real(z)+w,b=imag(z),d=dim(z))
        else:
            return MultiComplex(a=real(z)+real(w),b=imag(z)+imag(w),d=dim(z))

    def __getitem__(self, key=None):
        z = self
        if isinstance(key, tuple):
            if dim(z) == 1 and len(key) == 1:
                if key[0] == 0:
                    return real(z)
                elif key[0] == 1:
                    return imag(z)
            elif len(key) > 1 and key[0] == 0:
                return real(z)[key[1::]]
            elif len(key) > 1 and key[0] == 1:
                return imag(z)[key[1::]]
            else:
                raise IndexError("MultiComplex key length of {0:d} must match dim of {1:d}".format(len(key), dim(z)))
        else:
            raise TypeError("Input must be tuple")

    def __repr__(self):
        return str((self.a, self.b, self.d))

def real(z):
    return z.a
def imag(z):
    return z.b
def reim(z):
    return (z.a, z.b)
def dim(z):
    return z.d

rx = 0.1234
m1 = MultiComplex(v=[rx, 1e-100,0,0])
def f(x):
    return np.sin(x)
def dfdx(x):
    return np.cos(x)
m1s = f(m1)
ind_mocom = tuple([0,1]) # most complex index
print(dfdx(rx), m1s[ind_mocom]/1e-100)
print(m1s[(0,0)])
print(m1s[(0,1)])
print(m1s[(1,0)])
print(m1s[(1,1)])

rx = 0.1
m1 = MultiComplex(v=[rx, 1e-100]+[0]*14)

def f(x):
    return x*np.sin(x)

def dfdx(x):
    return np.sin(x) + x*np.cos(x)

m1s = f(m1)
ind_mocom = tuple([0]*(m1.d-1)+[1]) # most complex index
print(dfdx(rx), m1s[ind_mocom]/1e-100)