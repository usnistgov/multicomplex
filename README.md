# pymcx

pymcx is a Python-based multi complex algebra library.  This library is used to carry out numerical derivatives of (semi-)arbitrary numerical functions with nearly machine precision.

At its core is a minimal header-only C++11 library that does the heavy lifting, and pybind11 is used to make a 1-to-1 interface between the C++ code and the Python interface

## Examples:

Try it in your browser: [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ianhbell/mcx/master)

Here is how to calculate the first 10 derivatives of sin(x) to numerical precision:

``` Python
import pymcx
import numpy
derivs = pymcx.diff_mcx1(lambda x: np.sin(x), 0.1234, 10)
```

WARNING: Not all functions are implemented, due to the complications of multicomplex algebra.  The following functinons are implemented:

* Element-wise operations: *,+,/,-
* Trig functions: cos, sin, cosh, sinh
* Transcendental functions: ln, exp

## More Reading Material

Info from around the interwebs:

* http://doi.acm.org/10.1145/2168773.2168774
* http://adl.stanford.edu/hyperdual/Fike_AD2012_slides.pdf
* http://folk.ntnu.no/preisig/HAP_Specials/AdvancedSimulation_files/2014/AdvSim-2014__Verheule_Adrian_Complex_differenetiation.pdf

## License

*MIT licensed (see LICENSE for specifics), not subject to copyright in the USA. Foreign Rights Reserved, Secretary of Commerce.

## Contacts

The C++ and Python code was written by Ian Bell (of NIST), with special thanks to Ulrich Deiters (University of Cologne) and Bradley Alpert (of NIST)

## Dependencies

* Unmodified [pybind11](https://github.com/pybind/pybind11) for C++ <-> Python interfacing

## Contributing/Getting Help

If you would like to contribute to ``pymcx`` or report a problem, please open a pull request or submit an issue.  Especially welcome would be additional tests.

If you want to discuss or request assistance, please open an issue.

To get started, you should check out the Jupyter notebook ``Tutorial`` in the notebooks folder; they demonstrate some of the capabilities of the Python interface.

## Installation

### Prerequisites

You will need:

* git
* cmake (on windows, install from cmake, on linux ``sudo apt install cmake`` should do it, on OSX, ``brew install cmake``)
* Python (the anaconda distribution is used by the authors)
* a compiler (on windows, Visual Studio 2017+, g++ on linux/OSX)

If on linux you use Anaconda and end up with an error something like
```
ImportError: /home/theuser/anaconda3/bin/../lib/libstdc++.so.6: version `GLIBCXX_3.4.20' not found (required by /home/theuser/anaconda3/lib/python3.6/site-packages/pymcx.cpython-35m-x86_64-linux-gnu.so)
```
it can be sometimes fixed by installing ``libgcc`` with conda: ``conda install libgcc``.  [This is due to an issue in Anaconda](https://github.com/ContinuumIO/anaconda-issues/issues/483)

## To install in one line from github (easiest)

This will download the sources into a temporary directory and build and install the python extension so long as you have the necessary prerequisites:
```
pip install git+git://github.com/ianhbell/mcx.git
```

### From a cloned repository

Alternatively, you can clone (recursively!) and run the ``setup.py`` script

```
git clone --recursive https://github.com/ianhbell/mcx
cd mcx
pip -vvv install .
```

to install, or 

```
pip -v -e install .
```

to use a locally-compiled version for testing.  If you want to build a debug version, you can do so with

```
python setup.py build -g develop
```
With a debug build, you can step into the debugger to debug the C++ code, for instance.  

### Cmake build for C++

Starting in the root of the repo (a debug build with the default compiler, here on linux):

``` 
git clone --recursive https://github.com/ianhbell/mcx
cd mcx/pymcx
mkdir build
cd build
cmake ..
cmake --build .
```
For those using Anaconda on Linux, please use the following for cmake:
```
mkdir build
cd build
cmake .. -DPYTHON_EXECUTABLE=`which python`
cmake --build .
```
For Visual Studio 2015 (64-bit) in release mode, you would do:
``` 
git clone --recursive https://github.com/ianhbell/mcx
cd mcx/pymcx
mkdir build
cd build
cmake .. -G "Visual Studio 14 2015 Win64"
cmake --build . --config Release
```

If you need to update your submodules (pybind11 and friends)

```
git submodule update --init
```

For other options, see the cmake docs