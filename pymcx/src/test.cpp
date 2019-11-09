#include "MultiComplex/MultiComplex.hpp"

void time_one(int d) {
    double DELTA = 1e-100, rx = 0.1234;

    std::valarray<double> r(0.0, exp2i(d));
    r[0] = rx;
    r[1] = DELTA;

    double exact = rx* cos(rx) + sin(rx);
    MultiComplex<double> x1(r);

    double errs = 0;
    auto startTime = std::chrono::system_clock::now();
    for (auto ii = 0; ii < 1e6; ++ii) {
        // A real number and its complex step
        errs += (x1*sin(x1))[exp2i(d)-1]/DELTA - exact;
    }

    auto endTime = std::chrono::system_clock::now();
    double elap = std::chrono::duration<double>(endTime - startTime).count();
    std::cout << "run:" << elap << " us/evaluation\n";
    std::cout << "err/call:" << errs / 1e6 << " \n";
}

int main() {
    MultiComplex<double> n{{-5291.570474614246, 3.8335011232881095e-67, 3.8335011232881095e-67, -5.5543929473310135e-137}};
    auto nexp = exp(n);
    for (int d = 1; d <= 6; ++d) {
        std::cout << "************ " << d << " ******************\n";
        time_one(d);
    }
}
