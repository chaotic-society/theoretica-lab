
///
/// @file root-finding.cpp A comparison of univariate root-finding algorithms
/// @author M. Isgrò
///

#include <fstream>
#include <iomanip>

#define THEORETICA_LONG_DOUBLE_PRECISION
#include "theoretica.h"
using namespace th;


// A function for the search of fixed points of the Ising model is used here.
template<typename Number>
Number f(Number x) {

    return 2 * x * square(
        (th::exp(3 * x) + th::exp(-x)) /
        (th::exp(3 * x) + 3 * th::exp(-x))
    ) - x;
}


// Newton's method using automatic differentiation
real newton(dual(*f)(dual), real guess) {

    real x = guess;
        
    const dual f_x = f(dual(x, 1));
    x = x - f_x.Re() / f_x.Dual();

    return x;
}

// Halley's method using automatic differentiation
real halley(dual2(*f)(dual2), real guess) {

    real x = guess;

    const dual2 f_x = f(dual2(x, 1));

    x = x - (2 * f_x.Re() * f_x.Dual1())
        / (2 * square(f_x.Dual1()) - f_x.Re() * f_x.Dual2());

    return x;
}


// Chebyshev's method using automatic differentiation
real chebyshev(dual2(*f)(dual2), real guess) {

    real x = guess;

    const dual2 s = f(dual2(x, 1));
    const real f_x = s.Re();
    const real df_x = s.Dual1();
    const real d2f_x = s.Dual2();

    x = x - (f_x / df_x) * (1 - f_x / (2 * d2f_x));

    return x;
}


// Chun's method using automatic differentiation
real chun(dual(*f)(dual), real guess) {

    real x = guess;

    const dual f_x = f(dual(x, 1));

    const real x_s = x - f_x.Re() / f_x.Dual();
    const dual f_xs = f(dual(x_s, 1));

    x = x - f_x.Re() / f_x.Dual()
          - 2 * f_xs.Re() / f_x.Dual()
          + f_xs.Re() * f_xs.Dual() / square(f_x.Dual());

    return x;
}


// Test the convergence of an iterative root finding method
template<typename IterFunction, typename TestFunction>
void test(IterFunction map, TestFunction f, real guess, std::string filename, size_t iter) {

    std::cout << "Writing " << filename << "..." << std::endl;

    std::ofstream file (filename);
    file << 0 << "\t" << std::fixed << std::setprecision(20) << guess << std::endl;

    real x = guess;
    for (unsigned int i = 0; i < iter; ++i) {
        
        x = map(f, x);
        file << (i + 1) << "\t" << std::fixed << std::setprecision(20) << x << std::endl;
    }
}


int main() {

    // Initial guess
    const real guess = 0.5;

    // Maximum number of iterations
    const size_t iter = 10;

    test(newton, f<dual>, guess, "newton.csv", iter);
    test(halley, f<dual2>, guess, "halley.csv", iter);
    test(chebyshev, f<dual2>, guess, "chebyshev.csv", iter);
    test(chun, f<dual>, guess, "chun.csv", iter);
}
