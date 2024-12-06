
///
/// @file root-finding.cpp A comparison of univariate root-finding algorithms
/// @author M. Isgrò
///

#include <fstream>
#include <iomanip>

#define THEORETICA_LONG_DOUBLE_PRECISION
#include "theoretica.h"
#include "utility.h"
using namespace th;


template<typename Number>
Number ising(Number x) {

    return 2 * x * square(
        (th::exp(3 * x) + th::exp(-x)) /
        (th::exp(3 * x) + 3 * th::exp(-x))
    ) - x;
}


template<typename Number>
Number atanpoly(Number x) {

    return square(x - 1) * atan(x - 1);
}


template<typename Number>
Number noisyline(Number x) {

    return x + sin(100 * x) / 10;
}


// The function to find the root of. A function for the search
// of fixed points of the Ising model is used here.
template<typename Number>
Number f(Number x) {

    return atanpoly(x);
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
    const real u = s.Re() / s.Dual1();

    x = x - u - square(u) * s.Dual2() / (2.0 * s.Dual1());

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


// Ostrowski's method using automatic differentiation
real ostrowski(dual(*f)(dual), real guess) {

    real x = guess;

    const dual f_x = f(dual(x, 1));
    const real u = f_x.Re() / f_x.Dual();
    const real f_xu = f(dual(x - u, 0)).Re();

    x = x - u - (f_xu / f_x.Dual()) * (f_x.Re() / (f_x.Re() - 2 * f_xu));
    return x;
}


// Jarrat's method using automatic differentiation
real jarrat(dual(*f)(dual), real guess) {

    real x = guess;

    const dual f_x = f(dual(x, 1));
    const real u = f_x.Re() / f_x.Dual();
    const real f_xu = f(dual(x - 2.0 * u / 3.0, 1)).Dual();

    x = x - 0.5 * u + f_x.Re() / (f_x.Dual() - 3 * f_xu);
    return x;
}


// Bisection root finding
inline real bisection(real(*f)(real), real a, real b, real tol) {

    std::cout << "Running bisection method..." << std::endl;
    std::ofstream file ("bisection.dat");

    real x_avg = a;
    real x_min = a;
    real x_max = b;

    unsigned int iter = 0;

    while((x_max - x_min) > 2 * tol && iter <= OPTIMIZATION_BISECTION_ITER) {

        x_avg = (x_max + x_min) / 2.0;

        if(f(x_avg) * f(x_min) > 0)
            x_min = x_avg;
        else
            x_max = x_avg;

        std::cout << "[" << x_min << ", " << x_max << "]" << std::endl;
        file << iter << "\t" << std::setprecision(16)
             << log10(abs((x_min + x_max) / 2.0)) << std::endl;

        ++iter;
    }

    std::cout << "Bisection Iter. = " << iter << std::endl;

    if(iter > OPTIMIZATION_BISECTION_ITER) {
        return nan();
    }

    return x_avg;
}


// The Iterate-Truncate-Project method on a bracketing interval [a, b]
real itp(real(*f)(real), real a, real b, real eps, unsigned int n0 = 1, real k1 = 0) {

    std::cout << "Running ITP method..." << std::endl;

    std::ofstream file ("itp.dat");

    if (k1 == 0)
        k1 = 0.2 / (b - a);

    real y_a = f(a);
    real y_b = f(b);

    real x_t, x_new;

    const unsigned int n_half = th::floor(th::log2((b - a) / eps));
    const unsigned int n_max = n_half + n0;

    real eps_new = eps * (1 << n_max);
    unsigned int iter = 0;

    while((b - a) > (2 * eps) && iter <= OPTIMIZATION_BISECTION_ITER) {


        // Interpolation
        const real x_f = (a * y_b - b * y_a) / (y_b - y_a);
        const real x_half = (b + a) / 2.0;


        // Truncation
        const real sigma = th::sgn(x_half - x_f);
        const real delta = k1 * square(b - a); 

        if (delta <= abs(x_half - x_f))
            x_t = x_f + sigma * delta;
        else
            x_t = x_half;


        // Projection
        const real r = eps_new - (b - a) / 2.0;

        if (abs(x_t - x_half) <= r)
            x_new = x_t;
        else
            x_new = x_half - sigma * r;


        // Update
        const real y_new = f(x_new);

        if (y_new > 0) {
            b = x_new;
            y_b = y_new;
        } else if (y_new < 0) {
            a = x_new;
            y_a = y_new;
        } else {
            return x_new;
        }

        std::cout << "[" << a << ", " << b << "]" << std::endl;
        file << iter << "\t" << std::setprecision(16)
             << log10(abs((a + b) / 2.0)) << std::endl;

        eps_new /= 0.5;
        iter++;
    }

    std::cout << "ITP Iter. = " << iter << std::endl;

    return (b + a) / 2.0;
}


// Test the convergence of an iterative root finding method
template<typename IterFunction, typename TestFunction>
void test(IterFunction map, TestFunction f, real guess, std::string filename, size_t iter, real sol) {

    std::cout << "Writing " << filename << "..." << std::endl;

    std::ofstream file (filename);
    file << 0 << "\t" << std::fixed << std::setprecision(20);
    file << log10(abs(guess - sol)) << std::endl;

    real x = guess;
    for (unsigned int i = 0; i < iter; ++i) {
        
        x = map(f, x);
        file << (i + 1) << "\t" << std::fixed << std::setprecision(20);
        file << log10(abs(x - sol)) << std::endl;
    }
}


int main() {

    // Initial guess
    const real guess = 0.5;

    // Maximum number of iterations
    const size_t iter = 10;

    // Exact solution
    const real x0 = 1;

    std::cout.precision(16);

    test(newton, f<dual>, guess, "newton.dat", iter, x0);
    test(halley, f<dual2>, guess, "halley.dat", iter, x0);
    test(chebyshev, f<dual2>, guess, "chebyshev.dat", iter, x0);
    test(chun, f<dual>, guess, "chun.dat", iter, x0);
    test(ostrowski, f<dual>, guess, "ostrowski.dat", iter, x0);
    test(jarrat, f<dual>, guess, "jarrat.dat", iter, x0);

    std::cout << "Bisection Result = "
        << bisection(noisyline<real>, -1.0, +2.0, 1E-08) << std::endl;

    std::cout << "ITP Result = "
        << itp(noisyline<real>, -1.0, +2.0, 1E-08, 1) << std::endl;
}
