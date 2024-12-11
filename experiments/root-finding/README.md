# Numerical Experiment: Univariate Root-Finding
Many different numerical methods for univariate root-finding have been implemented in C++ and compared by studying convergence on different functions. Automatic differentiation was used to compute exact derivatives of the functions, as implemented in Theoretica. The general problem of root-finding corresponds to solving the nonlinear equation:

$$f(x^*) = 0$$

which in the context of numerical analysis consists in finding a $x_n$ so that $|x^* - x_n| < \epsilon$ or alternatively $|f(x_n)| < \epsilon$, where { $x_n$ } is generally a succession of values given by an iterative method, defined by a recursion relation.

## Methods
The following methods have been chosen for comparison.
#### Base methods:
- **Newton's method**: the most common method for root-finding when the first derivative is available.
- **Halley's method**: an extension of Newton's method which uses the second derivative.
- **Chebyshev's method**: a method similar to Halley's, which uses the second derivative.
#### Higher Order:
- **Ostrowski's method**: an advanced 4-th order method using 2 function evaluations.
- **Jarrat's method**: an advanced 4-th order method using 2 derivative evaluations.
- **Chun's method**: a higher order method making use of 2 function and derivative evaluations.
#### Derivative-Free:
- **Bisection method**: the most common derivative-free bracketing method.
- **ITP method**: a novel bracketing method with improved performance with respect to bisection.

The methods have been selected as good candidates for univariate root-finding of generic real functions, with derivative-free methods being a fallback for particularly ill-conditioned functions. The ITP method's parameters have been tuned to provide a good balance of performance and accuracy, as the generic method uses some computationally expensive functions such as _log2_ and _powf_. The fixed parameters were $k_2 = 2$ and $n_0 = 1$, which avoid the _powf_ and have good convergence for smooth functions.

## Experimental Approach
The following functions $f_1$ ("ising") and $f_2$ ("atanpoly") have been chosen for comparison of derivative methods, while "f_3" ("noisyline") has been used to test derivative-free methods. The initial guess has been chosen so that $x_g \in [x^* - 0.5, x^* + 0.5]$, to highlight the asymptotic behavior. The value of $|x^* - x_n|$ was then computed for an increasing number of iterations and plotted. In the case of bracketing, derivative-free methods, the value $|x^* - (x_l + x_r)/2|$ has been used instead.

## Results

![Root Finding - ising](https://github.com/chaotic-society/theoretica-lab/blob/main/experiments/root-finding/root-finding-ising.png?raw=true)

![Root Finding - atanpoly](https://github.com/chaotic-society/theoretica-lab/blob/main/experiments/root-finding/root-finding-atanpoly.png?raw=true)

As expected, Newton's method had the slowest convergence over both functions, but also the smallest number of function evaluations for each step of all methods.
Halley's and Chebyshev's method did similarly, with a faster convergence, meanwhile Chun's method did only slightly better then these two methods, while requiring an additional function evaluation. The higher order methods, Ostrowski's and Jarrat's, had also similar convergence and the fastest of all methods.

![Root Finding - NoisyLine](https://github.com/chaotic-society/theoretica-lab/blob/main/experiments/root-finding/root-finding-noisyline-bisect.png?raw=true)

![Root Finding - Ising](https://github.com/chaotic-society/theoretica-lab/blob/main/experiments/root-finding/root-finding-ising-bisect.png?raw=true)

The ITP method was found to follow the theoretical expectations, with an extremely faster convergence over smooth functions with respect to bisection. The choice of the $\epsilon$ parameters has been found to greatly influence the convergence of the algorithm, delaying its convergence.

## Conclusions
All of the methods have been selected for implementation in Theoretica with the exception of Chun's method, which in this implementation has been found to be sub-optimal with respect to its number of function evaluations. The ITP was found to be generally better than the bisection, as expected from theoretical results.
