#ifndef __INVGAUSS_HPP
#define __INVGAUSS_HPP

#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <tuple>

template < typename Real_Type >
class inverse_gaussian_distribution
{
public:
    typedef Real_Type result_type;
    typedef std::pair< Real_Type, Real_Type > param_type;

    inverse_gaussian_distribution(const param_type& p) : _p(p) {}
    inverse_gaussian_distribution(Real_Type mu, Real_Type lambda) : inverse_gaussian_distribution(std::make_pair(mu, lambda)) {}
    inverse_gaussian_distribution() : inverse_gaussian_distribution(1.0, 1.0) {}

    template < class Generator >
    Real_Type operator () (Generator& g, const param_type& p) { return generate(g, p.first, p.second); }
    template < class Generator >
    Real_Type operator () (Generator& g) { return generate(g, _p.first, _p.second); }

    void reset() const {}

    Real_Type mean() const { return _p.first; }
    Real_Type shape() const { return _p.second; }

    param_type param() const { return _p; }
    void param(const param_type& p) { _p = p; }

    Real_Type min() const { return 0.0; }
    Real_Type max() const { return std::numeric_limits< Real_Type >::max(); }

    friend bool operator == (const inverse_gaussian_distribution& lhs, const inverse_gaussian_distribution& rhs) { return lhs.params() == rhs.params(); }
    friend bool operator != (const inverse_gaussian_distribution& lhs, const inverse_gaussian_distribution& rhs) { return not (lhs == rhs); }

    friend std::ostream& operator << (std::ostream& os, const inverse_gaussian_distribution& d)
    {
        os << d.p.first << " " << d.p.second;
        return os;
    }
    friend std::istream& operator >> (std::istream& is, inverse_gaussian_distribution& d)
    {
        is >> d._p.first >> d._p.second;
        return is;
    }

private:
    param_type _p;

    template < class Generator >
    static Real_Type generate(Generator& g, Real_Type mu, Real_Type lambda)
    {
        auto v = std::normal_distribution< Real_Type >()(g);
        auto y = v * v;
        auto x = (mu + ( mu * mu * y ) / ( lambda * 2.0 )
                  - ( mu * std::sqrt( mu * y * lambda * 4.0 + mu * mu * y * y ) ) / ( lambda * 2.0 ) );
        auto z = std::uniform_real_distribution< Real_Type >()(g);
        if (z <= mu / (mu + x))
        {
            return x;
        }
        else
        {
            return ( mu * mu ) / x;
        }
    }
}; // class inverse_gaussian_distribution

#endif

#ifdef INVGAUSS_SAMPLE
/*

Compile:

g++ -std=c++11 -D INVGAUSS_SAMPLE -x c++ invgauss.hpp -o invgauss

Run:

./invgauss 1 40 1000 42 >ig.1.40.txt

Visualize with matplotlib:

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def inv_gauss(x, m, l):
    return sp.sqrt(l / (2 * sp.pi * x * x * x)) * sp.exp(- l * (x - m) * (x - m) / (2 * m * m * x))

x = np.loadtxt(open('ig.1.40.txt'))
b = np.linspace(.4, 3, 100)
cb = b[:-1] + (b[1] - b[0])/2
y = inv_gauss(cb, 1, 40)

plt.figure(0)
plt.hist(x, bins=b, normed=1, histtype='step')
plt.gca().set_color_cycle(None)
plt.plot(cb, y, '--')

*/

#include <cassert>
#include <chrono>
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        cerr << "use: " << argv[0] << " mu lambda n [seed]" << endl;
        exit(EXIT_FAILURE);
    }
    float mu = 0.0;
    float lambda = 0.0;
    size_t n = 0;
    size_t seed = 0;
    istringstream(argv[1]) >> mu;
    clog << "mu: " << mu << endl;
    assert(mu > 0.0);
    istringstream(argv[2]) >> lambda;
    clog << "lambda: " << lambda << endl;
    assert(lambda > 0.0);
    istringstream(argv[3]) >> n;
    clog << "n: " << n << endl;
    if (argc >= 5)
    {
        istringstream(argv[4]) >> seed;
    }
    if (seed == 0)
    {
        seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    clog << "seed: " << seed << endl;
    std::mt19937 rg(seed);
    inverse_gaussian_distribution< float > ig(mu, lambda);
    for (size_t i = 0; i < n; ++i)
    {
        cout << ig(rg) << endl;
    }
}

#endif
