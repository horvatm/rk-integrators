/*
  Example using integrator with the complex-number derivative defined via
  lambda function.

  Author: Martin Horvat, November 2021
*/

#include <iostream>
#include <cmath>
#include <complex>

#include "integrator_lambda.h"

using dcomplex = std::complex<double>;

int main(){

  int dim = 2;

  auto deriv = [& ] (double t, dcomplex y[2], dcomplex dydt[2]) {
    dydt[0] = y[1];
    dydt[1] = -y[0];
  };

  //using TMethod  = NonThreadSave::NonEmbedded::Hairer10<double, dcomplex, TDerivative>;
  //using TMethod  = NonThreadSave::NonEmbedded::Hairer8<double, dcomplex, TDerivative>;
  using TMethod = NonThreadSave::NonEmbedded::RKButcher<double, dcomplex>;
  //using TMethod  = NonThreadSave::NonEmbedded::RK4<double, dcomplex, TDerivative>;
  //using TMethod = NonThreadSave::NonEmbedded::RK2<double, dcomplex>;

  //
  // Testing stepping
  //

  if (1) {
    TMethod method(dim);

    double
        t = 0,
        h = 0.01;

    dcomplex
        y[2] = {1,0};

      for (int i = 0; i < 1000; ++i) {

        std::cout << t << '\t' << y[0] << '\t' << y[1] << '\n';
        method.step(h, t, y, deriv);
        t += h;
      }
  }
  //
  // Testing integrator
  //

  if (0) {

    using Tintegrator = NonThreadSave::NonEmbedded::TIntegrator <double, dcomplex, TMethod>;

    Tintegrator integrator(dim);

    int n = 1000;

    double
      h = 0.01,
      eps[2] = {1e-12, 1e-12};

    auto
      *y = new dcomplex [2*(n+1)];

    y[0] = 1;
    y[1] = 0;

    integrator.traj(eps, h, 0, n, y, deriv);

    for (int i = 0; i < n; ++i)
      std::cout << i*h << '\t' << y[2*i] << '\t' << y[2*i + 1] << '\n';

    delete [] y;
  }
  return 0;
}
