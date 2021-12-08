/*
  Example using integrator with the real-number derivative defined via
  lambda function.

  Author: Martin Horvat, November 2021
*/


#include <iostream>
#include <cmath>

#include "integrator_lambda.h"

int main(){

  int dim = 2;

  auto deriv = [& ] (double t, double y[2], double dydt[2]) {
    dydt[0] = y[1];
    dydt[1] = -y[0];
  };

  //using TMethod  = NonThreadSave::NonEmbedded::Hairer10<double, double, TDerivative>;
  //using TMethod  = NonThreadSave::NonEmbedded::Hairer8<double, double, TDerivative>;
  using TMethod = NonThreadSave::NonEmbedded::RKButcher<double, double>;
  //using TMethod  = NonThreadSave::NonEmbedded::RK4<double, double, TDerivative>;
  //using TMethod = NonThreadSave::NonEmbedded::RK2<double, double>;


  //
  // Testing stepping
  //

  if (1) {


    using TMethod = NonThreadSave::NonEmbedded::RK2<double, double>;

    TMethod method(dim);

    double
        t = 0,
        h = 0.01,
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
    using Tintegrator = NonThreadSave::NonEmbedded::TIntegrator <double, double, TMethod>;

    Tintegrator integrator(dim);

    int n = 1000;

    double
      h = 0.01,
      eps[2] = {1e-12, 1e-12},
      *y = new double [2*(n+1)];

    y[0] = 1;
    y[1] = 0;

    integrator.traj(eps, h, 0, n, y, deriv);

    for (int i = 0; i < n; ++i)
      std::cout << i*h << '\t' << y[2*i] << '\t' << y[2*i + 1] << '\n';

    delete [] y;
  }
  return 0;
}
