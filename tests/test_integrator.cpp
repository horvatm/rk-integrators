/*
  Testing the integrator via mathematical pendulum.

  Author: Martin Horvat, November 2013
*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

//#define EMBEDDED
#include "integrator.h"

#include "print_array.h"
#include "phys_model.h"

int main(int argc, char **argv ){

  if (argc != 4){
    std::cerr
      << "usage: ./test_integrator <q:p> <epsA:epsR> <dt:T> > <energy at steps>\n"
      << "note:\n"
      << "  q:p' -- initial point\n"
      << "  dt:T -- time step and whole interval\n"
      << "  epsA:epsR -- absolute and relative error\n";

    exit(EXIT_FAILURE);
  }


  int dim = 2;

  double eps_[2], dt_, T_, x_[dim];

  sscanf(argv[1], "%lf:%lf", x_, x_ + 1);

  sscanf(argv[2], "%lf:%lf", eps_, eps_ + 1);

  sscanf(argv[3], "%lf:%lf", &dt_, &T_);

  real t = 0, eps[2], dt, T, x[dim];

  eps[0] = eps_[0];
  eps[1] = eps_[1];

  x[0] = x_[0];
  x[1] = x_[1];

  dt = dt_;
  T = T_;

  perturbations::Perturbation p;

  typedef METHOD<real, real, perturbations::Perturbation> TMethod;

  TPertInt<real, real, TMethod, perturbations::Perturbation> integrate (&p);

  real H = energy(x);

  std::cout.precision(16);
  std::cout << std::scientific;

  #if 1     // point by point simulation

  int steps;

  real h = -1;

  std::cout
    << t << ' ' << W<real>(dim, x, " ") << ' '
    << 0 << ' ' << 0 << '\n';

  while (t <= T) {
    steps = integrate.step(eps, dt, h, t, x);
    std::cout
      << t << ' ' << W<real>(dim, x, " ") << ' '
      << energy(x) - H  << ' ' << steps << '\n';
  }

  #else   // calculating the whole trajectory

  int n = int(T/dt),
      *steps = new int [n+1];

  real *y = new real [(n+1)*dim];

  y[0] = x[0];
  y[1] = x[1];
  steps[0] = 0;

  // to optimize memory consumption
  // we could run:
  //   integrate.traj(eps, dt, t, n, y);

  integrate.traj(eps, dt, t, n, y, steps+1);

  real *q = y;

  for (int i = 0; i <= n; ++i) {
    std::cout << t << ' ' << W<real>(dim, q, " ") << ' '
              << energy(q) - H  << ' ' << steps[i] << '\n';
    q += dim;
    t += dt;
  }


  delete [] y;
  delete [] steps;

  #endif

  return EXIT_SUCCESS;
}
