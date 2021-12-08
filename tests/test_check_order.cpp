/*
  Checking the order of integrator methods via mathematical pendulum.

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

  if (argc != 3){
    std::cerr
      << "usage: ./test_check_order <q:p> <dt0:dt1:N> > <dt denergy>\n"
      << "note:\n"
      << "  q:p' -- initial point\n"
      << "  dt0:dt1:N -- time step interval\n";

    exit(EXIT_FAILURE);
  }

  int dim = 2, N;

  double dt[2], x[dim];

  sscanf(argv[1], "%lf:%lf", x, x + 1);

  sscanf(argv[2], "%lf:%lf:%d", dt, dt + 1, &N);

  perturbations::Perturbation p;

  METHOD<real, real, perturbations::Perturbation> method(&p);

  std::cout.precision(16);
  std::cout << std::scientific;

  int i, j;

  real t, y[2], z[2],
         h = dt[0], h_, H,
         fac = std::exp(std::log(dt[1]/dt[0])/N);

  for (i = 0; i < N; ++i){

    t = 0;

    for (j = 0; j < 2; ++j) z[j] = y[j] = x[j];

    H = energy(y);

    method.step(h, t, y);

    h_ = h/10;

    for (j = 0; j < 10; ++j) {
      method.step(h_, t, z);
      t += h_;
    }

    std::cout
      << h << ' '
      << energy(y) - H << ' '
      << std::max(std::abs(y[0] - z[0]), std::abs(y[1] -z[1])) << '\n';

    h *= fac;
  }

  return EXIT_SUCCESS;
}
