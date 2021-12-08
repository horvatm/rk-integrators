/*
  Testing if a integration step works. Storing points and derivatives.

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
      << "usage: ./test_method <q:p> <dt> <N> > <t point energy>\n"
      << "note:\n"
      << "  q:p' -- initial point\n"
      << "  dt -- time step interval\n"
      << "  N -- number of steps\n";

    exit(EXIT_FAILURE);
  }

  int dim = 2, N;

  double dt, x[dim], H;

  sscanf(argv[1], "%lf:%lf", x, x + 1);

  dt = atof(argv[2]);

  N = atoi(argv[3]);

  perturbations::Perturbation p;

  METHOD<real, real, perturbations::Perturbation> method(&p);

  std::cout.precision(16);
  std::cout << std::scientific;

  int i;

  real t = 0, H0 = energy(x), y[2*dim];

  for (i = 0; i < dim; ++i) y[i] = x[i];

  for (i = 0; i < N; ++i){

    H = energy(y);

    method.calc_derivative(t, y, y+dim);

    method.step_(dt, t, y, y+dim);

    std::cout << t << ' ' << W<real>(2*dim, y, " ") << ' ' << H-H0 << '\n';

    t += dt;
  }

  return EXIT_SUCCESS;
}
