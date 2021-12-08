#pragma once

/*
  Library for adaptive step evolution of Schwarzschield constants of motions

    (Q, P)  = QP[6]

  along lambda. Library provides

    class TPertInt

  using

  if defined(EMBEDDED)
    RKM4  -- Marson method
    Hairer8 -- Hairer (p=8) method
  #else
    RK4 -- standard RK4
    Hairer8 -- Hairer (p=8) method
    Hairer10 -- Hairer (p=10) method
  #endif

  For heuristic step-size control we use

  http://www.gnu.org/software/gsl/manual/html_node/Adaptive-Step_002dsize-Control.html

  where we set a_y = 1 and a_dydt=0

  Author: Martin Horvat, November 2013
*/

#include <cmath>
#include <limits>

// C-style matrix allocation and deallocation
#include "matrix.h"

// RK coefficients
#include "RKcoeff8d_1.h"
#include "RKcoeff10b_1.h"

//#define DEBUG

#if defined(EMBEDDED)

/* ============================================================================
  Calculating the evolution of Schwz. constants of motions gven by 1. order ODE
  using an adaptive step in a combination with some explicit EMBEDDED method of
  some order.

  The embedded give also an approximation of the error at each step.
============================================================================= */

/*
  General solver class:
    T  - numerical type for time
    F  - numerical type for states
    Derivative  - derivative in ODE (right side of the equation)
    Method  - explicit method of integration with given order
*/


template <class T, class F, class Method, class Perturbation>
class TPertInt {

  int dim;                      // dimension of phase space

  Method *method;               // pointer to a method

  F *s, *ds;                     // states needed in the steping

  public:

  // setup of adaptive algorithm
  T TINY, SUPERTINY;

  int PERIOD_OF_CHECKING, PERIOD_OF_NOTCHECKING;

  TPertInt(Perturbation *p) {

    dim = p->get_dim();

    method = new Method(p);
    s = new F [dim];
    ds = new F [dim];

    TINY = std::numeric_limits<F>::epsilon();
    SUPERTINY = TINY/1000;

    PERIOD_OF_CHECKING = 5;
    PERIOD_OF_NOTCHECKING = 20;
  }


  ~TPertInt (){
    delete method;
    delete [] s;
    delete [] ds;
  }

  /*
    Integrator of the ODE dx/dt = F(x). Making a controlled step over the
    length dt with suggested integration step h

    input:
      eps[2] -- abs. (accuracy) and relative (precision) permitted local error
      dt -- time step
      h -- suggested integration step
      t -- current time at which input state is given
      state -- position of the ODE at time t

    output:
      h -- last integration step
      t -- new time, t+dt
      state -- position of the ODE at t+dt

    return:
      number of steps done to integrate between [t,t+dt]
  */

  int step(T eps[2], const T &dt, T & h, T & t, F *state);

  /*
    Calculating the trajectory sampled with time step dt.

    input:
      eps[2] -- abs. (accuracy) and relative (precision) permitted local error
      dt -- time step
      t -- current time at which input state is given
      trajectory[dim] -- initial position x(t)

    output:
      trajectory[dim*(n+1)] - whole trajectory
                            = {x(t), x(t+dt), .., x(t+n*dt)}

      steps - number of steps done to integrate between [t+dt*i,t+dt*(i+1)]
              for i = 0, .., n-1
  */

  void traj(T eps[2], const T &dt, T t, const int &n, F *trajectory, int *steps=0);
};


template <class T, class F, class Method, class Perturbation>
int TPertInt <T, F, Method, Perturbation>::
    step(T eps[2], const T &dt, T & h, T & t, F *state) {

  int i, checked = 0, not_checked = 0, steps = 0;

  bool check = true;

  T A, E, D, tmp, h_out = h, t_end = t + dt;

  if (h <= 0)
    h_out = h = dt/4;
  else if (t + h > t_end)
    h = t_end - t;

  cp(s, state, dim);

  do {

    method->step(h, t, s, ds);

    if (check) {

      // max-norm of ds[i] -- estimate of abs. local error
      E = 0;
      for (i = 0; i < dim; ++i)
        if ((tmp = std::abs(ds[i])) > E) E = tmp;

      // max-norm of s[i]
      A = 0;
      for (i = 0; i < dim; ++i)
        if ((tmp = std::abs(s[i])) > A) A = tmp;

      // heuristic expression of permitted local error
      D = eps[0] + eps[1]*A;

      #if defined(DEBUG)
      std::cerr << t << '\t' << h << '\t' << E << '\t' << checked << '\n';
      #endif

      if (E > D) {
        h_out = (h *= 0.9*std::exp(log(D/E)/(method->order + 1)));
        cp(s, state, dim);
        checked = 0;
      } else {
        t += h;
        cp(state, s, dim);

        if (std::abs(t - t_end) < TINY*std::max(T(1), std::abs(t_end))) break;

        // is the step is too small but larger then zero
        if (SUPERTINY < E && E < D/2) {
          h_out = (h *= 0.9*std::exp(std::log(D/E)/(method->order+1)));
          checked = 0;
        } else if (++checked > PERIOD_OF_CHECKING) {
          check = false;
          not_checked = 0;
        }

        if (t + h > t_end) h = t_end - t;
      }

    } else {
      t += h;
      cp(state, s, dim);

      #if defined(DEBUG)
      std::cerr << t << '\t' << '\n';
      #endif

      if (std::abs(t - t_end) < TINY*std::max(T(1), std::abs(t_end))) break;

      if (++not_checked > PERIOD_OF_NOTCHECKING) {
        check = true;
        checked = 0;
      }

      if (t + h > t_end) h = t_end - t;
    }

    ++steps;
  } while (true);

  h = h_out;
  return steps;
}


template <class T, class F, class Method, class Perturbation>
void TPertInt<T, F, Method, Perturbation>::
  traj(T eps[2], const T &dt, T t, const int &n, F *trajectory, int *steps) {

  T h = -1;

  for (int i = 0; i < n; ++i){
    cp(trajectory + dim, trajectory, dim);
    trajectory += dim;
    if (steps)
      steps[i] = step(eps, dt, h, t, trajectory);
    else
      step(eps, dt, h, t, trajectory);
  }
}

/*
   Runge-Kutta-Merson method of 4th order.

   Ref: Bohte Z  Numericne metode  (DMFA, 1991)
*/

#define DERIVATIVE p->getQPdot

template <class T, class F, class Perturbation>
class RKM4 {

  int dim;

  Perturbation *p;

  public:

  F **k, *state_;

  static int order;

  RKM4 (Perturbation *p) : p(p) {

    dim = p->get_dim();

    k = matrix<F> (6, dim);
    state_ = k[5];
    /*
    k = new F* [5];
    for (int i = 0; i < 5; ++i) k[i] = new F [dim];
    state_ = new F [dim];
    */
  }

  ~RKM4(){
    free_matrix(k);
    /*
    for (int i = 0; i < 5; ++i) delete [] k[i];
    delete [] k;
    delete [] state_;
    */
  }

  void get_derivative(F *dy){ cp(dy, k[0], dim); }
  void calc_derivative(const F &t, F *y, F *dy){ DERIVATIVE(t, y, dy); }

  void step(const T &h, const T &t, F *state, F *dstate = 0){

    int i;

    DERIVATIVE(t, state, k[0]);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + (k[0][i] *= h)/3;

    DERIVATIVE(t + h/3, state_, k[1]);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + ((k[1][i] *= h) + k[0][i])/6;

    DERIVATIVE(t + h/3, state_, k[2]);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + (3*(k[2][i] *=h) + k[0][i])/8;

    DERIVATIVE(t + h/2, state_, k[3]);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + 2*(k[3][i] *= h) -3*k[2][i]/2 + k[0][i]/2;

    DERIVATIVE(t + h, state_, k[4]);

    if (dstate)
      for (i = 0; i < dim; ++i) {
        state[i] += (k[0][i] + 4*k[3][i] + (k[4][i] *= h))/6;
        dstate[i] = (2*k[0][i] - 9*k[2][i] + 8*k[3][i] - k[4][i])/30;
      }
    else
      for (i = 0; i < dim; ++i)
        state[i] += (k[0][i] + 4*k[3][i] + (k[4][i] *= h))/6;
  }

  void step_(const T &h, const T &t, F *state, F *k0 = 0){

    int i;

    if (k0 == 0)
      DERIVATIVE(t, state, k[0]);
    else if (k0 != k[0])
      cp(k[0], k0, dim);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + (k[0][i] *= h)/3;

    DERIVATIVE(t + h/3, state_, k[1]);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + ((k[1][i] *= h) + k[0][i])/6;

    DERIVATIVE(t + h/3, state_, k[2]);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + (3*(k[2][i] *=h) + k[0][i])/8;

    DERIVATIVE(t + h/2, state_, k[3]);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + 2*(k[3][i] *= h) -3*k[2][i]/2 + k[0][i]/2;

    DERIVATIVE(t + h, state_, k[4]);

    for (i = 0; i < dim; ++i)
      state[i] += (k[0][i] + 4*k[3][i] + (k[4][i] *= h))/6;
  }

};

template <class T, class F, class Perturbation>
int RKM4<T,F,Perturbation>::order = 4;

/*
   Hairer method (p=8)

   Ref: P.J. Prince and J.R. Dormand
*/

template <class T, class F, class Perturbation>
class Hairer8 {

  int dim;

  Perturbation *p;

  T **a, *b, *b_, *c;

  public:

  F **k, *state_;

  static int order;

  Hairer8 (Perturbation *p) : p(p) {
    int i;

    dim = p -> get_dim();

    k = matrix <F> (14, dim);

    state_ = k[13];

    /*
    k = new F* [13];
    for (i = 0; i < 13; ++i) k[i] = new F [dim];

    state_ = new F [dim];
    */

    a = new T* [12];

    for (i = 0; i < 12; ++i) a[i] = new T [i+1];

    b = new T [13];

    b_ = new T [13];

    c = new T [12];

    Hairer8_coef(c, a, b, b_);
  }

  ~Hairer8(){
    int i;

    delete [] b;
    delete [] b_;
    delete [] c;

    for (i = 0; i < 12; i++) delete [] a[i];
    delete [] a;

    free_matrix(k);
    /*
    delete [] state_;
    for (i = 0; i < 13; i++) delete [] k[i];
    delete [] k;
    */
  }

  void get_derivative(F *dy){ cp(dy, k[0], dim); }
  void calc_derivative(const F &t, F *y, F *dy){ DERIVATIVE(t, y, dy); }

  void step(const T &h, const T &t, F *state, F *dstate = 0){

    int i, j, l;

    F tmp, tmp_;

    DERIVATIVE(t, state, k[0]);

    for (i = 0; i < 12; ++i) {

      for (l = 0; l < dim; ++l) {

        tmp = 0;
        for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

        state_[l] = state[l] + h*tmp;
      }

      DERIVATIVE(t + c[i]*h, state_, k[i+1]);
    }

    if (dstate) {

      zero(dstate, dim);

      for (l = 0; l < dim; ++l) {
        tmp = tmp_ = 0.0;

        for (i = 0; i < 13; i++) {
          tmp  += b[i]*k[i][l];
          tmp_ += (b_[i]-b[i])*k[i][l];
        }

        state[l] += h*tmp;
        dstate[l] += h*tmp_;
      }
    } else
      for (l = 0; l < dim; ++l) {
        tmp = 0.0;

        for (i = 0; i < 13; i++) tmp  += b[i]*k[i][l];

        state[l] += h*tmp;
      }
  }

  void step_(const T &h, const T &t, F *state, F *k0 = 0){

    int i, j, l;

    F tmp, tmp_;

    if (k0 == 0)
      DERIVATIVE(t, state, k[0]);
    else if (k0 != k[0])
      cp(k[0], k0, dim);

    for (i = 0; i < 12; ++i) {

      for (l = 0; l < dim; ++l) {

        tmp = 0;
        for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

        state_[l] = state[l] + h*tmp;
      }

      DERIVATIVE(t + c[i]*h, state_, k[i+1]);
    }

    for (l = 0; l < dim; ++l) {
      tmp = 0.0;

      for (i = 0; i < 13; i++) tmp  += b[i]*k[i][l];

      state[l] += h*tmp;
    }
  }
};

template <class T, class F, class Perturbation>
int Hairer8<T,F,Perturbation>::order = 8;

#undef DERIVATIVE

#else // #if !defined(EMBEDDED)

/* ============================================================================
  Integrator of 1. order ODE defined as

    dx/dt = F(t, x)     x in T^d

  using an adaptive step in a combination with some explicit method of
  some order.
 ============================================================================ */

/*
  General solver class:
    T  - numerical type for time
    F  - numerical type for states
    Derivative  - derivative in ODE (right side of the equation)
    Method  - explicit method of integration with given order
*/

template <class T, class F, class Method, class Perturbation>
class TPertInt {

  int dim;    // dimension of phase space

  Method *method;  // pointer to a method

  F *s1, *s2;     // states needed in the steping

  public:

  // setup of adaptive algorithm
  T TINY, SUPERTINY;

  int PERIOD_OF_CHECKING, PERIOD_OF_NOTCHECKING;

  TPertInt(Perturbation *p) {

    dim = p->get_dim();

    method = new Method(p);
    s1 = new F [dim];
    s2 = new F [dim];

    TINY = std::numeric_limits<F>::epsilon();
    SUPERTINY = TINY/1000;

    PERIOD_OF_CHECKING = 5;
    PERIOD_OF_NOTCHECKING = 20;
  }


  ~TPertInt (){
    delete method;
    delete [] s1;
    delete [] s2;
  }

  /*
    Integrator of the ODE dx/dt = F(x). Making a controled step over the
    length dt with suggested integration step h

    input:
      eps[2] -- abs. (accuracy) and relative (precision) permitted local error
      dt -- time step
      h -- suggested integration step
      t -- current time at which input state is given
      state -- position of the ODE at time t

    output:
      h -- last integration step
      t -- new time, t+dt
      state -- position of the ODE at t+dt

    return:
      number of steps done to integrate between [t,t+dt]
  */

  int step(T eps[2], const T &dt, T & h, T & t, F *state);

   /*
    Calculating the trajectory sampled with time step dt.

    input:
      eps[2] -- abs. (accuracy) and relative (precision) permitted local error
      dt -- time step
      t -- current time at which input state is given
      trajectory[dim] -- initial position x(t)

    output:
      trajectory[dim*(n+1)] - whole trajectory
                            = {x(t), x(t+dt), .., x(t+n*dt)}

      steps - number of steps done to integrate between [t+dt*i,t+dt*(i+1)]
              for i = 0, .., n-1
  */

  void traj(T eps[2], const T &dt, T t, const int &n, F *trajectory, int *steps=0);
};

template <class T, class F, class Method, class Perturbation>
int TPertInt<T, F, Method, Perturbation>::
    step(T eps[2], const T &dt, T & h, T & t, F *state) {

  int i, checked = 0, not_checked = 0,
      a = int(1) << method->order, steps = 0;

  bool check = true;

  T A, D, E, tmp, h_out = h, t_end = t + dt;

  if (h <= 0)
    h_out = h = dt/4;
  else  {
    h_out = h;
    if (t + h > t_end) h = t_end - t;
  }

  do {
    cp(s1, state, dim);
    method -> step(h, t, s1);

    if (check) {

      cp(s2, state, dim);
      method -> step(h/2, t, s2, method -> k[0]);
      method -> step(h/2, t + h/2, s2);

      // max-norm of ds[i] -- estimate of abs. local error
      E = 0;
      for (i = 0; i < dim; ++i)
        if ((tmp = std::abs(s1[i] - s2[i])) > E) E = tmp;

      E *= a/(a - 1);

      // max-norm of s1[i], s2[i]
      A = 0;
      for (i = 0; i < dim; ++i)
        if ((tmp = std::max(std::abs(s1[i]), std::abs(s2[i]))) > A) A = tmp;

      // heuristic expression of permitted local error
      D = eps[0] + eps[1]*A;

      #if defined(DEBUG)
      std::cerr << t << '\t' << h << '\t' << E << '\t' << checked << '\n';
      #endif

      if (E > D) {
        h_out = (h *= 0.9*std::exp(std::log(D/E)/(method->order + 1)));
        checked = 0;
      } else {
        t += h;
        for (i = 0; i < dim; ++i) state[i] = (s2[i]*a - s1[i])/(a - 1);

        if (std::abs(t - t_end) < TINY*std::max(T(1), std::abs(t_end))) break;

        // is the step is too small but larger then zero
        if (SUPERTINY < E && E < D/2) {
          h_out = (h *= 0.9*std::exp(std::log(D/E)/(method->order + 1)));
          checked = 0;
        } else if (++checked > PERIOD_OF_CHECKING) {
          check = false;
          not_checked = 0;
        }

        if (t + h > t_end) h = t_end - t;
      }

    } else {
      #if defined(DEBUG)
      std::cerr << t << '\n';
      #endif

      t += h;
      cp(state, s1, dim);

      if (std::abs(t - t_end) < TINY*std::max(T(1), std::abs(t_end))) break;

      if (++not_checked > PERIOD_OF_NOTCHECKING) {
        check = true;
        checked = 0;
      }

      if (t + h > t_end) h = t_end - t;
    }

    ++steps;
  } while (true);

  h = h_out;
  return steps;
}

template <class T, class F, class Method, class Perturbation>
void TPertInt<T, F, Method, Perturbation>::
  traj(T eps[2], const T &dt, T t, const int &n, F *trajectory, int *steps) {

  T h = -1;

  for (int i = 0; i < n; ++i){
    cp(trajectory + dim, trajectory, dim);
    trajectory += dim;
    if (steps)
      steps[i] = step(eps, dt, h, t, trajectory);
    else
      step(eps, dt, h, t, trajectory);
  }
}


#define DERIVATIVE p->getQPdot


/*
   Trapezoidal (aka. two-stage Adams–Moulton) method class, p = 2
   Local error O(h^{p+1})
*/

template <class T, class F, class Perturbation>
class Trapezoidal {

  int dim;

  Perturbation *p;

  T epsA, epsR; // precision and accuracy

  int Nmax;     // max. steps in solving the implicit step

  public:

  F **k, *state_;

  static int order;

  Trapezoidal (Perturbation *p) : p(p) {

    dim = p -> get_dim();

    k = matrix<F> (3, dim);

    state_= k[2];

    /*
    k = new F* [2];
    for (int i = 0; i < 2; ++i) k[i] = new F [dim];
    state_ = new F [dim];
    */
    epsR = 10*std::numeric_limits<T>::epsilon();
    epsA = 1000*std::numeric_limits<T>::min();

    Nmax = 1000;
  }

  ~Trapezoidal(){
    free_matrix(k);
    /*
    for (int i = 0; i < 2; ++i) delete [] k[i];
    delete [] state_;
    delete [] k;*/
  }

  void get_derivative(F *dy){cp(dy, k[0], dim);}
  void calc_derivative(const F &t, F *y, F *dy){DERIVATIVE(t, y, dy);}

  // y_{n+1} = y_n + h( f(t_n, y_n) + f(t_{n+1}, y_{n+1}))/2
  void step(const T &h, const T &t, F *state, F *k0 = 0){

    int i;

    if (k0 == 0)
      DERIVATIVE(t, state, k[0]);
    else if (k0 != k[0])
      cp(k[0], k0, dim);

    // y_ = y_{n} + 1/2 h f(y_n)
    for (i = 0; i < dim; ++i) state_[i] = state[i] + h*k[0][i]/2;

    // Euler step: 1. approximation
    for (i = 0; i < dim; ++i) state[i] += h*k[0][i];

    // Solving
    //    y_{n+1} = y_ + h f(t_{n+1}, y_{n+1})/2
    // via simple forward iteration.
    T q;

    int nr = 0, j;

    do {
      DERIVATIVE(t + h, state, k[1]);

      j = 0;

      for (i = 0; i < dim; ++i) {
        q = state_[i] + h*k[1][i]/2;
        if (std::abs(q - state[i]) < epsA + epsR*std::max(std::abs(q), std::abs(state[i]))) ++j;
        state[i] = q;
      }

      if (++nr > Nmax){
        std::cerr << "Trapezoidal: Too many iterations\n";
        break;
      }

    } while (j < dim);
  }

  void step_(const T &h, const T &t, F *state, F *k0 = 0){step(h, t, state, k0);}
};

template <class T, class F, class Perturbation>
int Trapezoidal<T,F,Perturbation>::order = 2;

/*
   RK2 method class, p =2
*/

template <class T, class F, class Perturbation>
class RK2 {

  int dim;

  Perturbation *p;

  public:

  F **k, *state_;

  static int order;

  RK2 (Perturbation *p) : p(p) {
    dim = p -> get_dim();

    k = matrix <F> (3, dim);

    state_ = k[2];
    /*
    k = new F* [2];
    for (int i = 0; i < 2; ++i) k[i] = new F [dim];
    state_ = new F [dim];*/
  }

  ~RK2(){
    free_matrix(k);
    /*for (int i = 0; i < 2; ++i) delete [] k[i];
    delete [] state_;
    delete [] k;*/
  }

  void get_derivative(F *dy){cp(dy, k[0], dim);}
  void calc_derivative(const F &t, F *y, F *dy){DERIVATIVE(t, y, dy);}

  void step(const T &h, const T &t, F *state, F *k0 = 0){

    int i;

    if (k0 == 0)
      DERIVATIVE(t, state, k[0]);
    else if (k0 != k[0])
      cp(k[0], k0, dim);

    for (i = 0; i < dim; ++i) state_[i] = state[i] + h*k[0][i];

    DERIVATIVE(t + h, state_, k[1]);

    for (i = 0; i < dim; ++i)
      state[i] += h*(k[0][i] + k[1][i])/2;
  }

  void step_(const T &h, const T &t, F *state, F *k0 = 0){step(h, t, state, k0);}
};

template <class T, class F, class Perturbation>
int RK2<T,F,Perturbation>::order = 2;

/*
  Modified Euler or midpoint method  class, p =2
*/

template <class T, class F, class Perturbation>
class ModifiedEuler {

  int dim;

  Perturbation *p;

  public:

  F **k, *state_;

  static int order;

  ModifiedEuler (Perturbation *p) : p(p) {

    dim = p -> get_dim();

    k = matrix <F> (3, dim);

    state_ = k[2];

    /*
    k = new F* [2];
    for (int i = 0; i < 2; ++i) k[i] = new F [dim];
    state_ = new F [dim];
    */
  }

  ~ModifiedEuler(){
    free_matrix(k);
    /*
    for (int i = 0; i < 2; ++i) delete [] k[i];
    delete [] state_;
    delete [] k;
    */
  }

  void get_derivative(F *dy){cp(dy, k[0], dim);}
  void calc_derivative(const F &t, F *y, F *dy){DERIVATIVE(t, y, dy);}

  void step(const T &h, const T &t, F *state, F *k0 = 0){

    int i;

    if (k0 == 0)
      DERIVATIVE(t, state, k[0]);
    else if (k0 != k[0])
      cp(k[0], k0, dim);

    for (i = 0; i < dim; ++i) state_[i] = state[i] + h*k[0][i]/2;

    DERIVATIVE(t + h/2, state_, k[1]);

    for (i = 0; i < dim; ++i) state[i] += h*k[1][i];
  }

  void step_(const T &h, const T &t, F *state, F *k0 = 0){step(h, t, state, k0);}
};

template <class T, class F, class Perturbation>
int ModifiedEuler<T,F,Perturbation>::order = 2;



/*
   RK4 method class
*/


template <class T, class F, class Perturbation>
class RK4 {

  int dim;

  Perturbation *p;

  public:

  F **k, *state_;

  static int order;

  RK4 (Perturbation *p) : p(p) {

    dim = p -> get_dim();

    k = matrix <F> (5, dim);

    state_ = k[4];

    /*
    k = new F* [4];
    for (int i = 0; i < 4; ++i) k[i] = new F [dim];
    state_ = new F [dim];
    */
  }

  ~RK4(){
    free_matrix(k);
    /*
    for (int i = 0; i < 4; ++i) delete [] k[i];
    delete [] state_;
    delete [] k;
    */
  }

  void get_derivative(F *dy){cp(dy, k[0], dim);}
  void calc_derivative(const F &t, F *y, F *dy){DERIVATIVE(t, y, dy);}

  void step(const T &h, const T &t, F *state, F *k0 = 0){

    int i;

    if (k0 == 0)
      DERIVATIVE(t, state, k[0]);
    else if (k0 != k[0])
      cp(k[0], k0, dim);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + h*k[0][i]/2;

    DERIVATIVE(t + h/2, state_, k[1]);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + h*k[1][i]/2;

    DERIVATIVE(t + h/2, state_, k[2]);

    for (i = 0; i < dim; ++i)
      state_[i] = state[i] + h*k[2][i];

    DERIVATIVE(t + h, state_, k[3]);

    for (i = 0; i < dim; ++i)
      state[i] += h*(k[0][i] + k[3][i] + 2*(k[1][i] + k[2][i]))/6;
  }

  void step_(const T &h, const T &t, F *state, F *k0 = 0){step(h, t, state, k0);}
};

template <class T, class F, class Perturbation>
int RK4<T,F,Perturbation>::order = 4;

/*
 Hairer order (p=8), Peter Stone's version

 Ref:
 * High order embedded Runge-Kutta formulae, by P.J.Prince and J.R.Dormand,
Journal of Computational and Applied Mathematics, vol. 7, 1981, pages 67-75.
 * http://www.peterstone.name/Maplepgs/RKcoeff_8.html

*/
template <class T, class F, class Perturbation>
class Hairer8 {

  int dim;

  Perturbation *p;

  T **a, *b, *c;

  public:

  F **k, *state_;

  static int order;

  Hairer8 (Perturbation *p) : p(p) {
    int i;

    dim = p -> get_dim();

    k = matrix<F> (14, dim);

    state_ = k[13];

    /*
    k = new F* [13];
    for (i = 0; i < 13; ++i) k[i] = new F [dim];
    state_ = new F [dim];
    */

    a = new T* [12];

    for (i = 0; i < 12; ++i) a[i] = new T [i+1];

    b = new T [13];

    c = new T [12];

    Hairer8_coef(c, a, b);

  }

  ~Hairer8(){
    int i;

    delete [] b;
    delete [] c;

    for (i = 0; i < 12; i++) delete [] a[i];
    delete [] a;

    free_matrix(k);
    /*
    delete [] state_;
    for (i = 0; i < 13; i++) delete [] k[i];
    delete [] k;
    */

  }

  void get_derivative(F *dy){cp(dy, k[0], dim);}
  void calc_derivative(const F &t, F *y, F *dy){DERIVATIVE(t, y, dy);}

  void step(const T &h, const T &t, F *state, F *k0 = 0){

    int i, j, l;

    F tmp;

    if (k0 == 0)
      DERIVATIVE(t, state, k[0]);
    else if (k0 != k[0])
      cp(k[0], k0, dim);

    for (i = 0; i < 12; ++i) {

      for (l = 0; l < dim; ++l) {

        tmp = 0;
        for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

        state_[l] = state[l] + h*tmp;
      }

      DERIVATIVE(t + c[i]*h, state_, k[i+1]);
    }

    for (l = 0; l < dim; ++l) {
      tmp = 0;
      for (i = 0; i < 13; i++) tmp += b[i]*k[i][l];
      state[l] += h*tmp;
    }
  }

  void step_(const T &h, const T &t, F *state, F *k0 = 0){step(h, t, state, k0);}
};

template <class T, class F, class Perturbation>
int Hairer8 <T,F,Perturbation>::order = 8;

/*
  Hairer order (p=10), Peter Stone's version

  Ref:
  * "A Runge-Kutta Method of Order 10", by E. Hairer, Journal of the Institute of Mathematics and its Applications (1975) 16, pages 35 to 55.
  *  http://www.peterstone.name/Maplepgs/RKcoeff_10.html
*/
template <class T, class F, class Perturbation>
class Hairer10 {

  int dim;

  Perturbation *p;

  T **a, *b, *c;

  public:

  F **k, *state_;

  static int order;

  Hairer10 (Perturbation *p) : p(p) {

    dim = p -> get_dim();

    k = matrix<F>(18, dim);

    state_ = k[17];

    /*
    k = new F* [17];
    for (i = 0; i < 17; ++i) k[i] = new F [dim];
    state_ = new F [dim];
    */

    c = new T [16];

    a = new T* [16];

    for (int i = 0; i < 16; ++i) a[i] = new T [i+1];

    b = new T [17];

    Hairer10_coef(c, a, b);
  }

  ~Hairer10(){
    delete [] b;
    delete [] c;

    for (int i = 0; i < 16; i++) delete [] a[i];
    delete [] a;

    free_matrix(k);
    /*
    delete [] state_;
    for (i = 0; i < 17; i++) delete [] k[i];
    delete [] k;
    */

  }

  void get_derivative(F *dy){cp(dy, k[0], dim);}
  void calc_derivative(const F &t, F *y, F *dy){DERIVATIVE(t, y, dy);}

  void step(const T &h, const T &t, F *state, F* k0 = 0){

    int i, j, l;

    F tmp;

    if (k0 == 0)
      DERIVATIVE(t, state, k[0]);
    else if (k0 != k[0])
      cp(k[0], k0, dim);

    for (i = 0; i < 16; ++i) {

      for (l = 0; l < dim; ++l) {

        tmp = 0;
        for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

        state_[l] = state[l] + h*tmp;
      }

      DERIVATIVE(t + c[i]*h, state_, k[i+1]);
    }

    for (l = 0; l < dim; ++l) {
      tmp = 0;
      for (i = 0; i < 17; i++) tmp += b[i]*k[i][l];
      state[l] += h*tmp;
    }
  }
  void step_(const T &h, const T &t, F *state, F *k0 = 0){step(h, t, state, k0);}
};

template <class T, class F, class Perturbation>
int Hairer10 <T,F,Perturbation>::order = 10;

#undef  DERIVATIVE

#endif // #if defined(EMBEDDED)


