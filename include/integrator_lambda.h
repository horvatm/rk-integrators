#pragma once

/*
  Library for adaptive step evolution motions.

    dx/dt = derivative(t, x)

    x - state
    t - time

  Library provides

    class TIntegrator

  in two namespaces

  Thredsave and NonTreadsave:
    Embedded:
      RKM4  -- Marson method
      Hairer8 -- Hairer (p=8) method

    NonEmbedded:
      RK4 -- standard RK4
      Hairer8 -- Hairer (p=8) method
      Hairer10 -- Hairer (p=10) method

  For heuristic step-size control we use

  http://www.gnu.org/software/gsl/manual/html_node/Adaptive-Step_002dsize-Control.html

  where we set a_y = 1 and a_dydt=0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Derivative of the form  derivative(t, y, dydt) is given by a lambda
    function and functins are not thread save
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Author: Martin Horvat, September 2021
*/

#include <cmath>
#include <limits>
#include <complex>

// C-style matrix allocation and deallocation
#include "matrix.h"

// RK coefficients
#include "RKcoeff8d_1.h"
#include "RKcoeff10b_1.h"

// enabling complex*int, int*complex, complex/int, int/complex operations
template< typename T, typename S > inline
typename std::enable_if< !std::is_same<T,S>::value, std::complex<T> >::type
operator* ( const std::complex<T>& c, S n ) { return c * T(n) ; }

template< typename T, typename S > inline
typename std::enable_if< !std::is_same<T,S>::value, std::complex<T> >::type
operator* ( S n, const std::complex<T>& c ) { return T(n) * c ; }


template< typename T, typename S > inline
typename std::enable_if< !std::is_same<T,S>::value, std::complex<T> >::type
operator/ ( const std::complex<T>& c, S n ) { return c / T(n) ; }

template< typename T, typename S > inline
typename std::enable_if< !std::is_same<T,S>::value, std::complex<T> >::type
operator/ ( S n, const std::complex<T>& c ) { return T(n) / c ; }

namespace ThreadSave {

  namespace Embedded {

    /*
      General solver class:
        T  - numerical type for time
        F  - numerical type for states
        TDerivative  - derivative in ODE (right side of the equation)
        TMethod  - explicit method of integration with given order
    */

    template <class T, class F, class TMethod >
    class TIntegrator {

      int dim;                      // dimension of phase space

      TMethod *method;               // pointer to a method

      public:

      // setup of adaptive algorithm
      T TINY, SUPERTINY;

      int PERIOD_OF_CHECKING, PERIOD_OF_NOTCHECKING;

      TIntegrator(int dim): dim(dim) {

        method = new TMethod(dim);

        TINY = std::numeric_limits<T>::epsilon();
        SUPERTINY = TINY/1000;

        PERIOD_OF_CHECKING = 5;
        PERIOD_OF_NOTCHECKING = 20;
      }

      ~TIntegrator (){
        delete method;
      }

      /*
        Integrator of the ODE dx/dt = F(t, x). Making a controlled step over the
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

      template <typename Derivative>
      int step(T eps[2], const T &dt, T & h, T & t, F *state, Derivative deriv){

        int i, checked = 0, not_checked = 0, steps = 0;

        bool check = true;

        F  *s = new F [2*dim],
           *ds = s + dim;

        T A, E, D, tmp, h_out = h, t_end = t + dt;

        if (h <= 0)
          h_out = h = dt/4;
        else if (t + h > t_end)
          h = t_end - t;

        cp(s, state, dim);

        do {

          method->step_embedded(h, t, s, ds, deriv);

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

        delete [] s;

        return steps;
      }
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
      template <typename Derivative>
      void traj(T eps[2], const T &dt, T t, const int &n, F *trajectory, Derivative deriv, int *steps=0){

        T h = -1;

        for (int i = 0; i < n; ++i){
          cp(trajectory + dim, trajectory, dim);
          trajectory += dim;

          if (steps)
            steps[i] = step(eps, dt, h, t, trajectory, deriv);
          else
            step(eps, dt, h, t, trajectory, deriv);
        }
      }
    };

    /*
       Runge-Kutta-Merson method of 4th order.

       Ref: Bohte Z  Numericne metode  (DMFA, 1991)
    */

    template <class T, class F>
    class RKM4 {

      int dim;

      public:

      const int order = 4;

      RKM4 (int dim): dim(dim) { }

      ~RKM4(){ }

      template <typename D>
      void step_embedded(const T &h, const T &t, F *state, F *dstate, D deriv){

        int i;

        F **k = matrix<F> (6, dim),
          *state_ = k[5];

        deriv(t, state, k[0]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + (k[0][i] *= h)/3;

        deriv(t + h/3, state_, k[1]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + ((k[1][i] *= h) + k[0][i])/T(6);

        deriv(t + h/3, state_, k[2]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + (T(3)*(k[2][i] *= h) + k[0][i])/T(8);

        deriv(t + T(0.5)*h, state_, k[3]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + T(2)*(k[3][i] *= h) - T(1.5)*k[2][i] + T(0.5)*k[0][i];

        deriv(t + h, state_, k[4]);

        if (dstate)
          for (i = 0; i < dim; ++i) {
            state[i] += (k[0][i] + T(4)*k[3][i] + (k[4][i] *= h))/T(6);
            dstate[i] = (T(2)*k[0][i] - T(9)*k[2][i] + T(8)*k[3][i] - k[4][i])/T(30);
          }
        else
          for (i = 0; i < dim; ++i)
            state[i] += (k[0][i] + T(4)*k[3][i] + (k[4][i] *= h))/T(6);

        free_matrix(k);
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;

         F  **k = matrix<F> (6, dim),
            *state_ = k[5];

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + (k[0][i] *= h)/3;

        deriv(t + h/3, state_, k[1]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + ((k[1][i] *= h) + k[0][i])/T(6);

        deriv(t + h/3, state_, k[2]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + (3*(k[2][i] *= h) + k[0][i])/T(8);

        deriv(t + T(0.5)*h, state_, k[3]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + T(2)*(k[3][i] *= h) - T(1.5)*k[2][i] + T(0.5)*k[0][i];

        deriv(t + h, state_, k[4]);

        for (i = 0; i < dim; ++i)
          state[i] += (k[0][i] + T(4)*k[3][i] + (k[4][i] *= h))/T(6);

        free_matrix(k);
      }

    };

    /*
       Hairer method (p=8)

       Ref: P.J. Prince and J.R. Dormand
    */

    template <class T, class F>
    class Hairer8 {

      int dim;

      T **a, *b, *b_, *c;

      public:

      const int order = 8;

      Hairer8 (int dim): dim(dim) {

        a = new T* [12];

        for (int i = 0; i < 12; ++i) a[i] = new T [i+1];

        b = new T [13];

        b_ = new T [13];

        c = new T [12];

        Hairer8_coef(c, a, b, b_);
      }

      ~Hairer8(){

        delete [] b;
        delete [] b_;
        delete [] c;

        for (int i = 0; i < 12; i++) delete [] a[i];
        delete [] a;
      }

      template <typename D>
      void step_embedded(const T &h, const T &t, F *state, F *dstate, D deriv){

        int i, j, l;

        F tmp, tmp_,
          **k = matrix <F> (14, dim),
          *state_ = k[13];

        deriv(t, state, k[0]);

        for (i = 0; i < 12; ++i) {

          for (l = 0; l < dim; ++l) {

            tmp = 0;
            for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

            state_[l] = state[l] + h*tmp;
          }

          deriv(t + c[i]*h, state_, k[i+1]);
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

        free_matrix(k);
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i, j, l;

        F tmp, tmp_,
          **k = matrix <F> (14, dim),
          *state_ = k[13];

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < 12; ++i) {

          for (l = 0; l < dim; ++l) {

            tmp = 0;
            for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

            state_[l] = state[l] + h*tmp;
          }

          deriv(t + c[i]*h, state_, k[i+1]);
        }

        for (l = 0; l < dim; ++l) {
          tmp = 0.0;

          for (i = 0; i < 13; i++) tmp  += b[i]*k[i][l];

          state[l] += h*tmp;
        }
        free_matrix(k);
      }
    };
  }

  namespace NonEmbedded{

  /* ============================================================================
    Integrator of 1. order ODE defined as

      dx/dt = F(t, x)     x in T^d

    using an adaptive step in a combination with some explicit method of
    some order.

    General solver class:
      T  - numerical type for time
      F  - numerical type for states
      Derivative  - derivative in ODE (right side of the equation)
      TMethod  - explicit method of integration with given order
   ============================================================================ */

    template <class T, class F, class TMethod>
    class TIntegrator {

      int dim;    // dimension of phase space

      TMethod *method;  // pointer to a method

      public:

      // setup of adaptive algorithm
      T TINY, SUPERTINY;

      int PERIOD_OF_CHECKING, PERIOD_OF_NOTCHECKING;

      TIntegrator(int dim): dim(dim) {

        method = new TMethod(dim);

        TINY = std::numeric_limits<T>::epsilon();
        SUPERTINY = TINY/1000;

        PERIOD_OF_CHECKING = 5;
        PERIOD_OF_NOTCHECKING = 20;
      }

      ~TIntegrator (){
        delete method;
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
      template <typename Derivative >
      int step(T eps[2], const T &dt, T & h, T & t, F *state, Derivative deriv){

        int i, checked = 0, not_checked = 0,
            a = int(1) << method->order, steps = 0;

        bool check = true;

        T A, D, E, tmp, h_out = h, t_end = t + dt;

        F *s1 = new F [2*dim],
          *s2 = s1 + dim;

        if (h <= 0)
          h_out = h = dt/4;
        else  {
          h_out = h;
          if (t + h > t_end) h = t_end - t;
        }

        do {
          cp(s1, state, dim);
          method -> step(h, t, s1, deriv);

          if (check) {

            cp(s2, state, dim);
            method -> step(T(0.5)*h, t, s2, deriv, method->k[0]);
            method -> step(T(0.5)*h, t + T(0.5)*h, s2, deriv);

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

        delete [] s1;

        return steps;
      }

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
      template <typename Derivative>
      void traj(T eps[2], const T &dt, T t, const int &n, F *trajectory, Derivative deriv, int *steps=0) {

        T h = -1;

        for (int i = 0; i < n; ++i){
          cp(trajectory + dim, trajectory, dim);
          trajectory += dim;
          if (steps)
            steps[i] = step(eps, dt, h, t, trajectory, deriv);
          else
            step(eps, dt, h, t, trajectory, deriv);
        }
      }
    };

    /*
       Trapezoidal (aka. two-stage Adams–Moulton) method class, p = 2
       Local error O(h^{p+1})
    */

    template <class T, class F>
    class Trapezoidal {

      int dim;

      T epsA, epsR; // precision and accuracy

      int Nmax;     // max. steps in solving the implicit step

      public:

      const int order = 2;

      Trapezoidal (int dim): dim(dim) {

        epsR = 10*std::numeric_limits<T>::epsilon();
        epsA = 1000*std::numeric_limits<T>::min();

        Nmax = 1000;
      }

      ~Trapezoidal(){ }

      // y_{n+1} = y_n + h( f(t_n, y_n) + f(t_{n+1}, y_{n+1}))/2
      template <typename D> void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;

        F **k = matrix<F> (3, dim),
          *state_ = k[2];

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        // y_ = y_{n} + 1/2 h f(y_n)
        for (i = 0; i < dim; ++i) state_[i] = state[i] + T(0.5)*h*k[0][i];

        // Euler step: 1. approximation
        for (i = 0; i < dim; ++i) state[i] += h*k[0][i];

        // Solving
        //    y_{n+1} = y_ + h f(t_{n+1}, y_{n+1})/2
        // via simple forward iteration.
        T q;

        int nr = 0, j;

        do {
          deriv(t + h, state, k[1]);

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

        free_matrix(k);
      }


    };

    /*
       RK2 method class, p =2
    */

    template <class T, class F>
    class RK2 {

      int dim;

      public:

      const int order = 2;

      RK2 (int dim): dim(dim) { }

      ~RK2(){ }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;

        F **k = matrix <F> (3, dim),
          *state_ = k[2];

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < dim; ++i) state_[i] = state[i] + h*k[0][i];

        deriv(t + h, state_, k[1]);

        for (i = 0; i < dim; ++i)  state[i] += h*(k[0][i] + k[1][i])/2;

        free_matrix(k);
      }


    };

    /*
      Modified Euler or midpoint method  class, p = 2
    */

    template <class T, class F>
    class ModifiedEuler {

      int dim;

      public:

      const int order = 2;

      ModifiedEuler (int dim): dim(dim) { }

      ~ModifiedEuler(){ }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;

        F **k  = matrix <F> (3, dim),
          *state_ = k[2];

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < dim; ++i) state_[i] = state[i] + T(0.5)*h*k[0][i];

        deriv(t + T(0.5)*h, state_, k[1]);

        for (i = 0; i < dim; ++i) state[i] += h*k[1][i];

        free_matrix(k);
      }


    };

    /*
       RK4 method class
    */


    template <class T, class F>
    class RK4 {

      int dim;

      public:

      const int order = 4;

      RK4 (int dim): dim(dim) { }

      ~RK4(){ }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;

        F **k = matrix <F> (5, dim),
          *state_ = k[4];

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + T(0.5)*h*k[0][i];

        deriv(t + T(0.5)*h, state_, k[1]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + T(0.5)*h*k[1][i];

        deriv(t + T(0.5)*h, state_, k[2]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + h*k[2][i];

        deriv(t + h, state_, k[3]);

        for (i = 0; i < dim; ++i)
          state[i] += h*(k[0][i] + k[3][i] + T(2)*(k[1][i] + k[2][i]))/T(6);

        free_matrix(k);
      }


    };

    /*
     Hairer order (p=8), Peter Stone's version

     Ref:
     * High order embedded Runge-Kutta formulae, by P.J.Prince and J.R.Dormand,
    Journal of Computational and Applied Mathematics, vol. 7, 1981, pages 67-75.
     * http://www.peterstone.name/Maplepgs/RKcoeff_8.html

    */
    template <class T, class F>
    class Hairer8 {

      int dim;

      T **a, *b, *c;

      public:

      int order = 8;

      Hairer8 (int dim): dim(dim) {

        a = new T* [12];

        for (int i = 0; i < 12; ++i) a[i] = new T [i+1];

        b = new T [13];

        c = new T [12];

        Hairer8_coef(c, a, b);

      }

      ~Hairer8(){

        delete [] b;
        delete [] c;

        for (int i = 0; i < 12; i++) delete [] a[i];
        delete [] a;
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i, j, l;

        F tmp,
          **k = matrix<F> (14, dim),
          *state_ = k[13];

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < 12; ++i) {

          for (l = 0; l < dim; ++l) {

            tmp = 0;
            for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

            state_[l] = state[l] + h*tmp;
          }

          deriv(t + c[i]*h, state_, k[i+1]);
        }

        for (l = 0; l < dim; ++l) {
          tmp = 0;
          for (i = 0; i < 13; i++) tmp += b[i]*k[i][l];
          state[l] += h*tmp;
        }

        free_matrix(k);
      }


    };

    /*
      Hairer order (p=10), Peter Stone's version

      Ref:
      * "A Runge-Kutta TMethod of Order 10", by E. Hairer, Journal of the Institute of
        Mathematics and its Applications (1975) 16, pages 35 to 55.
      *  http://www.peterstone.name/Maplepgs/RKcoeff_10.html
    */
    template <class T, class F>
    class Hairer10 {

      int dim;

      T **a, *b, *c;

      public:

      int order = 10;

      Hairer10 (int dim): dim(dim) {

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
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F* k0 = 0){

        int i, j, l;

        F tmp,
          **k = matrix <F> (18, dim),
          *state_ = k[17];

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < 16; ++i) {

          for (l = 0; l < dim; ++l) {

            tmp = 0;
            for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

            state_[l] = state[l] + h*tmp;
          }

          deriv(t + c[i]*h, state_, k[i+1]);
        }

        for (l = 0; l < dim; ++l) {
          tmp = 0;
          for (i = 0; i < 17; i++) tmp += b[i]*k[i][l];
          state[l] += h*tmp;
        }

        free_matrix(k);
      }

    };

  } // namespace ThreadSave::NonEmbedded
}  // namespace ThreadSave

namespace NonThreadSave {

  namespace Embedded {

    /* ============================================================================
      Integrator of 1. order ODE defined as

        dx/dt = F(t, x)     x in T^d

      using an adaptive step in a combination with some explicit method of
      some order with embedded error estimate.

      General solver class:
        T  - numerical type for time
        F  - numerical type for states
        TDerivative  - derivative in ODE (right side of the equation)
        TMethod  - explicit method of integration with given order
     ============================================================================ */


    template <class T, class F, class TMethod>
    class TIntegrator {

      int dim;                      // dimension of phase space

      TMethod *method;               // pointer to a method

      F *s, *ds;                     // states needed in the steping

      public:

      // setup of adaptive algorithm
      T TINY, SUPERTINY;

      int PERIOD_OF_CHECKING, PERIOD_OF_NOTCHECKING;

      TIntegrator(int dim): dim(dim) {

        method = new TMethod(dim);
        s = new F [dim];
        ds = new F [dim];

        TINY = std::numeric_limits<T>::epsilon();
        SUPERTINY = TINY/1000;

        PERIOD_OF_CHECKING = 5;
        PERIOD_OF_NOTCHECKING = 20;
      }


      ~TIntegrator (){
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
      template <typename Derivative>
      int step(T eps[2], const T &dt, T & h, T & t, F *state, Derivative deriv) {

        int i, checked = 0, not_checked = 0, steps = 0;

        bool check = true;

        T A, E, D, tmp, h_out = h, t_end = t + dt;

        if (h <= 0)
          h_out = h = dt/4;
        else if (t + h > t_end)
          h = t_end - t;

        cp(s, state, dim);

        do {

          method->step_embedded(h, t, s, ds, deriv);

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
      template <typename Derivative>
      void traj(T eps[2], const T &dt, T t, const int &n, F *trajectory, Derivative deriv, int *steps=0){

        T h = -1;

        for (int i = 0; i < n; ++i){
          cp(trajectory + dim, trajectory, dim);
          trajectory += dim;
          if (steps)
            steps[i] = step(eps, dt, h, t, trajectory, deriv);
          else
            step(eps, dt, h, t, trajectory, deriv);
        }
      }
    };

    /*
       Runge-Kutta-Merson method of 4th order.

       Ref: Bohte Z  Numericne metode  (DMFA, 1991)
    */

    template <class T, class F>
    class RKM4 {

      int dim;

      public:

      F **k, *state_;

      const int order = 4;

      RKM4 (int dim): dim(dim)  {

        k = matrix<F> (6, dim);
        state_ = k[5];
      }

      ~RKM4(){
        free_matrix(k);
      }

      template <typename D>
      void step_embedded(const T &h, const T &t, F *state, F *dstate, D deriv){

        int i;

        deriv(t, state, k[0]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + (k[0][i] *= h)/3;

        deriv(t + h/3, state_, k[1]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + ((k[1][i] *= h) + k[0][i])/6;

        deriv(t + h/3, state_, k[2]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + (T(3)*(k[2][i] *=h) + k[0][i])/T(8);

        deriv(t + T(0.5)*h, state_, k[3]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + T(2)*(k[3][i] *= h) - T(1.5)*k[2][i] + T(0.5)*k[0][i];

        deriv(t + h, state_, k[4]);

        if (dstate)
          for (i = 0; i < dim; ++i) {
            state[i] += (k[0][i] + T(4)*k[3][i] + (k[4][i] *= h))/6;
            dstate[i] = (T(2)*k[0][i] - T(9)*k[2][i] + T(8)*k[3][i] - k[4][i])/T(30);
          }
        else
          for (i = 0; i < dim; ++i)
            state[i] += (k[0][i] + T(4)*k[3][i] + (k[4][i] *= h))/T(6);
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + (k[0][i] *= h)/T(3);

        deriv(t + h/3, state_, k[1]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + ((k[1][i] *= h) + k[0][i])/T(6);

        deriv(t + h/3, state_, k[2]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + (T(3)*(k[2][i] *=h) + k[0][i])/T(8);

        deriv(t + T(0.5)*h, state_, k[3]);

        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + T(2)*(k[3][i] *= h) - T(1.5)*k[2][i] + T(0.5)*k[0][i];

        deriv(t + h, state_, k[4]);

        for (i = 0; i < dim; ++i)
          state[i] += (k[0][i] + T(4)*k[3][i] + (k[4][i] *= h))/T(6);
      }

    };

    /*
       Hairer method (p=8)
       Ref: P.J. Prince and J.R. Dormand
    */

    template <class T, class F>
    class Hairer8 {

      int dim;

      T **a, *b, *b_, *c;

      public:

      F **k, *state_;

      const int order = 8;

      Hairer8 (int dim): dim(dim)  {

        k = matrix <F> (14, dim);

        state_ = k[13];

        a = new T* [12];

        for (int i = 0; i < 12; ++i) a[i] = new T [i+1];

        b = new T [13];

        b_ = new T [13];

        c = new T [12];

        Hairer8_coef(c, a, b, b_);
      }

      ~Hairer8(){

        delete [] b;
        delete [] b_;
        delete [] c;

        for (int i = 0; i < 12; i++) delete [] a[i];
        delete [] a;

        free_matrix(k);
      }

      template <typename D>
      void step_embedded(const T &h, const T &t, F *state, F *dstate, D deriv){

        int i, j, l;

        F tmp, tmp_;

        deriv(t, state, k[0]);

        for (i = 0; i < 12; ++i) {

          for (l = 0; l < dim; ++l) {

            tmp = 0;
            for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

            state_[l] = state[l] + h*tmp;
          }

          deriv(t + c[i]*h, state_, k[i+1]);
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

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i, j, l;

        F tmp, tmp_;

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < 12; ++i) {

          for (l = 0; l < dim; ++l) {

            tmp = 0;
            for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

            state_[l] = state[l] + h*tmp;
          }

          deriv(t + c[i]*h, state_, k[i+1]);
        }

        for (l = 0; l < dim; ++l) {
          tmp = 0.0;

          for (i = 0; i < 13; i++) tmp  += b[i]*k[i][l];

          state[l] += h*tmp;
        }
      }
    };

  } // namespace ThreadSave::Embedded

  namespace NonEmbedded {

    /* ============================================================================
      Integrator of 1. order ODE defined as

        dx/dt = F(t, x)     x in T^d

      using an adaptive step in a combination with some explicit method of
      some order.

      General solver class:
        T  - numerical type for time
        F  - numerical type for states
        TDerivative  - derivative in ODE (right side of the equation)
        TMethod  - explicit method of integration with given order
     ============================================================================ */


    template <class T, class F, class TMethod>
    class TIntegrator {

      int dim;    // dimension of phase space

      TMethod *method;  // pointer to a method

      F *s1, *s2;     // states needed in the steping

      public:

      // setup of adaptive algorithm
      T TINY, SUPERTINY;

      int PERIOD_OF_CHECKING, PERIOD_OF_NOTCHECKING;

      TIntegrator(int dim): dim(dim)  {

        method = new TMethod(dim);
        s1 = new F [dim];
        s2 = new F [dim];

        TINY = std::numeric_limits<T>::epsilon();
        SUPERTINY = TINY/1000;

        PERIOD_OF_CHECKING = 5;
        PERIOD_OF_NOTCHECKING = 20;
      }

      ~TIntegrator (){
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

      template <typename Derivative>
      int step(T eps[2], const T &dt, T & h, T & t, F *state, Derivative deriv){

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
          method -> step(h, t, s1, deriv);

          if (check) {

            cp(s2, state, dim);
            method -> step(T(0.5)*h, t, s2, deriv, method->k[0]);
            method -> step(T(0.5)*h, t + T(0.5)*h, s2, deriv);

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
              for (i = 0; i < dim; ++i) state[i] = (s2[i]*T(a) - s1[i])/T(a - 1);

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

      template <typename Derivative>
      void traj(T eps[2], const T &dt, T t, const int &n, F *trajectory, Derivative deriv, int *steps=0){

        T h = -1;

        for (int i = 0; i < n; ++i){
          cp(trajectory + dim, trajectory, dim);
          trajectory += dim;
          if (steps)
            steps[i] = step(eps, dt, h, t, trajectory, deriv);
          else
            step(eps, dt, h, t, trajectory, deriv);
        }
      }
    };

    /*
       Trapezoidal (aka. two-stage Adams–Moulton) method class, p = 2
       Local error O(h^{p+1})
    */

    template <class T, class F>
    class Trapezoidal {

      int dim;

      T epsA, epsR; // precision and accuracy

      int Nmax;     // max. steps in solving the implicit step

      public:

      F **k, *state_;

      const int order = 2;

      Trapezoidal (int dim): dim(dim) {

        k = matrix<F> (3, dim);

        state_= k[2];

        epsR = 10*std::numeric_limits<T>::epsilon();
        epsA = 1000*std::numeric_limits<T>::min();

        Nmax = 1000;
      }

      ~Trapezoidal(){
        free_matrix(k);
      }

      // y_{n+1} = y_n + h( f(t_n, y_n) + f(t_{n+1}, y_{n+1}))/2
      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        // y_ = y_{n} + 1/2 h f(y_n)
        for (i = 0; i < dim; ++i) state_[i] = state[i] + T(0.5)*h*k[0][i];

        // Euler step: 1. approximation
        for (i = 0; i < dim; ++i) state[i] += h*k[0][i];

        // Solving
        //    y_{n+1} = y_ + h f(t_{n+1}, y_{n+1})/2
        // via simple forward iteration.
        T q;

        int nr = 0, j;

        do {
          deriv(t + h, state, k[1]);

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
    };

    /*
       RK2 method class, p =2
    */

    template <class T, class F>
    class RK2 {

      int dim;

      public:

      F **k, *state_;

      const int order = 2;

      RK2 (int dim): dim(dim)  {

        k = matrix <F> (3, dim);

        state_ = k[2];
      }

      ~RK2(){
        free_matrix(k);
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < dim; ++i) state_[i] = state[i] + h*k[0][i];

        deriv(t + h, state_, k[1]);

        for (i = 0; i < dim; ++i)
          state[i] += T(0.5)*h*(k[0][i] + k[1][i]);
      }
    };

    /*
      Modified Euler or midpoint method  class, p =2
    */

    template <class T, class F>
    class ModifiedEuler {

      int dim;

      public:

      F **k, *state_;

      const int order = 2;

      ModifiedEuler (int dim): dim(dim)  {

        k = matrix <F> (3, dim);

        state_ = k[2];
      }

      ~ModifiedEuler(){
        free_matrix(k);
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < dim; ++i) state_[i] = state[i] + T(0.5)*h*k[0][i];

        deriv(t + T(0.5)*h, state_, k[1]);

        for (i = 0; i < dim; ++i) state[i] += h*k[1][i];
      }
    };

    /*
       RK4 method class
    */
    template <class T, class F>
    class RK4 {

      int dim;

      public:

      F **k, *state_;

      const int order = 4;

      RK4 (int dim): dim(dim) {

        k = matrix <F> (5, dim);

        state_ = k[4];
      }

      ~RK4(){
        free_matrix(k);
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < dim; ++i) state_[i] = state[i] + T(0.5)*h*k[0][i];
        deriv(t + T(0.5)*h, state_, k[1]);

        for (i = 0; i < dim; ++i)  state_[i] = state[i] + T(0.5)*h*k[1][i];

        deriv(t + T(0.5)*h, state_, k[2]);

        for (i = 0; i < dim; ++i) state_[i] = state[i] + h*k[2][i];

        deriv(t + h, state_, k[3]);

        for (i = 0; i < dim; ++i)
          state[i] += h*(k[0][i] + k[3][i] + T(2)*(k[1][i] + k[2][i]))/T(6);

      }
    };

    /*
       Runge-Kutta-Butcher method - Sixth order method

       The 7 stage 6th order Runge-Kutta method for solving a differential
       equation y'(x) = f(x,y) with initial condition y = c when x = x0
       evaluates f(x,y) six times per step. For step i+1,
        y[i+1] = y[i] + (11 k1 + 81 k3 + 81 k4 - 32 k5 - 32 k6 + 11 k7) / 120
       where,
         k1 = h * f( x[i], y[i] ),
         k2 = h * f( x[i]+h/3, y[i]+k1/3 ),
         k3 = h * f( x[i]+2h/3, y[i] + 2 k2 / 3 ),
         k4 = h * f( x[i]+h/3, y[i] + ( k1 + 4 k2 - k3 ) / 12 ),
         k5 = h * f( x[i]+h/2, y[i] + (-k1 + 18 k2 - 3 k3 -6 k4)/16 ),
         k6 = h * f( x[i]+h/2, y[i] + (9 k2 - 3 k3 - 6 k4 + 4 k5)/8 ),
         k7 = h * f( x[i]+h, y[i] + (9 k1 - 36 k2 + 63 k3 + 72 k4 - 64 k5)/44 ),

       Ref:
        * Butcher, J. (1963). On the integration processes of A. Huťa. Journal
          of the Australian Mathematical Society, 3(2), 202-206.
          doi:10.1017/S1446788700027944, p. 192.

        * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_butcher.html
        * http://www.mymathlib.com/c_source/diffeq/runge_kutta/runge_kutta_butcher.c

    */
    template <class T, class F>
    class RKButcher {

      int dim;

      public:

      F **k, *state_;

      const int order = 6;

      RKButcher (int dim): dim(dim) {

        k = matrix <F> (8, dim);

        state_ = k[7];
      }

      ~RKButcher(){
        free_matrix(k);
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i;
        // k1 = h * f( x[i], y[i] )

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        //  k2 = h * f( x[i]+h/3, y[i]+k1/3 )
        for (i = 0; i < dim; ++i) state_[i] = state[i] + h*k[0][i]/T(3);
        deriv(t + h/T(3), state_, k[1]);

        // k3 = h * f( x[i]+2h/3, y[i] + 2 k2 / 3 ),
        for (i = 0; i < dim; ++i) state_[i] = state[i] + T(2)*h*k[1][i]/T(3);
        deriv(t + T(2)*h/T(3), state_, k[2]);

        //k4 = h * f( x[i] + h/3, y[i] + ( k1 + 4 k2 - k3 ) / 12 ),
        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + h*(k[0][i] + T(4)*k[1][i] - k[2][i])/T(12);
        deriv(t + h/T(3), state_, k[3]);

        // k5 = h * f( x[i]+h/2, y[i] + (-k1 + 18 k2 - 3 k3 -6 k4)/16 ),
        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + h*(-k[0][i] + T(18)*k[1][i] - T(3)*k[2][i] - T(6)*k[3][i])/T(16);
        deriv(t + h/T(2), state_, k[4]);

        // k6 = h * f( x[i]+h/2, y[i] + (9 k2 - 3 k3 - 6 k4 + 4 k5)/8 ),
        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + h*(T(9)*k[1][i] - T(3)*k[2][i] - T(6)*k[3][i] + T(4)*k[4][i])/T(8);
        deriv(t + h/T(2), state_, k[5]);

        // k7 = h * f( x[i]+h, y[i] + (9 k1 - 36 k2 + 63 k3 + 72 k4 - 64 k5)/44 ),
        for (i = 0; i < dim; ++i)
          state_[i] = state[i] + h*(T(9)*k[0][i] - T(36)*k[1][i] + T(63)*k[2][i] + T(72)*k[3][i] - T(64)*k[4][i])/T(44);
        deriv(t + h, state_, k[6]);


        // y[i+1] = y[i] + (11 k1 + 81 k3 + 81 k4 - 32 k5 - 32 k6 + 11 k7) / 120
        for (i = 0; i < dim; ++i)
          state[i] += h*(T(11)*(k[0][i] + k[6][i]) + T(81)*(k[2][i] + k[3][i]) - T(32)*(k[4][i] + k[5][i]))/T(120);
      }
    };


    /*
     Hairer order (p=8), Peter Stone's version

     Ref:
     * High order embedded Runge-Kutta formulae, by P.J.Prince and J.R.Dormand,
    Journal of Computational and Applied Mathematics, vol. 7, 1981, pages 67-75.
     * http://www.peterstone.name/Maplepgs/RKcoeff_8.html

    */
    template <class T, class F>
    class Hairer8 {

      int dim;

      T **a, *b, *c;

      public:

      F **k, *state_;

      const int order = 8;

      Hairer8 (int dim): dim(dim) {

        k = matrix<F> (14, dim);

        state_ = k[13];

        a = new T* [12];

        for (int i = 0; i < 12; ++i) a[i] = new T [i+1];

        b = new T [13];

        c = new T [12];

        Hairer8_coef(c, a, b);
      }

      ~Hairer8(){

        delete [] b;
        delete [] c;

        for (int i = 0; i < 12; ++i) delete [] a[i];
        delete [] a;

        free_matrix(k);
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F *k0 = 0){

        int i, j, l;

        F tmp;

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < 12; ++i) {

          for (l = 0; l < dim; ++l) {

            tmp = 0;
            for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

            state_[l] = state[l] + h*tmp;
          }

          deriv(t + c[i]*h, state_, k[i+1]);
        }

        for (l = 0; l < dim; ++l) {
          tmp = 0;
          for (i = 0; i < 13; ++i) tmp += b[i]*k[i][l];
          state[l] += h*tmp;
        }
      }
    };

    /*
      Hairer order (p=10), Peter Stone's version

      Ref:
      * "A Runge-Kutta Method of Order 10", by E. Hairer,
         Journal of the Institute of Mathematics and its Applications (1975) 16, pages 35 to 55.
      *  http://www.peterstone.name/Maplepgs/RKcoeff_10.html
    */
    template <class T, class F>
    class Hairer10 {

      int dim;

      T **a, *b, *c;

      public:

      F **k, *state_;

      const int order = 10;

      Hairer10 (int dim): dim(dim) {

        k = matrix<F>(18, dim);

        state_ = k[17];

        c = new T [16];

        a = new T* [16];

        for (int i = 0; i < 16; ++i) a[i] = new T [i+1];

        b = new T [17];

        Hairer10_coef(c, a, b);
      }

      ~Hairer10(){
        delete [] b;
        delete [] c;

        for (int i = 0; i < 16; ++i) delete [] a[i];
        delete [] a;

        free_matrix(k);
      }

      template <typename D>
      void step(const T &h, const T &t, F *state, D deriv, F* k0 = 0){

        int i, j, l;

        F tmp;

        if (k0 == 0)
          deriv(t, state, k[0]);
        else if (k0 != k[0])
          cp(k[0], k0, dim);

        for (i = 0; i < 16; ++i) {

          for (l = 0; l < dim; ++l) {

            tmp = 0;
            for (j = 0; j <= i; ++j) tmp += a[i][j]*k[j][l];

            state_[l] = state[l] + h*tmp;
          }

          deriv(t + c[i]*h, state_, k[i+1]);
        }

        for (l = 0; l < dim; ++l) {
          tmp = 0;
          for (i = 0; i < 17; i++) tmp += b[i]*k[i][l];
          state[l] += h*tmp;
        }
      }
    };

  } // namespace NonThreadSave::NonEmbedded

} // namespace NonThreadSave
