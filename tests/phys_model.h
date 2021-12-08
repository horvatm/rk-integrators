#if !defined(__phys_model_h)
#define __phys_model_h

#if defined(EMBEDDED)
  //#define METHOD RKM4
  #define METHOD Hairer8
#else
  //#define METHOD RK4
  #define METHOD Hairer8
  //#define METHOD Hairer10
#endif

typedef double real;

// Energy of the pendulum, x[0] = q, x[1] = p
real energy(real *x){
  return 1 - cos(x[0]) + 0.5*x[1]*x[1];
}

// Time derivative of the trajectory of the pendulum
// x[0] = q, x[1] = p
namespace perturbations {
  class Perturbation {
    
    public:
    
    int get_dim() const {
      return 2;
    }
      
    void getQPdot(const real &lambda, const real* const x, real* dxdt){
      dxdt[0] = x[1];
      dxdt[1] = -std::sin(x[0]);
    }
  };
}
#endif
