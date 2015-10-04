#include "header.cuh"
#include <math.h>
#include <float.h>
#include "dydt.cuh"
#include "gpu_macros.cuh"

#define FD_ORD 1

// Finite difference coefficients
#if FD_ORD == 2
  __constant__ double x_coeffs[FD_ORD] = {-1.0, 1.0};
  __constant__ double y_coeffs[FD_ORD] = {-0.5, 0.5};
#elif FD_ORD == 4
  __constant__ double x_coeffs[FD_ORD] = {-2.0, -1.0, 1.0, 2.0};
  __constant__ double y_coeffs[FD_ORD] = {1.0 / 12.0, -2.0 / 3.0, 2.0 / 3.0, -1.0 / 12.0};
#elif FD_ORD == 6
  __constant__ double x_coeffs[FD_ORD] = {-3.0, -2.0, - 1.0, 1.0, 2.0, 3.0};
  __constant__ double y_coeffs[FD_ORD] = {-1.0 / 60.0, 3.0 / 20.0, -3.0 / 4.0, 3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0};
#endif

__device__
void eval_jacob (const double t, const double pres, double * y, double * jac) {
  double dy[NN];
  double error[NN];

  dydt (t, pres, y, dy);
  
  #pragma unroll
  for (int i = 0; i < NN; ++i) {
    error[i] = ATOL + (RTOL * fabs(y[i]));
  }
  
  // unit roundoff of machine
  double srur = sqrt(DBL_EPSILON);
  
  double sum = 0.0;
  #pragma unroll
  for (int i = 0; i < NN; ++i) {
    sum += (error[i] * dy[i]) * (error[i] * dy[i]);
  }
  double fac = sqrt(sum / ((double)(NN)));
  double r0 = 1000.0 * RTOL * DBL_EPSILON * ((double)(NN)) * fac;
  
#ifndef GLOBAL_MEM
  double f_temp[NN];
#endif
  
  #pragma unroll
  for (int j = 0; j < NN; ++j) {
    double yj_orig = y[j];
    double r = fmax(srur * fabs(yj_orig), r0 / error[j]);
    
    #if FD_ORD == 1
      y[j] = yj_orig + r;
      dydt (t, pres, y, f_temp);
        
      #pragma unroll
      for (int i = 0; i < NN; ++i) {
        jac[i + NN*j] = (f_temp[i] - dy[i]) / r;
      }
    #else
      #pragma unroll
      for (int i = 0; i < NN; ++i) {
        jac[i + NN*j] = 0.0;
      }
      #pragma unroll
      for (int k = 0; k < FD_ORD; ++k) {
        y[j] = yj_orig + x_coeffs[k] * r;
        dydt (t, pres, y, f_temp);
        
        #pragma unroll
        for (int i = 0; i < NN; ++i) {
          jac[i + NN*j] += y_coeffs[k] * f_temp[i];
        }
      }
      #pragma unroll
      for (int i = 0; i < NN; ++i) {
        jac[i + NN*j] /= r;
      }
    #endif
    
    y[j] = yj_orig;
  }
  
}