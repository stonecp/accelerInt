#ifndef __cklib_simd_h
#define __cklib_simd_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>
#include <string>
#include <typeinfo>

#ifdef __cplusplus
  #define RESTRICT __restrict__
#endif

// Normal object header.
#include "cklib.h"

namespace SIMD
{

// Master VCL header
#ifndef  MAX_VECTOR_SIZE
# define MAX_VECTOR_SIZE (256)
#endif
#define VCL_NAMESPACE VCL
#include "vcl/vectorclass.h"
#include "vcl/vectormath_exp.h"

using namespace VCL;

// Internal utility functions ...

template <typename ValueType>
inline ValueType sqr( const ValueType& x )
{
   return x*x;
}

template <typename ValueType>
inline ValueType __powu( const ValueType& x, const int p )
{
   if      (p == 0) return ValueType(1);
   else if (p == 1) return x;
   else if (p == 2) return x*x;
   else if (p == 3) return x*x*x;
   else return pow(x, double(p) );
}
template <typename ValueType>
inline ValueType __powi( const ValueType& x, const int p )
{
   if (p >= 0)
      return __powu(x, p);
   else
      return __powu((ValueType(1)/x), -p);
}

template <typename ValueType>
inline ValueType fmax( const ValueType& a, const ValueType& b )
{
   auto mask = a > b;
   return select( mask, a, b );
}
template <typename ValueType>
inline ValueType fmax( const ValueType& a, const double b )
{
   auto mask = a > ValueType(b);
   return select( mask, a, b );
}
template <typename ValueType>
inline ValueType fmin( const ValueType& a, const ValueType& b )
{
   auto mask = a < b;
   return select( mask, a, b );
}
template <typename ValueType>
inline ValueType fmin( const ValueType& a, const double b )
{
   auto mask = a < ValueType(b);
   return select( mask, a, b );
}

template <typename ValueType>
inline ValueType compute_H_RT (const int k,
                               const ValueType& T,
                               const ckdata_t *RESTRICT ck)
{
   auto f = [&] ( const double * RESTRICT a ) {
      return a[0] + a[5] / T + T * (a[1] / 2.0 + T * (a[2] / 3.0 + T * (a[3] / 4.0 + T * a[4] / 5.0)));
   };

   auto mask = ( T < ck->th_tmid[k] );

   return select( mask, f( ck->th_alo[k] ), f( ck->th_ahi[k] ) );
}

template <typename ValueType>
inline ValueType compute_Cp_R (const int k,
                               const ValueType& T,
                               const ckdata_t *RESTRICT ck)
{
   // Cp / R = Sum_(i=1)^(5){ a_i * T^(i-1) }

   auto f = [&] ( const double * RESTRICT a ) {
      return a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4])));
   };

   auto mask = ( T < ck->th_tmid[k] );

   return select( mask, f( ck->th_alo[k] ), f( ck->th_ahi[k] ) );
}

// Mixture Cp in mass units given mass fractions and temperature ... erg / (g * k)
template <typename ValueType>
inline ValueType ckcpbs (const ValueType& T,
                         const ValueType *RESTRICT y,
                         const ckdata_t  *RESTRICT ck)
{
   ValueType cp_mix ( 0.);

   #pragma ivdep
   for (int k = 0; k < ck->n_species; ++k)
   {
      ValueType cp_k = compute_Cp_R(k, T, ck);
      cp_k *= (__RU__ / ck->sp_mwt[k]);
      cp_k *= y[k];
      cp_mix += cp_k;
   }

   return cp_mix;
}
// Species enthalpies in mass units given temperature ... erg / g
template <typename ValueType>
inline void ckhms (const ValueType& T,
                         ValueType *RESTRICT h,
                   const ckdata_t  *RESTRICT ck)
{
   const ValueType RUT = __RU__ * T;

   #pragma ivdep
   for (int k = 0; k < ck->n_species; ++k)
   {
      ValueType h_k = compute_H_RT(k, T, ck);
      h_k *= (RUT / ck->sp_mwt[k]);

      h[k] = h_k;
   }
}

// Species S/R - H/RT ... special function.
template <typename ValueType>
inline void cksmh (const ValueType& T,
                   const ValueType& logT,
                         ValueType *RESTRICT smh,
                   const ckdata_t  *RESTRICT ck)
{
   //const double logTm1 = log(T) - 1.;
   const ValueType logTm1 = logT - 1.;
   const ValueType invT   = 1. / T;
   const ValueType T1     = T / 2.;
   const ValueType T2     = T*T / 6.;
   const ValueType T3     = T*T*T / 12.;
   const ValueType T4     = T*T*T*T / 20.;

   auto f = [&] ( const double * RESTRICT a ) {
      return a[0] * logTm1 + T1 * a[1] + T2 * a[2] + T3 * a[3] + T4 * a[4] - a[5] * invT + a[6];
   };

   #pragma ivdep
   for (int k = 0; k < ck->n_species; ++k)
   {
      auto mask = ( T < ck->th_tmid[k] );

      smh[k] = select( mask, f( ck->th_alo[k] ), f( ck->th_ahi[k] ) );
   }
}

// Mean molecular weight given mass fractions ... g / mol
template <typename ValueType>
inline ValueType ckmmwy (const ValueType *RESTRICT y,
                         const ckdata_t *RESTRICT ck)
{  
   // <W> = 1 / Sum_k { y_k / w_k }
   ValueType sumyow (0.0);
   for (int k = 0; k < ck->n_species; ++k)
      sumyow += (y[k] / ck->sp_mwt[k]);

   return ValueType(1) / sumyow;
}

// Mixture density given pressure, temperature and mass fractions ... g / cm^3
template <typename ValueType>
inline ValueType ckrhoy (const ValueType& p,
                         const ValueType& T,
                         const ValueType *RESTRICT y,
                         const ckdata_t  *RESTRICT ck)
{
   ValueType mean_mwt = ckmmwy(y, ck);

   // rho = p / (<R> * T) = p / (RU / <W> * T)

   //return p / (T * RU / mean_mwt);
   mean_mwt *= p;
   mean_mwt /= T;
   mean_mwt /= __RU__;
   return mean_mwt;
}

// Compute the molar concentration given p/T/y_k
template <typename ValueType>
inline void ckytcp (const ValueType& p,
                    const ValueType& T,
                    const ValueType *RESTRICT y,
                          ValueType *RESTRICT c,
                    const ckdata_t  *RESTRICT ck)
{
   // [c]_k = rho * y_k / w_k
   const ValueType rho = ckrhoy (p, T, y, ck);

   #pragma ivdep
   for (int k = 0; k < ck->n_species; ++k)
   {
      //c[k] = rho * y[k] / ck->sp_mwt[k];
      c[k] = rho * y[k] / ck->sp_mwt[k];
   }
}

// Compute temperature-dependent molar forward/reverse reaction rates
// ... utility function ...
template <typename ValueType>
inline void ckratt (const ValueType& T,
                          ValueType *RESTRICT smh,
                        //ValueType *RESTRICT eqk,
                          ValueType *RESTRICT rkf,
                          ValueType *RESTRICT rkr,
                     const ckdata_t *RESTRICT ck)
{
   const int kk = ck->n_species;
   const int ii = ck->n_reactions;

   const ValueType logT = log(T);
   const ValueType invT = 1.0 / T;
   const ValueType pfac = __PA__ / (__RU__ * T); // (dyne / cm^2) / (erg / mol / K) / (K)

   // I. Temperature-dependent rates ...

   // S/R - H/RT ... only needed for equilibrium.
   cksmh (T, logT, smh, ck);

   #pragma ivdep
   for (int i = 0; i < ii; ++i)
   {
      // Basic Arrhenius rates: A * exp( logT * b - E_R / T)
      rkf[i] = ck->rx_A[i] * exp( ck->rx_b[i] * logT - ck->rx_E[i] * invT );
   }

   #pragma ivdep
   for (int i = 0; i < ii; ++i)
   {
      // Irreversible reaction ...
      if (__is_enabled(ck->rx_info[i], __rx_flag_irrev))
      {
         rkr[i] = 0.0;
         //eqk[i] = 0.0;
      }
      // Reversible parameters ...
      else if (__is_enabled(ck->rx_info[i], __rx_flag_rparams))
      {
         rkr[i] = ck->rx_rev_A[i] * exp(ck->rx_rev_b[i] * logT - ck->rx_rev_E[i] * invT);
         //eqk[i] = rkf[i] / rkr[i];
      }
      // Use equilibrium for reversible rate ...
      else
      {
         // Sum_k { nu_k * (S/R - H/RT)_k }

         #define __nu (ck->rx_nu[i])
         #define __nuk (ck->rx_nuk[i])

         ValueType           sumsmh  = (__nu[0] * smh[__nuk[0]]);
         if (__nuk[1] != -1) sumsmh += (__nu[1] * smh[__nuk[1]]);
         if (__nuk[2] != -1) sumsmh += (__nu[2] * smh[__nuk[2]]);
                             sumsmh += (__nu[3] * smh[__nuk[3]]);
         if (__nuk[4] != -1) sumsmh += (__nu[4] * smh[__nuk[4]]);
         if (__nuk[5] != -1) sumsmh += (__nu[5] * smh[__nuk[5]]);

         #undef __nu
         #undef __nuk

         //eqk[__getIndex(i)] = exp(fmin(sumsmh, __exparg__));
         ValueType eqk_ = exp( fmin(sumsmh, __exparg__) );

         if (ck->rx_sumnu[i] != 0)
            eqk_ *= __powi(pfac, ck->rx_sumnu[i]);
            //eqk[__getIndex(i)] *= __powi(pfac,ck->rx_sumnu[i]);

         //if (!(ck->rx_info[i] & __rx_flag_irrev))
            //rkr[__getIndex(i)] = rkf[__getIndex(i)] / fmax(eqk[__getIndex(i)],__small__);
            rkr[i] = rkf[i] / fmax( eqk_, __small__);
      }
   }
}

template <typename ValueType>
inline void ckratc (const ValueType& T,
                    const ValueType *RESTRICT c,
                          ValueType *RESTRICT ctb,
                          ValueType *RESTRICT rkf,
                          ValueType *RESTRICT rkr,
                    const ckdata_t  *RESTRICT ck)
{
   const int kk = ck->n_species;
   const int ii = ck->n_reactions;

   const ValueType logT = log(T);
   const ValueType invT = 1.0 / T;

   // II. Concentration-dependent rates ...

   #pragma ivdep
   for (int i = 0; i < ii; ++i)
      ctb[i] = 1.0;

   // Third-body reactions ...
   if (ck->n_thdbdy > 0)
   {
      ValueType ctot( 0.0 );
      for (int k = 0; k < kk; ++k)
         ctot += c[k];

      for (int n = 0; n < ck->n_thdbdy; ++n)
      {
         const int rxn_idx = ck->rx_thdbdy_idx[n];

         ctb[rxn_idx] = ctot;

         // Add in the specific efficiencies ...

         for (int m = ck->rx_thdbdy_offset[n]; m < ck->rx_thdbdy_offset[n+1]; ++m)
         {
            const int k = ck->rx_thdbdy_spidx[m];
            ctb[rxn_idx] += (ck->rx_thdbdy_alpha[m] - 1.0) * c[k];
         }
      }
   }

   // Fall-off pressure dependencies ...
   if (ck->n_falloff > 0)
   {
      #pragma ivdep
      for (int n = 0; n < ck->n_falloff; ++n)
      {
         const int rxn_idx = ck->rx_falloff_idx[n];

         // Concentration of the third-body ... could be a specific species, too.
         ValueType cthb;
         if (ck->rx_falloff_spidx[n] != -1)
         {
            cthb = ctb[rxn_idx];
            ctb[rxn_idx] = 1.0;
         }
         else
            cthb = c[ ck->rx_falloff_spidx[n] ];

         #define __fpar (ck->rx_falloff_params[n])

         // Low-pressure limit rate ...
         ValueType rklow = __fpar[0] * exp(__fpar[1] * logT - __fpar[2] * invT);

         // Reduced pressure ...
         ValueType pr    = rklow * cthb / rkf[rxn_idx];

         // Correction ... k_infty (pr / (1+pr)) * F()
         ValueType p_cor;

         // Different F()'s ...
         //if (ck->rx_info[rxn_idx] & __rx_flag_falloff_sri)
         //{
         //   printf("SRI fall-off rxn not ready\n");
         //   exit(-1);
         //}
         //else if (ck->rx_info[rxn_idx] & __rx_flag_falloff_troe)
         if (__is_enabled(ck->rx_info[rxn_idx], __rx_flag_falloff_troe))
         {
            // 3-parameter Troe form ...
            ValueType Fcent = (1.0 - __fpar[3]) * exp(-T / __fpar[4]) + __fpar[3] * exp(-T / __fpar[5]);

            // Additional 4th (T**) parameter ...
            if (__is_enabled(ck->rx_info[rxn_idx], __rx_flag_falloff_troe4))
               Fcent += exp(-__fpar[6] * invT);

            ValueType log_Fc = log10( fmax(Fcent,__small__) );
            ValueType eta    = 0.75 - 1.27 * log_Fc;
            ValueType log_pr = log10( fmax(pr,__small__) );
            ValueType plus_c = log_pr - (0.4 + 0.67 * log_Fc);
          //ValueType _tmp   = plus_c / (eta - 0.14 * plus_c);
            ValueType log_F  = log_Fc / (1.0 + sqr(plus_c / (eta - 0.14 * plus_c)));
            ValueType Fc     = exp10(log_F);

            p_cor = Fc * (pr / (1.0 + pr));
         }
         else // Lindermann form
         {
            p_cor = pr / (1.0 + pr);
         }

         #undef __fpar

         rkf[rxn_idx] *= p_cor;
         rkr[rxn_idx] *= p_cor;
         //printf("%3d, %3d, %e, %e\n", n, rxn_idx, ck->rx_info[rxn_idx], p_cor, cthb);
      }

   } // fall-off's

   // II. Stoichiometry rates ...

   #pragma ivdep
   for (int i = 0; i < ii; ++i)
   {
      #define __nu (ck->rx_nu[i])
      #define __nuk (ck->rx_nuk[i])

      ValueType rkf_ = rkf[i] * ctb[i];
      ValueType rkr_ = rkr[i] * ctb[i];

                             rkf_ *= __powu( c[__nuk[0]],-__nu[0]);
      if (__nuk[1] != -1) {  rkf_ *= __powu( c[__nuk[1]],-__nu[1]);
         if (__nuk[2] != -1) rkf_ *= __powu( c[__nuk[2]],-__nu[2]);
      }

                             rkr_ *= __powu( c[__nuk[3]], __nu[3]);
      if (__nuk[4] != -1) {  rkr_ *= __powu( c[__nuk[4]], __nu[4]);
         if (__nuk[5] != -1) rkr_ *= __powu( c[__nuk[5]], __nu[5]);
      }

      #undef __nu
      #undef __nuk

      rkf[i] = rkf_;
      rkr[i] = rkr_;
   }
}

template <typename ValueType>
void ckwyp (const ValueType& p,
            const ValueType& T,
            const ValueType *RESTRICT y,
                  ValueType *RESTRICT wdot,
            const ckdata_t  *RESTRICT ck,
            ValueType rwk[])
{
   const int kk = ck->n_species;
   const int ii = ck->n_reactions;

   ValueType *RESTRICT rkf = rwk;
   ValueType *RESTRICT rkr = rkf + (ii);
   ValueType *RESTRICT ctb = rkr + (ii);
   ValueType *RESTRICT c   = ctb + (ii);
   ValueType *RESTRICT smh = c;
 //ValueType *RESTRICT eqk = ctb;

   // Compute temperature-dependent forward/reverse rates ... mol / cm^3 / s
   ckratt (T, smh, /*NULL,*/ rkf, rkr, ck);

   // Convert to molar concentrations ... mol / cm^3
   ckytcp (p, T, y, c, ck);

   // Compute concentration-dependent forward/reverse rates ... mol / cm^3 / s
   ckratc (T, c, ctb, rkf, rkr, ck);

   // Compute species net production rates ... mol / cm^3 / s

   for (int k = 0; k < kk; ++k)
      wdot[k] = 0.0;

   for (int i = 0; i < ii; ++i)
   {
      const ValueType rop = rkf[i] - rkr[i];

      #define __nu (ck->rx_nu[i])
      #define __nuk (ck->rx_nuk[i])

                             wdot[__nuk[0]] += (rop * __nu[0]);
      if (__nuk[1] != -1) {  wdot[__nuk[1]] += (rop * __nu[1]);
         if (__nuk[2] != -1) wdot[__nuk[2]] += (rop * __nu[2]);
      }

                             wdot[__nuk[3]] += (rop * __nu[3]);
      if (__nuk[4] != -1) {  wdot[__nuk[4]] += (rop * __nu[4]);
         if (__nuk[5] != -1) wdot[__nuk[5]] += (rop * __nu[5]);
      }

      #undef __nu
      #undef __nuk
   }
}

// Meta-function to compute RHS for constant-pressure reaction system
template <typename ValueType>
void ckrhs (const ValueType& p,
            const ValueType& T,
                  ValueType& Tdot,
            const ValueType  yk[],
                  ValueType  ykdot[],
            const ckdata_t *RESTRICT ck,
                  ValueType  rwk[])
{
   const int kk = ck->n_species;

   /* Compute local density given p/T/y_k */
   //const double rho = ckrhoy(p, T, y, ck);

   /* Compute the molar concentration ( mol / cm^3 ) */
   //ckytcr (rho, y, ydot /* as c[] */, ck);

   /* Compute molar reaction rate. (mol / (s*cm^3) */
   //ckwc (T, ydot /* as c[]*/, ydot, ck);
   ckwyp (p, T, yk, ykdot, ck, rwk);

   /* Compute mixture Cp (ergs / gm*K) */
   const ValueType cp_mix = ckcpbs(T, yk, ck);

   /* Compute species enthalpy (ergs / K) */
   ckhms(T, rwk, ck);

   /* Extract the molecular weights of the species ... this could just be a pointer. */
   const ValueType rho = ckrhoy(p, T, yk, ck);

   Tdot = 0.0;
   for (int k = 0; k < kk; ++k)
   {
      /* Convert from molar to mass units. */
      ykdot[k] *= ck->sp_mwt[k];
      ykdot[k] /= rho;

      /* Sum up the net enthalpy change. */
      Tdot -= (rwk[k] * ykdot[k]);
   }

   Tdot /= cp_mix;

   return;
}

template <typename Func>
void test_simd_rhs ( const int numProblems, const double *u_in, Func& func, const ckdata_t *RESTRICT ck )
{
   const int kk = ck->n_species;
   const int neq = kk+1;

   typedef VCL::Vec4d SimdType;
   const int VectorLength = sizeof(SimdType) / sizeof(double);
   
   printf("Instruction Set= %d %s %d\n", INSTRSET, typeid(SimdType).name(), VectorLength);

   const double p = func.getPressure();

   VectorType<double,64> scalar_out( neq * numProblems );
   VectorType<double,64> vector_out( neq * numProblems );

   SimdType v_p(p);

   alignas(64) double T0[VectorLength];
   alignas(64) double Tdot[VectorLength];
   VectorType<double,64> y0(VectorLength*kk);
   VectorType<double,64> ykdot(VectorLength*kk);
   VectorType<SimdType,64> v_rwk( ck_lenrwk(ck) );

   double sum_v = 0, sum_s = 0;

   for (int iter = 0; iter < 1; ++iter)
   {

   for (int i = 0; i < numProblems; ++i)
      func.rhs( neq, 0.0, const_cast<double*>(&u_in[neq*i]), &scalar_out[neq*i] );

   for (int i0 = 0; i0 < numProblems; i0 += VectorLength)
   {
      if ( i0 + VectorLength < numProblems )
      {
         for (int i = 0; i < VectorLength; ++i)
         {
            const double *ui = u_in + neq*(i0+i);
            T0[i] = ui[ getTempIndex(neq) ];
            for (int k = 0; k < kk; ++k)
               y0[k*VectorLength+i] = ui[ getFirstSpeciesIndex(neq)+k ];
         }

         SimdType *v_T = (SimdType *) &T0;
         SimdType *v_Tdot = (SimdType *) &Tdot;
         SimdType *v_y0 = (SimdType *) y0.getPointer();
         SimdType *v_ykdot = (SimdType *) ykdot.getPointer();

         SIMD::ckrhs( v_p, *v_T, *v_Tdot, v_y0, v_ykdot, ck, v_rwk.getPointer() );

         for (int i = 0; i < VectorLength; ++i)
         {
            double *v_out = vector_out.getPointer() + neq*(i0+i);
            v_out[ getTempIndex(neq) ] = Tdot[i];
            for (int k = 0; k < kk; ++k)
               v_out[ getFirstSpeciesIndex(neq)+k ] = v_ykdot[k].extract(i);
         }
      }
      else
      {
         for (int i = i0; i < numProblems; ++i)
            func.rhs( neq, 0.0, const_cast<double*>(&u_in[neq*i]), &vector_out[neq*i] );
      }
   }

   sum_v += vector_out[iter];
   sum_s += scalar_out[iter];

   }

   {
      double err2, ref2 = 0;
      for (int i = 0; i < numProblems; ++i)
      {
         const double *v_out = vector_out.getPointer() + neq*i;
         const double *s_out = scalar_out.getPointer() + neq*i;
         double diff = s_out[ getTempIndex(neq) ]
                     - v_out[ getTempIndex(neq) ];
         err2 += sqr( diff );
         ref2 += sqr( s_out[ getTempIndex(neq) ] );
      }

      printf("err2= %e %e %e %d\n", err2, ref2, std::sqrt(err2)/std::sqrt(ref2), numProblems % VectorLength);
   }
   for (int k = 0; k < kk; ++k)
   {
      double err2, ref2 = 0;
      for (int i = 0; i < numProblems; ++i)
      {
         const double *v_out = vector_out.getPointer() + neq*i;
         const double *s_out = scalar_out.getPointer() + neq*i;
         double diff = s_out[ getFirstSpeciesIndex(neq)+k ]
                     - v_out[ getFirstSpeciesIndex(neq)+k ];
         err2 += sqr( diff );
         ref2 += sqr( s_out[ getFirstSpeciesIndex(neq)+k ] );
      }

      printf("err2= %e %e %e\n", err2, ref2, std::sqrt(err2)/std::sqrt(ref2));
   }

   return;
}

} // namespace

#endif // ifndef
