#ifndef __cklib_simd_h
#define __cklib_simd_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <cmath>
#include <string>
#include <typeinfo>

#ifdef __cplusplus
  #define RESTRICT __restrict__
#endif

// Local header(s).
#include "cklib.h"
#include "clock.h"

namespace SIMD
{

const int Alignment = 64;

// Master VCL header
#ifndef  MAX_VECTOR_SIZE
# define MAX_VECTOR_SIZE (512)
#endif
#define VCL_NAMESPACE VCL
#include "vcl/vectorclass.h"
#include "vcl/vectormath_exp.h"

using namespace VCL;

template <typename T, int VL>
struct VCL_TypeSelector;

template<> struct VCL_TypeSelector<double,2> { typedef VCL::Vec2d value_type; };
template<> struct VCL_TypeSelector<double,4> { typedef VCL::Vec4d value_type; };
template<> struct VCL_TypeSelector<double,8> { typedef VCL::Vec8d value_type; };
template<> struct VCL_TypeSelector<int64_t,2> { typedef VCL::Vec2q value_type; };
template<> struct VCL_TypeSelector<int64_t,4> { typedef VCL::Vec4q value_type; };
template<> struct VCL_TypeSelector<int64_t,8> { typedef VCL::Vec8q value_type; };

template <typename V>
struct VCL_MaskSelector;

template<> struct VCL_MaskSelector<VCL::Vec2d> { typedef VCL::Vec2db mask_type; };
template<> struct VCL_MaskSelector<VCL::Vec4d> { typedef VCL::Vec4db mask_type; };
template<> struct VCL_MaskSelector<VCL::Vec8d> { typedef VCL::Vec8db mask_type; };
template<> struct VCL_MaskSelector<VCL::Vec2q> { typedef VCL::Vec2qb mask_type; };
template<> struct VCL_MaskSelector<VCL::Vec4q> { typedef VCL::Vec4qb mask_type; };
template<> struct VCL_MaskSelector<VCL::Vec8q> { typedef VCL::Vec8qb mask_type; };

template <typename V>
struct VCL_Length;

template<> struct VCL_Length<VCL::Vec2d> { enum { length = 2 }; };
template<> struct VCL_Length<VCL::Vec4d> { enum { length = 4 }; };
template<> struct VCL_Length<VCL::Vec8d> { enum { length = 8 }; };

template <typename MaskType>
inline bool all( const MaskType& mask ) { return horizontal_and( mask ); }

template <typename MaskType>
inline bool any( const MaskType& mask ) { return horizontal_or( mask ); }

template <typename SimdType>
std::string toString (const SimdType& x)
{
   const int len = SimdType::size();
   std::ostringstream oss;
   oss << "[";
   for (int i = 0; i < len; ++i) {
      oss << x[i];
      if (i != len-1) oss << ",";
   }
   oss << "]";
   
   return std::string( oss.str() );
}

std::ostream& operator<< ( std::ostream& os, const VCL::Vec2d& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec4d& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec8d& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec4i& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec8i& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec2q& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec4q& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec8q& obj) { return ( os << toString(obj) ); }

std::ostream& operator<< ( std::ostream& os, const VCL::Vec2db& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec4db& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec8db& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec4ib& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec8ib& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec2qb& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec4qb& obj) { return ( os << toString(obj) ); }
std::ostream& operator<< ( std::ostream& os, const VCL::Vec8qb& obj) { return ( os << toString(obj) ); }

template <typename ValueType, int VectorLength>
inline
std::ostream& operator<< ( std::ostream& os,
                           const typename VCL_TypeSelector<ValueType,VectorLength>::value_type& obj )
{
   // stream obj's data into os
   os << toString(obj);
   /*const int len = SimdType::size();
   std::ostringstream oss;
   os << "[";
   for (int i = 0; i < len; ++i) {
      os << x[i];
      if (i != len-1) oss << ",";
   }
   os << "]";*/
   return os;
}

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

template <typename SimdType>
struct simd_cklib_functor
{
   const ckdata_t *m_ckptr;
   const SimdType m_pres;
   VectorType<SimdType,Alignment> m_rwk;

   simd_cklib_functor( const ckdata_t *ckptr, const double p = __PA__ )
      : m_ckptr(ckptr), m_pres(p)
      {
         this->m_rwk.resize( VCL_Length<SimdType>::length * ck_lenrwk(m_ckptr) );
      }

   SimdType getPressure(void) const { return m_pres; }

   int operator() (const int &neq, const SimdType& time, SimdType y[], SimdType f[])
   {
      return this->rhs( neq, time, y, f );
   }

   int rhs (const int &neq, const SimdType& time, SimdType y[], SimdType f[])
   {
      const int kk = this->m_ckptr->n_species;

      const SimdType T = y[ getTempIndex(neq) ];
      const SimdType *yk = y + getFirstSpeciesIndex(neq);

      SimdType *ykdot = f + getFirstSpeciesIndex(neq);
      SimdType *Tdot = f + getTempIndex(neq);

      const SimdType &p = this->m_pres;

      ckrhs ( p, T, *Tdot, yk, ykdot, this->m_ckptr, this->m_rwk.getPointer() );

      return 0;
   }

   int jac (const int &neq, double t, double y[], double fy[]);
};


template <typename Func>
void test_simd_rhs ( const int numProblems, const double *u_in, Func& func, const ckdata_t *RESTRICT ck )
{
   const int kk = ck->n_species;
   const int neq = kk+1;

   typedef typename VCL_TypeSelector<double,4>::value_type SimdType;
   typedef typename VCL_MaskSelector<SimdType>::mask_type MaskType;
   //const int VectorLength = sizeof(SimdType) / sizeof(double);
   const int VectorLength = VCL_Length<SimdType>::length;

   printf("Instruction Set= %d %s %d %s\n", INSTRSET, typeid(SimdType).name(), VectorLength, typeid( MaskType).name());

   const double p = func.getPressure();

   VectorType<double,Alignment> scalar_out( neq * numProblems );
   VectorType<double,Alignment> vector_out( neq * numProblems );

   alignas(Alignment) double T0[VectorLength];
   alignas(Alignment) double Tdot[VectorLength];
   VectorType<double,Alignment> y0(VectorLength*kk);
   VectorType<double,Alignment> ykdot(VectorLength*kk);
   VectorType<SimdType,Alignment> v_rwk( ck_lenrwk(ck) );

   VectorType<double,Alignment> _u(VectorLength*neq);
   VectorType<double,Alignment> _f(VectorLength*neq);

   simd_cklib_functor<SimdType> simd_func( ck );

   double sum_v = 0, sum_s = 0;

   double time_scalar = 0, time_vector = 0;

   for (int iter = 0; iter < 20; ++iter)
      if (iter % 2 == 0)
      {
         double time_start = WallClock();

         for (int i = 0; i < numProblems; ++i)
            func.rhs( neq, 0.0, const_cast<double*>(&u_in[neq*i]), &scalar_out[neq*i] );

         sum_s += scalar_out[iter/2];
         time_scalar += (WallClock() - time_start);
      }
      else
      {
         double time_start = WallClock();
         for (int i0 = 0; i0 < numProblems; i0 += VectorLength)
         {
            if ( i0 + VectorLength < numProblems )
            {
               for (int i = 0; i < VectorLength; ++i)
               {
                  const double *ui = u_in + neq*(i0+i);
                  for (int j = 0; j < neq; ++j)
                     _u[j*VectorLength+i] = ui[j];
               }

               SimdType *v_u = (SimdType *) _u.getPointer();
               SimdType *v_f = (SimdType *) _f.getPointer();

               simd_func.rhs ( neq, SimdType(0), v_u, v_f );

               for (int i = 0; i < VectorLength; ++i)
               {
                  double *v_out = vector_out.getPointer() + neq*(i0+i);
                  for (int j = 0; j < neq; ++j)
                     v_out[j] = _f[j*VectorLength+i];
               }
            }
            else
            {
               for (int i = i0; i < numProblems; ++i)
                  func.rhs( neq, 0.0, const_cast<double*>(&u_in[neq*i]), &vector_out[neq*i] );
            }
         }

         sum_v += vector_out[iter/2];
         time_vector += (WallClock() - time_start);
      }

   printf("SIMD timer: %f %f %.1f %d %e\n", 1000.*time_vector, 1000.*time_scalar, time_scalar/time_vector, sum_s == sum_v, fabs(sum_s-sum_v)/fabs(sum_s));

   {
      double err2 = 0, ref2 = 0;
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
   {
      double err2 = 0, ref2 = 0;
      for (int k = 0; k < kk; ++k)
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

enum SolverTags { RKF45 = 1, ROS4 };

enum ErrorFlags { ERR_SUCCESS         = 0,
                  ERR_TOO_MUCH_WORK   = 1,
                  ERR_TDIST_TOO_SMALL = 2,
                  ERR_HIN_MAX_ITERS   = 3 };

const ErrorFlags GetErrorInt ( const int errflag )
{
   switch ( errflag )
   {
      case ( int(ERR_SUCCESS) ) :
         return ERR_SUCCESS;
      case ( ERR_TOO_MUCH_WORK ) :
         return ERR_TOO_MUCH_WORK;
      case ( ERR_TDIST_TOO_SMALL ) :
         return ERR_TDIST_TOO_SMALL;
      case ( ERR_HIN_MAX_ITERS ) :
         return ERR_HIN_MAX_ITERS;
      default:
      {
         fprintf(stderr,"Unknown error value: %d\n", errflag);
         exit(2);
      }
   }
}
const char* GetErrorString ( const ErrorFlags& errflag )
{
   switch ( errflag )
   {
      case ( ERR_SUCCESS ) :
      {
         static const char* flag = "Success (None)";
         return flag;
      }
      case ( ERR_TOO_MUCH_WORK ) :
      {
         static const char* flag = "Too much work in solver";
         return flag;
      }
      case ( ERR_TDIST_TOO_SMALL ) :
      {
         static const char* flag = "Time too small; not resolvable";
         return flag;
      }
      case ( ERR_HIN_MAX_ITERS ) :
      {
         static const char* flag = "Exceeded maximum hin() iterations";
         return flag;
      }
   }
}
const char* GetErrorString ( const int errflag )
{
   return GetErrorString( GetErrorInt(errflag) );
}

template <typename T>
struct CommonSolverType
{
   typedef T ValueType;
   typedef typename VCL_MaskSelector< ValueType >::mask_type MaskType;

   enum { vector_length = VCL_Length<ValueType>::length,
          vlen = vector_length };

   struct CountersType
   {
      typedef typename VCL_TypeSelector<int64_t,vlen>::value_type IndexType;
      int nit;
      IndexType nst;
      IndexType nfe;
      IndexType nje;
      IndexType nlu;
      IndexType errflag;
   };

   int neq;

   int itol;
   double s_rtol, s_atol;
 //double v_rtol[], v_atol[];

   int max_iters, min_iters;

   double adaption_limit;

   double h_min, h_max;
   double t_round;
   double t_stop;

   inline double uround(void) const { return DBL_EPSILON; }

   CommonSolverType (const int _neq)
      : neq(_neq),
        min_iters(1), max_iters(100),
        h_min(0), h_max(0),
        adaption_limit(4),
        itol(1),
        s_rtol(1.0e-11), s_atol(1.0e-9),
        t_stop(-1)
   {}

   virtual ~CommonSolverType()
   {}

private:
   CommonSolverType (void);

public:
   int init (const double t0, const double t_stop)
   {
      this->t_stop = t_stop;

      this->h_min = 0.0;
      this->h_max = 0.0;

      const double t_dist = this->t_stop - t0;
      this->t_round = t_dist * this->uround();

      if (t_dist < (this->t_round * 2.0))
      {
         //fprintf(stderr, "error: tdist < 2.0 * tround %e\n", tdist);
         return RK_TDIST_TOO_SMALL;
      }

      if (this->h_min < this->t_round) this->h_min = this->t_round * 100.0;
      if (this->h_max < this->t_round) this->h_max = t_dist / (double)this->min_iters;

      return ERR_SUCCESS;
   }

   ValueType wnorm (const ValueType *RESTRICT x, const ValueType *RESTRICT y)
   {
      const int neq = this->neq;
      ValueType sum = 0;
      for (int k = 0; k < neq; k++)
      {
         ValueType ewt = (this->s_rtol * abs(y[k])) + this->s_atol;
         ValueType prod = x[k] / ewt;
         sum += (prod*prod);
      }

      const double denom = 1.0 / neq;
      return sqrt(sum * denom);
   }

   template <class Functor>
   int hin ( const ValueType t, ValueType *h0, ValueType *RESTRICT y, ValueType *RESTRICT rwk, Functor& func)
   {
      //value_type tround = tdist * this->uround();
      //double tdist = t_stop - t;
      //double tround = tdist * rk_uround();

      // Set lower and upper bounds on h0, and take geometric mean as first trial value.
      // Exit with this value if the bounds cross each other.

      //rk->h_min = fmax(tround * 100.0, rk->h_min);
      //rk->h_max = fmin(tdist, rk->h_max);

      //const int neq = this->neq;

      ValueType *RESTRICT ydot  = rwk;
      ValueType *RESTRICT y1    = ydot + neq;
      ValueType *RESTRICT ydot1 = y1 + neq;

      int need_ydot = 1;

      // Adjust upper bound based on ydot ...
   /*    if (0)
         {
            need_ydot = false;

            // compute ydot at t=t0
            func (neq, y, ydot);
            ++this->nfe;

            for (int k = 0; k < neq; k++)
            {
               value_type dely = 0.1 * fabs(y[k]) + this->atol;
               value_type hub0 = hub;
               if (hub * fabs(ydot[k]) > dely) hub = dely / fabs(ydot[k]);
               //printf("k=%d, hub0 = %e, hub = %e\n", k, hub0, hub);
            }
         }*/

      double hlb = this->h_min;
      double hub = this->h_max;

      ValueType hg = std::sqrt(hlb*hub);

      if (hub < hlb)
      {
         *h0 = hg;
         return ERR_SUCCESS;
      }

      // Start iteration to find solution to ... {WRMS norm of (h0^2 y'' / 2)} = 1

      const int miters = 10;
      MaskType hnew_is_ok( false );
      ValueType hnew = hg;
      int iter = 0;
      int ierr = ERR_SUCCESS;

      // compute ydot at t=t0
      if (need_ydot)
      {
         func(neq, 0.0, y, ydot);
         //++rk->nfe;
         need_ydot = 0;
      }

      while(1)
      {
         // Estimate y'' with finite-difference ...
         //double t1 = hg;

         #pragma ivdep
         for (int k = 0; k < neq; k++)
            y1[k] = y[k] + hg * ydot[k];

         // compute y' at t1
         func (neq, 0.0, y1, ydot1);
         //++rk->nfe;

         // Compute WRMS norm of y''
         #pragma ivdep
         for (int k = 0; k < neq; k++)
            y1[k] = (ydot1[k] - ydot[k]) / hg;

         ValueType yddnrm = this->wnorm ( y1, y );

         //std::cout << "iter " << iter << " hg " << hg << " y'' " << yddnrm << std::endl;
         //std::cout << "ydot " << ydot[neq-1] << std::endl;

         // should we accept this?
         hnew = select( hnew_is_ok | MaskType( iter == miters ), hg, hnew );
         if ( all(hnew_is_ok) )
         {
            ierr = ERR_SUCCESS;
            break;
         }
         else if (iter == miters)
         {
            ierr = ERR_HIN_MAX_ITERS;
            break;
         }

         // Get the new value of h ...
         {
            auto mask = (yddnrm*hub*hub > ValueType(2.0));
            hnew = select( mask, sqrt(2.0 / yddnrm), sqrt(hg * hub) );
         }

         // test the stopping conditions.
         ValueType hrat = hnew / hg;

         // Accept this value ... the bias factor should bring it within range.
         hnew_is_ok = (hrat > 0.5) && (hrat < 2.0);

         // If y'' is still bad after a few iterations, just accept h and give up.
         if ( iter > 1 ) {
            //hnew_is_ok = (hrat > 2.0);
            //hnew = select( hnew_is_ok, hg, hnew );
            auto mask = (hrat > 2.0);
            hnew = select( mask, hg, hnew );
            //hnew_is_ok = select( mask, mask, hnew_is_ok );
            hnew_is_ok |= mask;
         }

         //printf("iter=%d, yddnrw=%e, hnew=%e, hlb=%e, hub=%e\n", iter, yddnrm, hnew, hlb, hub);

         hg = hnew;
         iter ++;
      }

      // bound and bias estimate
      *h0 = hnew * 0.5;
      *h0 = max(*h0, hlb);
      *h0 = min(*h0, hub);

      //printf("h0=%e, hlb=%e, hub=%e\n", h0, hlb, hub);

      return ierr;
   }
};

template <typename ValueType>//, SolverTags SolverTag >
struct SimdERKSolverType : public CommonSolverType<ValueType>
{
   typedef CommonSolverType<ValueType> Parent;

   enum { vlen = Parent::vlen };

   typedef typename Parent::MaskType MaskType;
   typedef typename Parent::CountersType CountersType;

   VectorType<ValueType> m_rwk;

   SimdERKSolverType(const int _neq)
      : Parent(_neq)
   {
      m_rwk.resize( get_lenrwk() );
   }

   constexpr size_t get_lenrwk(const int _neq) const { return (8 * _neq); }
   constexpr size_t get_lenrwk(void) const { return (8 * this->neq); }

   template <class Functor>
   int oneStep (const ValueType& h,
                      ValueType *RESTRICT y,
                      ValueType *RESTRICT y_out,
                      ValueType *RESTRICT rwk,
                      Functor&   func)
   {
      static
      const double c20 = 0.25,
                   c21 = 0.25,
                   c30 = 0.375,
                   c31 = 0.09375,
                   c32 = 0.28125,
                   c40 = 0.92307692307692,
                   c41 = 0.87938097405553,
                   c42 =-3.2771961766045,
                   c43 = 3.3208921256258,
                   c51 = 2.0324074074074,
                   c52 =-8.0,
                   c53 = 7.1734892787524,
                   c54 =-0.20589668615984,
                   c60 = 0.5,
                   c61 =-0.2962962962963,
                   c62 = 2.0,
                   c63 =-1.3816764132554,
                   c64 = 0.45297270955166,
                   c65 =-0.275,
                   a1 = 0.11574074074074,
                   a2 = 0.0,
                   a3 = 0.54892787524366,
                   a4 = 0.5353313840156,
                   a5 =-0.2,
                   b1 = 0.11851851851852,
                   b2 = 0.0,
                   b3 = 0.51898635477583,
                   b4 = 0.50613149034201,
                   b5 =-0.18,
                   b6 = 0.036363636363636;

      const int neq = this->neq;

      // local dependent variables (5 total)
      ValueType *RESTRICT f1   = rwk ;
      ValueType *RESTRICT f2   = rwk + (  neq) ;
      ValueType *RESTRICT f3   = rwk + (2*neq) ;
      ValueType *RESTRICT f4   = rwk + (3*neq) ;
      ValueType *RESTRICT f5   = rwk + (4*neq) ;
      ValueType *RESTRICT f6   = rwk + (5*neq) ;
      ValueType *RESTRICT ytmp = rwk + (6*neq) ;

      // 1)
      func(neq, 0.0, y, f1);

      for (int k = 0; k < neq; k++)
      {
         //f1[k] = h * ydot[k];
         f1[k] *= h;
         ytmp[k] = y[k] + c21 * f1[k];
      }

      // 2)
      func(neq, 0.0, ytmp, f2);

      for (int k = 0; k < neq; k++)
      {
         //f2[k] = h * ydot[k];
         f2[k] *= h;
         ytmp[k] = y[k] + c31 * f1[k] + c32 * f2[k];
      }

      // 3)
      func(neq, 0.0, ytmp, f3);

      for (int k = 0; k < neq; k++) {
         //f3[k] = h * ydot[k];
         f3[k] *= h;
         ytmp[k] = y[k] + c41 * f1[k] + c42 * f2[k] + c43 * f3[k];
      }

      // 4)
      func(neq, 0.0, ytmp, f4);

      for (int k = 0; k < neq; k++) {
         //f4[k] = h * ydot[k];
         f4[k] *= h;
         ytmp[k] = y[k] + c51 * f1[k] + c52 * f2[k] + c53 * f3[k] + c54 * f4[k];
      }

      // 5)
      func(neq, 0.0, ytmp, f5);

      for (int k = 0; k < neq; k++) {
         //f5[k] = h * ydot[k];
         f5[k] *= h;
         ytmp[k] = y[k] + c61*f1[k] + c62*f2[k] + c63*f3[k] + c64*f4[k] + c65*f5[k];
      }

      // 6)
      func(neq, 0.0, ytmp, f6);

      for (int k = 0; k < neq; k++)
      {
         //const T f6 = h * ydot[k];
         f6[k] *= h;

         // 5th-order RK value.
         const ValueType r5 = b1*f1[k] + b3*f3[k] + b4*f4[k] + b5*f5[k] + b6*f6[k];

         // 4th-order RK residual.
         const ValueType r4 = a1*f1[k] + a3*f3[k] + a4*f4[k] + a5*f5[k];

         // Trucation error: difference between 4th and 5th-order RK values.
         rwk[k] = abs(r5 - r4);

         // Update solution.
         y_out[k] = y[k] + r5; // Local extrapolation
      }

      return ERR_SUCCESS;
   }

   template <typename Functor, typename Counters>
   int solve (ValueType* tcur, ValueType *hcur, Counters* counters, ValueType y[], Functor& func)
   {
      int ierr = ERR_SUCCESS;

      ValueType *rwk = m_rwk.getPointer();

      // Estimate the initial step size ...
      {
         auto mask = (*hcur < this->h_min);
         if ( any(mask) )
         {
            ValueType h0 = *hcur;
            //std::cout << "Before hin: " << h0 << ", " << mask << "\n";
            ierr = this->hin ( *tcur, &h0, y, rwk, func );
            if (ierr != ERR_SUCCESS)
            {
               fprintf(stderr,"Failure solve(): %d %s\n", ierr, GetErrorString(ierr));
               return ierr;
            }
            //std::cout << "After hin: " << h0 << "\n";

            *hcur = select( mask, h0, *hcur );
         }
         //printf("hin = %s %s %e %e\n", toString(*hcur).c_str(), toString(mask).c_str(), this->h_min, this->h_max);
      }

      #define t (*tcur)
      #define h (*hcur)
      #define iter (counters->nit)
      #define nst (counters->nst)

      nst = 0;
      iter = 0;

      MaskType not_done = abs(t - this->t_stop) > ValueType( this->t_round );

      while (  any(not_done) )
      {
         ValueType *ytry = rwk + this->neq*7;

         // Take a trial step over h_cur ...
         this->oneStep ( h, y, ytry, rwk, func );

         ValueType herr = max(1.0e-20, this->wnorm ( rwk, y ));

         // Is there error acceptable?
         MaskType accept = ((herr <= 1.0) | (h <= this->h_min)) & not_done;

         // update solution ...
         if ( any(accept) )
         {
            t   = select ( accept, t + h, t   );
            nst = select ( accept, nst+1, nst );

            for (int k = 0; k < this->neq; k++)
               y[k] = select( accept, ytry[k], y[k] );

            not_done = abs(t - this->t_stop) > this->t_round;
         }

         ValueType fact = sqrt( sqrt(1.0 / herr) ) * (0.840896415);

         // Restrict the rate of change in dt
         fact = max(fact, 1.0 / this->adaption_limit);
         fact = min(fact,       this->adaption_limit);

#if defined(VERBOSE) && (VERBOSE > 0)
         if (iter % 10 == 0)
            printf("iter = %d: accept=%s, not_done=%s t=%s, fact=%s %s %s\n", iter, toString(accept).c_str(), toString(not_done).c_str(), toString(t).c_str(), toString(fact).c_str(), toString(y[this->neq-1]).c_str(), toString(h).c_str());
#endif

         // Apply grow/shrink factor for next step.
         h = select( not_done, h * fact, h);

         // Limit based on the upper/lower bounds
         h = fmin(h, this->h_max);
         h = fmax(h, this->h_min);

         // Stretch the final step if we're really close and we didn't just fail ...
         h = select( accept & ( abs((t + h) - this->t_stop) < this->h_min), this->t_stop - t, h );

         // Don't overshoot the final time ...
         h = select( not_done & ((t + h) > this->t_stop), this->t_stop - t, h );

         ++iter;
         if ( this->max_iters && iter > this->max_iters )
         {
            ierr = ERR_TOO_MUCH_WORK;
            //printf("(iter > max_iters)\n");
            break;
         }
      }

      return ierr;

      #undef t
      #undef h
      #undef iter
      #undef nst
   }
};

template <typename Functor, typename RHSptr>
void simd_rk_driver ( const int numProblems, const double *u_in, const double t_stop, const Functor& func, const RHSptr rhs_func, const ckdata_t *RESTRICT ck )
{
   const int kk = ck->n_species;
   const int neq = kk+1;
   const double p = func.getPressure();

   typedef typename VCL_TypeSelector<double,4>::value_type SimdType;
   typedef typename VCL_MaskSelector<SimdType>::mask_type MaskType;
   const int VectorLength = VCL_Length<SimdType>::length;

   printf("Instruction Set= %d %s %d %s\n", INSTRSET, typeid(SimdType).name(), VectorLength, typeid( MaskType).name());

   VectorType<double,Alignment> scalar_out( neq * numProblems );
   VectorType<double,Alignment> vector_out( neq * numProblems );

   rk_t rk;

   rk_create (&rk, neq);

   rk.max_iters = 1000;
   rk.min_iters = 1;

   int lenrwk = rk_lenrwk (&rk);
   VectorType<double,Alignment> rwk( VectorLength * lenrwk );
   VectorType<double,Alignment> u( VectorLength * neq );

   int nst = 0, nit = 0, nfe = 0;

   auto scalar_solver = [&](const int i, VectorType<double,Alignment>& out, rk_counters_t *counters)
   {
      for (int k = 0; k < neq; ++k)
         u[k] = u_in[ i*neq + k ];

      const double T0 = u_in[ i*neq + getTempIndex(neq) ];

      double t = 0, h = 0;

      rk_init (&rk, t, t_stop);

      double t_begin = WallClock();

      int ierr = rk_solve (&rk, &t, &h, counters, u.getPointer(), rwk.getPointer(), rhs_func, (void*)&func);
      if (ierr != RK_SUCCESS)
      {
         fprintf(stderr,"%d: rk_solve error %d %d %d\n", i, ierr, counters->niters, rk.max_iters);
         exit(-1);
      }

      double t_end = WallClock();

      const int _nst = counters->nsteps;
      const int _nit = counters->niters;
      const int _nfe = _nit * 6;

      nst += _nst;
      nit += _nit;
      nfe += _nfe;

      for (int k = 0; k < neq; ++k)
         out[i*neq + k] = u[k];

      if (i % 1 == 0)
         printf("%d: %d %d %d %e %e %f\n", i, _nst, _nit, _nfe, u[ getTempIndex(neq) ], T0, 1000*(t_end-t_begin));
   };

   double time_scalar = WallClock();

   for (int i = 0; i < numProblems; ++i)
   {
      rk_counters_t counters;
      scalar_solver(i, scalar_out, &counters);
      //if (i % 10 == 0)
      //   printf("%d: %d %d %f\n", i, counters.nsteps, counters.niters, scalar_out[i*neq+getTempIndex(neq)] );
   }

   rk_destroy(&rk);

   time_scalar = WallClock() - time_scalar;

   double time_vector = WallClock();

   simd_cklib_functor<SimdType> simd_func( ck, p );
   //SimdERKSolverType<SimdType> simd_solver( neq );
   typedef SimdERKSolverType<SimdType> SimdSolverType;
   SimdSolverType simd_solver( neq );
   simd_solver.max_iters=200;

   for (int i0 = 0; i0 < numProblems; i0 += VectorLength)
   {
      if ( i0 + VectorLength <= numProblems )
      {
         for (int i = 0; i < VectorLength; ++i)
         {
            const double *ui = u_in + neq*(i0+i);
            for (int j = 0; j < neq; ++j)
               u[j*VectorLength+i] = ui[j];
         }

         SimdType *v_u = (SimdType *) u.getPointer();

         SimdType T0 = v_u[ getTempIndex(neq) ];

         SimdType t(0), h(0);
         typename SimdSolverType::CountersType counters;

         simd_solver.init( 0.0, t_stop );
         int ierr = simd_solver.solve ( &t, &h, &counters, v_u, simd_func );
         if (ierr != ERR_SUCCESS)
         {
            fprintf(stderr,"%d: simd_solver error %d %s\n", i0, ierr, GetErrorString(ierr));
            if (ierr == ERR_TOO_MUCH_WORK)
               fprintf(stderr,"--: simd_solver nit= %d %s\n", counters.nit, toString(counters.nst).c_str());
            exit(-1);
         }

         for (int i = 0; i < VectorLength; ++i)
         {
            double *v_out = vector_out.getPointer() + neq*(i0+i);
            for (int j = 0; j < neq; ++j)
               v_out[j] = u[j*VectorLength+i];
         }

         printf("i0: %d %s %d %s %s\n", i0, toString(counters.nst).c_str(), counters.nit, toString( v_u[getTempIndex(neq)] ).c_str(), toString(T0).c_str());
      }
      else
      {
         rk_counters_t counters;
         for (int i = i0; i < numProblems; ++i)
            scalar_solver( i, vector_out, &counters );
      }
   }

   time_vector = WallClock() - time_vector;

   printf("SIMD timer: %f %f %.1f\n", 1000.*time_vector, 1000.*time_scalar, time_scalar/time_vector);

   {
      double err2 = 0, ref2 = 0;
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
   {
      double err2 = 0, ref2 = 0;
      for (int k = 0; k < kk; ++k)
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
