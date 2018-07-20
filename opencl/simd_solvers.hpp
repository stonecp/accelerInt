#ifndef __simd_solvers_h
#define __simd_solvers_h

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

#include "simd_dtypes.hpp"
#include "simd_cklib.hpp"

namespace SIMD
{

enum SolverTags { RKF45 = 1, ROS4 };

enum ErrorFlags { ERR_SUCCESS         = 0,
                  ERR_TOO_MUCH_WORK   = 1,
                  ERR_TDIST_TOO_SMALL = 2,
                  ERR_HIN_MAX_ITERS   = 3 };

ErrorFlags GetErrorInt ( const int errflag )
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
   typedef typename SIMD_MaskSelector< ValueType >::mask_type   MaskType;

   enum { vector_length = SIMD_Length<ValueType>::length,
          vlen = vector_length };

   typedef int64_t BaseIndexType;
   typedef typename SIMD_TypeSelector<BaseIndexType,vlen>::value_type IndexType;

   struct CountersType
   {
      int nit;
      IndexType nst;
      IndexType nfe;
      IndexType nje;
      IndexType nlu;
      IndexType nni;
      IndexType errflag;

      CountersType(void)
         : nit(0),
           nst(0),
           nfe(0),
           nje(0),
           nlu(0),
           nni(0),
           errflag(ERR_SUCCESS)
      {}
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
        min_iters(1), max_iters(1000),
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

   /*ValueType wnorm (const ValueType *RESTRICT x, const ValueType *RESTRICT y)
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
   }*/

   template <class Functor>
   int hin ( const ValueType t, ValueType *h0, ValueType *RESTRICT y, ValueType *RESTRICT rwk, Functor& func) const
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
      *h0 = fmax(*h0, hlb);
      *h0 = fmin(*h0, hub);

      //printf("h0=%e, hlb=%e, hub=%e\n", h0, hlb, hub);

      return ierr;
   }

   ValueType getErrorWeight ( const int k, const ValueType *RESTRICT y ) const
   {
      //return (this->s_rtol * abs( y[k] )) + this->s_atol;
      return 1.0/((this->s_rtol * fabs( y[k] )) + this->s_atol);
   }

   ValueType wnorm (const ValueType *RESTRICT x, const ValueType *RESTRICT y) const
   {
      const int neq = this->neq;
      ValueType sum(0);
      for (int k = 0; k < neq; k++)
      {
         //ValueType prod = x[k] / this->getErrorWeight( k, y );
         ValueType prod = x[k] * this->getErrorWeight( k, y );
         sum += (prod*prod);
      }

      return sqrt(sum / (double)neq);
   }
   inline void dzero ( const int len, ValueType *RESTRICT x) const
   {
      const ValueType zero(0);
      for (int k = 0; k < len; ++k)
         x[k] = zero;
   }

   inline void dset ( const int len, ValueType *RESTRICT x, const double sval) const
   {
      const ValueType val(sval);
      for (int k = 0; k < len; ++k)
         x[k] = sval;
   }

   inline void dcopy ( const int len, const ValueType *RESTRICT src, ValueType *RESTRICT dst )
   {
      for (int k = 0; k < len; ++k)
         dst[k] = src[k];
   }
   /*inline void dcopy_if (const int len, const MaskType &mask, const __global __ValueType src[], __global __ValueType dst[])
   {
      for (int k = 0; k < len; ++k)
         dst[k] = if_then_else (mask, src[k], dst[k]);
   }*/

   /*inline void dscal (const int len, const double alpha, ValueType *RESTRICT y ) const
   {
      // Alpha is scalar type ... and can be easily checked.
      if (alpha == 1.0)
         return;
      else
      {
         #pragma ivdep
         for (int k = 0; k < len; ++k)
            y[(k)] *= alpha;
      }
   }*/
   inline void dscal (const int len, const ValueType& alpha, ValueType *RESTRICT y ) const
   {
      #pragma ivdep
      for (int k = 0; k < len; ++k)
         y[(k)] *= alpha;
   }

   inline void daxpy1 ( const int len, const double alpha, const ValueType *RESTRICT x, ValueType *RESTRICT y )
   {
      // Alpha is scalar type ... and can be easily checked.
      if (alpha == 1.0)
      {
         for (int k = 0; k < len; ++k)
            y[k] += x[k];
      }
      else if (alpha == -1.0)
      {
         for (int k = 0; k < len; ++k)
            y[k] -= x[k];
      }
      else if (alpha != 0.0)
      {
         for (int k = 0; k < len; ++k)
            y[k] += alpha * x[k];
      }
   }

   inline void daxpy ( const int len, const ValueType& alpha, const ValueType *RESTRICT x, ValueType *RESTRICT y ) const
   {
      // Alpha is vector type ... tedious to switch.
      for (int k = 0; k < len; ++k)
         y[k] += alpha * x[k];
   }

   // SIMD version.
   template <typename _ValueType>
   typename std::enable_if< (SIMD_isScalar<_ValueType>::value == 0), IndexType >::type
      lu_pivot ( const int k, const int n, _ValueType *RESTRICT A ) const
   {
      enum { nelems = vector_length };

      alignas(Alignment) BaseIndexType all_pivk[nelems];

      ValueType *A_k = A + (k*n); // pointer to this column

      /* find pivot row number */
      for (int el = 0; el < nelems; ++el)
      {
         int pivk = k;
         double Akp = A_k[pivk].extract(el);
         for (int i = k+1; i < n; ++i)
         {
            //const double Aki = __read_from( A_k[__getIndex(i)], el);
            double Aki = A_k[i].extract(el);
            if ( std::fabs(Aki) > std::fabs(Akp) )
            {
               pivk = i;
               Akp = Aki;
            }
         }

         // Test for singular value ...
         if (Akp == 0.0)
         {
            //ierr = (k+1);
            //printf("Singular value %d %d\n", k, el);
            break;
         }

         /* swap a(k,1:N) and a(piv,1:N) if necessary */
         if ( pivk != k )
         {
            ValueType *A_i = A; // pointer to the first column
            for (int i = 0; i < n; ++i, A_i += (n))
            {
               const double Aik = A_i[k   ].extract(el);
               const double Aip = A_i[pivk].extract(el);
               A_i[k   ].insert( el, Aip );
               A_i[pivk].insert( el, Aik );
            }
         }

         all_pivk[el] = pivk;

      } // End scalar section

      return IndexType().load_a( all_pivk );
   }

   // Scalar version.
   template <typename _ValueType>
   typename std::enable_if< SIMD_isScalar<_ValueType>::value, IndexType >::type
      lu_pivot ( const int k, const int n, _ValueType *RESTRICT A ) const
   {
      /* find pivot row number */

      ValueType *A_k = A + (k*n); // pointer to this column

      int pivk = k;
      double Akp = A_k[pivk];
      for (int i = k+1; i < n; ++i)
      {
         //const double Aki = __read_from( A_k[__getIndex(i)], el);
         double Aki = A_k[i];
         if ( std::fabs(Aki) > std::fabs(Akp) )
         {
            pivk = i;
            Akp = Aki;
         }
      }

      // Test for singular value ...
      if (Akp == 0.0)
      {
         return (k+1);
         //printf("Singular value %d %d\n", k, el);
      }

      /* swap a(k,1:N) and a(piv,1:N) if necessary */
      if ( pivk != k )
      {
         ValueType *A_i = A; // pointer to the first column
         for (int i = 0; i < n; ++i, A_i += (n))
         {
            const double Aik = A_i[k   ];
            const double Aip = A_i[pivk];
            A_i[k   ] = Aip;
            A_i[pivk] = Aik;
         }
      }

      return pivk;
   }

   int ludec ( const int n, ValueType *RESTRICT A, IndexType *RESTRICT ipiv ) const
   {
      int ierr = ERR_SUCCESS;

      const int nelems = this->vector_length;

      //typedef typename decltype( ipiv[0][0] ) IndexBaseType;
      //typedef int64_t BaseIndexType;
      //alignas(Alignment) BaseIndexType all_pivk[nelems];

      /* k-th elimination step number */
      for (int k = 0; k < n; ++k)
      {
        ValueType *A_k = A + (k*n); // pointer to this column

        ipiv[k] = lu_pivot( k, n, A );

        /* Scale the elements below the diagonal in
         * column k by 1.0/a(k,k). After the above swap
         * a(k,k) holds the pivot element. This scaling
         * stores the pivot row multipliers a(i,k)/a(k,k)
         * in a(i,k), i=k+1, ..., M-1.
         */
        const ValueType mult = 1.0 / A_k[k];
        for (int i = k+1; i < n; ++i)
          A_k[i] *= mult;

        /* row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., m-1 */
        /* row k is the pivot row after swapping with row l.      */
        /* The computation is done one column at a time,          */
        /* column j=k+1, ..., n-1.                                */

        for (int j = k+1; j < n; ++j)
        {
          ValueType *A_j = A + (j*n);
          const ValueType a_kj = A_j[k];

          /* a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j)  */
          /* a_kj = a(k,j), col_k[i] = - a(i,k)/a(k,k) */
          //if (any(a_kj != 0.0)) {
            for (int i = k+1; i < n; ++i) {
              A_j[i] -= a_kj * A_k[i];
            }
          //}
        }
      }

      return ierr;
      //if (ierr)
      //{
      //  fprintf(stderr,"Singular pivot j=%d\n", ierr-1);
      //  exit(-1);
      //}
   }

   template <typename _ValueType>
      typename std::enable_if< SIMD_isScalar<_ValueType>::value, void >::type
   lu_permute ( const int k, const IndexType& pivk, _ValueType *RESTRICT b) const
   {
      if ( pivk != IndexType( k)  )
      {
         const double bk = b[k   ];
         const double bp = b[pivk];
         b[k   ] = bp;
         b[pivk] = bk;
      }
   }

   // SIMD version.
   template <typename _ValueType>
      typename std::enable_if< (SIMD_isScalar<_ValueType>::value == 0), void >::type
   lu_permute ( const int k, const IndexType& ipiv, _ValueType *RESTRICT b) const
   {
      if ( any( ipiv != IndexType( k ) ) )
      {
         for (int el = 0; el < vlen; ++el)
         {
            const int pivk = ipiv.extract(el);
            if ( pivk != k )
            {
               const double bk = b[k   ].extract( el);
               const double bp = b[pivk].extract( el);
               b[k   ].insert( el, bp );
               b[pivk].insert( el, bk );
            }
         }
      }
   }

   void lusol ( const int n, ValueType *RESTRICT A, IndexType *RESTRICT ipiv, ValueType *RESTRICT b) const
   {
      /* Permute b, based on pivot information in p */
      for (int k = 0; k < n; ++k)
      {
         lu_permute (k, ipiv[k], b );
      }

      /* Solve Ly = b, store solution y in b */
      for (int k = 0; k < n-1; ++k)
      {
         ValueType *A_k = A + (k*n);
         const ValueType bk = b[k];
         for (int i = k+1; i < n; ++i)
            b[i] -= A_k[i] * bk;
      }
      /* Solve Ux = y, store solution x in b */
      for (int k = n-1; k > 0; --k)
      {
         ValueType *A_k = A + (k*n);
         b[k] /= A_k[k];
         const ValueType bk = b[k];
         for (int i = 0; i < k; ++i)
            b[i] -= A_k[i] * bk;
      }
      b[0] /= A[0];
   }

   template <class Functor>
   void fdjac ( const ValueType& tcur,
                const ValueType& hcur,
                      ValueType *RESTRICT y,
                      ValueType *RESTRICT fy,
                      ValueType *RESTRICT Jy,
                      Functor& func ) const
   {
      const int neq = this->neq;

      // Norm of fy(t) ...
      const ValueType fnorm = this->wnorm( fy, y );

      // Safety factors ...
      const double sround = std::sqrt( this->uround() );
      ValueType r0 = (1000. * this->uround() * neq) * (hcur * fnorm);
      //if (r0 == 0.) r0 = 1.;
      r0 = select( (r0 == 0.0), ValueType(1), r0 );

      // Build each column vector ...
      for (int j = 0; j < neq; ++j)
      {
         const ValueType ysav = y[j];
         const ValueType ewtj = this->getErrorWeight( j, y);
         //const ValueType dely = fmax( sround * abs(ysav), r0 * ewtj );
         const ValueType dely = fmax( sround * fabs(ysav), r0 / ewtj );
         y[j] += dely;

         ValueType *jcol = Jy + (j*neq);

         func ( neq, tcur, y, jcol );

         const ValueType delyi = 1. / dely;
         for (int i = 0; i < neq; ++i)
            jcol[i] = (jcol[i] - fy[i]) * delyi;

         y[j] = ysav;
      }
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

   size_t get_lenrwk(const int _neq) const { return (8 * _neq); }
   size_t get_lenrwk(void) const { return (8 * this->neq); }

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
         rwk[k] = fabs(r5 - r4);

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

      MaskType not_done = fabs(t - this->t_stop) > ValueType( this->t_round );

      while (  any(not_done) )
      {
         ValueType *ytry = rwk + this->neq*7;

         // Take a trial step over h_cur ...
         this->oneStep ( h, y, ytry, rwk, func );

         const ValueType errmin (1.0e-20);
         const ValueType err = this->wnorm ( rwk, y );
         ValueType herr = fmax(/*1.0e-20*/ errmin, err );

         // Is there error acceptable?
         MaskType accept = ((herr <= 1.0) | (h <= this->h_min)) & not_done;

         // update solution ...
         if ( any(accept) )
         {
            t   = select ( accept, t + h, t   );
            nst = select ( accept, nst+1, nst );

            for (int k = 0; k < this->neq; k++)
               y[k] = select( accept, ytry[k], y[k] );

            not_done = fabs(t - this->t_stop) > this->t_round;
         }

         ValueType fact = sqrt( sqrt(1.0 / herr) ) * (0.840896415);

         // Restrict the rate of change in dt
         fact = fmax(fact, 1.0 / this->adaption_limit);
         fact = fmin(fact,       this->adaption_limit);

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
         h = select( accept & ( fabs((t + h) - this->t_stop) < this->h_min), this->t_stop - t, h );

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

template <typename _ValueType>//, SolverTags SolverTag >
struct SimdRosSolverType : public CommonSolverType<_ValueType>
{
   enum SolverTags { Ros3, Ros4, Rodas3, Rodas4 };

   enum { _maxStages = 6 };

   typedef CommonSolverType<_ValueType> Parent;

   enum { vlen = Parent::vlen };

   typedef typename Parent::ValueType    ValueType;
   typedef typename Parent::MaskType     MaskType;
   typedef typename Parent::CountersType CountersType;
   typedef typename Parent::IndexType    IndexType;

   VectorType<ValueType> m_rwk;
   VectorType<IndexType> m_iwk;

   SolverTags solverTag;
   int numStages;
   int ELO;
   double A[_maxStages*(_maxStages-1)/2];
   double C[_maxStages*(_maxStages-1)/2];
   int newFunc[_maxStages];
   double E[_maxStages];
   double M[_maxStages];
   double alpha[_maxStages];
   double gamma[_maxStages];

   SimdRosSolverType( const int _neq, const SolverTags& _tag = Ros4 )
      : Parent(_neq),
        solverTag( _tag )
   {
      if (solverTag == Ros4)
         this->setRos4();
      else
      {
         fprintf(stderr,"Invalid solverTag = %d\n", solverTag);
         exit(1);
      }

      // lenrwk function uses subspace lengths.
      m_rwk.resize( this->get_lenrwk() );
      m_iwk.resize( this->get_leniwk() );
   }

   size_t get_lenrwk(const int _neq) const
   {
      int lenrwk = 0;
      lenrwk +=  _neq;		// fy
      lenrwk +=  _neq;		// ynew & yerr
      lenrwk += (_neq * _neq);	// Jy
    //lenrwk +=  _neq;		// ewt
      lenrwk +=  _neq * this->numStages;// ktmp

      return lenrwk;
   }

   size_t get_leniwk(const int _neq) const
   {
      return (_neq); // ipiv
   }

   size_t get_lenrwk(void) const { return this->get_lenrwk(this->neq); }
   size_t get_leniwk(void) const { return this->get_leniwk(this->neq); }

   // 4th/3rd-order L-stable Rosenbrock method with 4 stages.
   // -- E. Hairer and G. Wanner, "Solving ordinary differential equations II:
   //    stiff and differential-algebraic problems," Springer series in
   //    computational mathematics, Springer-Verlag (1990).
   void setRos4 (void)
   {
      this->solverTag = Ros4;
      this->numStages = 4;
      this->ELO = 4;

      // A and C are strictly lower-triangular matrices in row-major order!!!!
      // -- A(i,j) = [(i)*(i-1)/2 + j] ... A(1,0) = A[0], A(2,0) = A[1]
      this->A[0] = 2.0;
      this->A[1] = 1.867943637803922;
      this->A[2] = 0.2344449711399156;
      this->A[3] = this->A[1];
      this->A[4] = this->A[2];
      this->A[5] = 0.0;

      this->C[0] =-7.137615036412310;
      this->C[1] = 2.580708087951457;
      this->C[2] = 0.6515950076447975;
      this->C[3] =-2.137148994382534;
      this->C[4] =-0.3214669691237626;
      this->C[5] =-0.6949742501781779;

      // Does the stage[i] need a new function eval or can it reuse the
      // prior one from stage[i-1]?
      this->newFunc[0] = 1;
      this->newFunc[1] = 1;
      this->newFunc[2] = 1;
      this->newFunc[3] = 0;

      // M_i = Coefficients for new step solution
      this->M[0] = 2.255570073418735;
      this->M[1] = 0.2870493262186792;
      this->M[2] = 0.4353179431840180;
      this->M[3] = 1.093502252409163;

      // E_i = Coefficients for error estimator
      this->E[0] =-0.2815431932141155;
      this->E[1] =-0.07276199124938920;
      this->E[2] =-0.1082196201495311;
      this->E[3] =-1.093502252409163;

      // Y( T + h*alpha_i )
      this->alpha[0] = 0.0;
      this->alpha[1] = 1.14564;
      this->alpha[2] = 0.65521686381559;
      this->alpha[3] = this->alpha[2];

      // gamma_i = \Sum_j  gamma_{i,j}
      this->gamma[0] = 0.57282;
      this->gamma[1] =-1.769193891319233;
      this->gamma[2] = 0.7592633437920482;
      this->gamma[3] =-0.104902108710045;
   }

   // ROS internal routines ...
   template <class Functor, typename Counters>
   int solve ( ValueType *tcur, ValueType *hcur, Counters* counters, ValueType y[], Functor& func)
   {
      int ierr = ERR_SUCCESS;

      const int neq = this->neq;

      #define nst (counters->nst)
      #define nfe (counters->nfe)
      #define nje (counters->nje)
      #define nlu (counters->nlu)
      #define iter (counters->nit)
      #define h (*hcur)
      #define t (*tcur)
      #define A(_i,_j) (this->A[ (((_i)-1)*(_i))/2 + (_j) ] )
      #define C(_i,_j) (this->C[ (((_i)-1)*(_i))/2 + (_j) ] )

      ValueType *rwk = m_rwk.getPointer();
      IndexType *iwk = m_iwk.getPointer();

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

      // Zero the counters ...
      nst = 0;
      nfe = 0;
      //nlu = 0;
      //nje = 0;
      iter = 0;

      // Set the work arrays ...
      ValueType *RESTRICT fy   = rwk;
      ValueType *RESTRICT ynew = fy +   (neq);
      ValueType *RESTRICT Jy   = ynew + (neq);
      ValueType *RESTRICT ktmp = Jy +   (neq*neq);
      ValueType *RESTRICT yerr = ynew;
      //__global double *ewt  = &Jy[neq*neq];

      MaskType not_done = fabs(t - this->t_stop) > ValueType( this->t_round );
      while ( any(not_done) )
      {
         // Compute the RHS and Jacobian matrix.
         func (neq, t, y, fy);
         nfe++;

         //if (jac == NULL)
         {
            this->fdjac ( t, h, y, fy, Jy, func );
            nfe += neq;
         }
         //else
         //{
         //   jac (neq, t, y, Jy, user_data);
         //}

         //nje++;

         // Construct iteration matrix J' := 1/(gamma*h) - J
         {
            const ValueType one_hgamma = 1.0 / ( h * this->gamma[0] );

            for (int j = 0; j < neq; ++j)
            {
               ValueType *jcol = Jy + (j*neq);

               for (int i = 0; i < neq; ++i)
                  jcol[i] = -jcol[i];

               jcol[j] += one_hgamma;
            }
         }

         // Factorization J'
         this->ludec( neq, Jy, iwk ); // simd variant
         //nlu++;

         for (int s = 0; s < this->numStages; s++)
         {
            // Compute the function at this stage ...
            if (s == 0)
            {
               //func (neq, y, fy.getPointer());
            }
            else if ( this->newFunc[s] )
            {
               this->dcopy (neq, y, ynew);

               for (int j = 0; j < s; ++j)
               {
                  const double Asj = A(s,j);
                  //if (Asj != 0.0)
                  {
                     //printf("Asj = %f %d %d\n", Asj, s, j);
                     ValueType *k_j = ktmp + (j*neq);

                     this->daxpy1 ( neq, Asj, k_j, ynew );
                  }
               }

               func (neq, t, ynew, fy);
               nfe++;

               //printf("newF=%d\n", s);
               //for (int k = 0; k < neq; ++k)
               //   printf("ynew[%d] = %e %e\n", k, ynew[k], fy[k]);
            }

            //printf("stage=%d\n", s);
            //for (int k = 0; k < neq; ++k)
            //   printf("fy[%d] = %e\n", k, fy[k]);

            // Build the sub-space vector K
            ValueType *k_s = ktmp + (s*neq);
            this->dcopy (neq, fy, k_s);

            for (int j = 0; j < s; j++)
            {
               //if (C(s,j) != 0.0)
               {
                  const ValueType hCsj = C(s,j) / h;
                  //printf("C/h = %f %d %d\n", hCsj, s, j);

                  ValueType *k_j = ktmp + (j*neq);
                  this->daxpy (neq, hCsj, k_j, k_s);
               }
            }

            //printf("k before=%d\n", s);
            //for (int k = 0; k < neq; ++k)
            //   printf("k[%d] = %e\n", k, ks[k]);

            // Solve the current stage ..
            this->lusol (neq, Jy, iwk, k_s);

            //printf("k after=%d\n", s);
            //for (int k = 0; k < neq; ++k)
            //   printf("k[%d] = %e\n", k, ks[k]);
         }

         // Compute the error estimation of the trial solution
         this->dzero (neq, yerr);

         for (int j = 0; j < this->numStages; ++j)
         {
            //if (this->E[j] != 0.0)
            {
               ValueType *k_j = ktmp + (j*neq);
               this->daxpy1 (neq, this->E[j], k_j, yerr);
            }
         }

         const ValueType herr = fmax( 1.0e-20, this->wnorm (yerr, y));

         // Is there error acceptable?
         MaskType accept = ((herr <= 1.0) | (h <= this->h_min)) & not_done;

         // update solution ...
         if ( any(accept) )
         {
            t   = select ( accept, t + h, t   );
            nst = select ( accept, nst+1, nst );

            not_done = fabs(t - this->t_stop) > this->t_round;

            // Need to actually compute the new solution since it was delayed from above.
            this->dcopy (neq, y, ynew);
            for (int j = 0; j < this->numStages; ++j)
            {
               //if (this->M[j] != 0.0)
               {
                  ValueType *k_j = ktmp + (j*neq);
                  this->daxpy1 (neq, this->M[j], k_j, ynew);
               }
            }

            for (int k = 0; k < this->neq; k++)
               y[k] = select( accept, ynew[k], y[k] );
         }

         ValueType fact = 0.9 * pow( 1.0 / herr, (1.0/ this->ELO));

         // Restrict the rate of change in dt
         fact = fmax(fact, 1.0 / this->adaption_limit);
         fact = fmin(fact,       this->adaption_limit);

#if defined(VERBOSE) && (VERBOSE > 0)
         if (iter % VERBOSE == 0)
         {
            std::cout << "iter= " << iter;
            std::cout << " accept= " << accept;
            std::cout << " done= " << not_done;
            std::cout << " t= " << t;
            std::cout << " h= " << h;
            std::cout << " fact= " << fact;
            std::cout << " T= " << y[getTempIndex(neq)] << "\n";
         }
#endif

         // Apply grow/shrink factor for next step.
         h = select( not_done, h * fact, h);

         // Limit based on the upper/lower bounds
         h = fmin(h, this->h_max);
         h = fmax(h, this->h_min);

         // Stretch the final step if we're really close and we didn't just fail ...
         h = select( accept & ( fabs((t + h) - this->t_stop) < this->h_min), this->t_stop - t, h );

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

      #undef nst
      #undef nfe
      #undef nje
      #undef nlu
      #undef iter
      #undef h
      #undef t
      #undef A
      #undef C
   }

};

#define __matrix_index(_i,_j,_var) (this-> _var [(_i)][(_j)])
#define A_(_i,_j)     ( __matrix_index(_i,_j,A) )
#define Theta_(_i,_j) ( __matrix_index(_i,_j,Theta) )
#define Alpha_(_i,_j) ( __matrix_index(_i,_j,Alpha) )

template <typename _ValueType, bool _ReplicateSerial = false >
struct SimdSdirkSolverType : public CommonSolverType<_ValueType>
{
   enum { _maxStages = 5 };
   enum { _replicate = (_ReplicateSerial) ? 1 : 0 };

   typedef CommonSolverType<_ValueType> Parent;

   enum { vlen = Parent::vlen };

   typedef typename Parent::ValueType     ValueType;
   typedef typename Parent::MaskType      MaskType;
   typedef typename Parent::CountersType  CountersType;
   typedef typename Parent::IndexType     IndexType;
   typedef typename Parent::BaseIndexType BaseIndexType;

   typedef sdirk_solverTags solverTags;

   VectorType<ValueType> m_rwk;
   VectorType<IndexType> m_iwk;

   solverTags solverTag;
   int numStages;
   int ELO;
   double A[_maxStages][_maxStages];
   double B[_maxStages];
   double Bhat[_maxStages];
   double C[_maxStages];
   double D[_maxStages];
   double E[_maxStages];
   double Theta[_maxStages][_maxStages];
   double Alpha[_maxStages][_maxStages];
   double gamma;

   SimdSdirkSolverType( const int _neq, const solverTags _tag = S4a )
      : Parent(_neq),
        solverTag( _tag)
   {
      this->adaption_limit = 5.0; // This is higher than the default to match sdirk.c!!!

      if (solverTag == S4a)
         this->set_S4a();
      else
      {
         fprintf(stderr,"Invalid Sdirk solverTag = %d\n", solverTag);
         exit(1);
      }

      // lenrwk function uses subspace lengths.
      m_rwk.resize( this->get_lenrwk() );
      m_iwk.resize( this->get_leniwk() );
   }

   size_t get_lenrwk(const int _neq) const
   {
      int lenrwk = 0;
      lenrwk +=  _neq;			 // fy
      lenrwk +=  _neq;			 // del & yerr
      lenrwk += (_neq * _neq);		 // Jy
      lenrwk += (_neq * _neq);		 // M
      lenrwk += (_neq * this->numStages);// z
      lenrwk +=  _neq;			 // g

      if ( _replicate ) lenrwk += (_neq*_neq); // An extra matrix.

      return lenrwk;
   }

   size_t get_leniwk(const int _neq) const
   {
      size_t leniwk = (_neq); // ipiv
      if ( _replicate ) leniwk += (_neq); // copy of ipiv
      return leniwk;
   }

   size_t get_lenrwk(void) const { return this->get_lenrwk(this->neq); }
   size_t get_leniwk(void) const { return this->get_leniwk(this->neq); }

   // SDIRK internal routines ...

   // 4th/3rd-order L-stable SDIRK method with 5 stages.
   // -- E. Hairer and G. Wanner, "Solving ordinary differential equations II:
   //    stiff and differential-algebraic problems," Springer series in
   //    computational mathematics, Springer-Verlag (1990).
   void set_S4a (void)
   {
      this->solverTag = S4a;
      this->numStages = 5;
      this->ELO = 4;

      // Constant diagonal
      this->gamma = 8.0 / 30.0; //0.2666666666666666666666666666666667;

      // A and C are lower-triangular matrices in column-major order!!!!
      // -- A(i,j) = [(i)*(i+1)/2 + j] ... A(1,0) = A[0], A(2,0) = A[1]
      for (int i = 0; i < _maxStages; ++i)
         for (int j = 0; j < _maxStages; ++j)
            A_(i,j) = 0.0;
      A_(0,0) = this->gamma;
      A_(1,0) = 0.5;
      A_(1,1) = this->gamma;
      A_(2,0) = 0.3541539528432732316227461858529820;
      A_(2,1) =-0.05415395284327323162274618585298197;
      A_(2,2) = this->gamma;
      A_(3,0) = 0.08515494131138652076337791881433756;
      A_(3,1) =-0.06484332287891555171683963466229754;
      A_(3,2) = 0.07915325296404206392428857585141242;
      A_(3,3) = this->gamma;
      A_(4,0) = 2.100115700566932777970612055999074;
      A_(4,1) =-0.7677800284445976813343102185062276;
      A_(4,2) = 2.399816361080026398094746205273880;
      A_(4,3) =-2.998818699869028161397714709433394;
      A_(4,4) = this->gamma;

      this->B[0]    = 2.100115700566932777970612055999074;
      this->B[1]    =-0.7677800284445976813343102185062276;
      this->B[2]    = 2.399816361080026398094746205273880;
      this->B[3]    =-2.998818699869028161397714709433394;
      this->B[4]    = this->gamma;

      this->Bhat[0] = 2.885264204387193942183851612883390;
      this->Bhat[1] =-0.1458793482962771337341223443218041;
      this->Bhat[2] = 2.390008682465139866479830743628554;
      this->Bhat[3] =-4.129393538556056674929560012190140;
      this->Bhat[4] = 0.;

      this->C[0]    = 8.0  / 30.0; //0.2666666666666666666666666666666667;
      this->C[1]    = 23.0 / 30.0; // 0.7666666666666666666666666666666667;
      this->C[2]    = 17.0 / 30.0; // 0.5666666666666666666666666666666667;
      this->C[3]    = 0.3661315380631796996374935266701191;
      this->C[4]    = 1.;

      // Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      this->D[0] = 0.;
      this->D[1] = 0.;
      this->D[2] = 0.;
      this->D[3] = 0.;
      this->D[4] = 1.;

      // Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      this->E[0] =-0.6804000050475287124787034884002302;
      this->E[1] = 1.558961944525217193393931795738823;
      this->E[2] =-13.55893003128907927748632408763868;
      this->E[3] = 15.48522576958521253098585004571302;
      this->E[4] = 1.;

      // h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      for (int i = 0; i < _maxStages; ++i)
         for (int j = 0; j < _maxStages; ++j)
            Theta_(i,j) = 0.0;
      Theta_(1,0) = 1.875;
      Theta_(2,0) = 1.708847304091539528432732316227462;
      Theta_(2,1) =-0.2030773231622746185852981969486824;
      Theta_(3,0) = 0.2680325578937783958847157206823118;
      Theta_(3,1) =-0.1828840955527181631794050728644549;
      Theta_(3,2) = 0.2968246986151577397160821594427966;
      Theta_(4,0) = 0.9096171815241460655379433581446771;
      Theta_(4,1) =-3.108254967778352416114774430509465;
      Theta_(4,2) = 12.33727431701306195581826123274001;
      Theta_(4,3) =-11.24557012450885560524143016037523;

      // Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      for (int i = 0; i < _maxStages; ++i)
         for (int j = 0; j < _maxStages; ++j)
            Alpha_(i,j) = 0.0;
      Alpha_(1,0) = 2.875000000000000000000000000000000;
      Alpha_(2,0) = 0.8500000000000000000000000000000000;
      Alpha_(2,1) = 0.4434782608695652173913043478260870;
      Alpha_(3,0) = 0.7352046091658870564637910527807370;
      Alpha_(3,1) =-0.09525565003057343527941920657462074;
      Alpha_(3,2) = 0.4290111305453813852259481840631738;
      Alpha_(4,0) =-16.10898993405067684831655675112808;
      Alpha_(4,1) = 6.559571569643355712998131800797873;
      Alpha_(4,2) =-15.90772144271326504260996815012482;
      Alpha_(4,3) = 25.34908987169226073668861694892683;
   }

   #define nst (counters->nst)
   #define nfe (counters->nfe)
   #define nje (counters->nje)
   #define nlu (counters->nlu)
   #define nni (counters->nni)
   #define iter (counters->nit)
   #define h (*hcur)
   #define t (*tcur)

   template <class Functor, typename Counters>
   int solve_anyall ( ValueType *tcur, ValueType *hcur, Counters* counters, ValueType y[], Functor& func)
   {
      int ierr = ERR_SUCCESS;

      const int neq = this->neq;

      const int InterpolateNewton  = 1;	// Start at zero (0) or interpolate a starting guess (1)
      const int MaxNewtonIterations= 8;	// Max number of newton iterations
      const double NewtonThetaMin  = 0.005; // Minimum convergence rate for the Newton Iteration (0.001)
      const double NewtonThetaMax  = 0.999; // Maximum residual drop acceptable
      const double NewtonTolerance = 0.03;  // Convergence criteria
      const double Qmax = 1.2;		// Max h-adaption to recycle M
      const double Qmin = 1.;		// Min ""

      ValueType *rwk = m_rwk.getPointer();
      IndexType *iwk = m_iwk.getPointer();

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

      // Zero the counters ...
      nst = 0;
      nfe = 0;
      nlu = 0;
      nje = 0;
      nni = 0;
      iter = 0;

      // Set the work arrays ...
      ValueType *RESTRICT fy   = rwk;
      ValueType *RESTRICT del  = fy + (neq);
      ValueType *RESTRICT Jy   = del + (neq);
      ValueType *RESTRICT M    = Jy + (neq*neq);
      ValueType *RESTRICT z    = M + (neq*neq);
      ValueType *RESTRICT g    = z + (neq*this->numStages);
      ValueType *RESTRICT yerr = del;

      bool ComputeJ = 1;
      bool ComputeM = 1;

      MaskType not_done = fabs(t - this->t_stop) > ValueType( this->t_round );
      while ( any(not_done) )
      {
         //std::cout << "Step: " << iter << nst << " " << ComputeM << " " << ComputeJ << not_done << t << "\n";

         // Construct the Iteration matrix ... if needed.
         if (ComputeM)
         {
            // Compute the Jacobian matrix or recycle an old one.
            if (ComputeJ)
            {
               //if (jac)
               //   jac (neq, t, y, Jy, user_data);
               //else
               {
                  // Compute the RHS ... the fd algorithm expects it.
                  //if (func)
                     func (neq, t, y, fy);
                  //else
                  //   cklib_callback (neq, t, y, fy, user_data);
                  //nfe++;

                  this->fdjac ( t, h, y, fy, Jy, func );
                  //nfe += neq;
                  nfe = select( not_done, nfe + (neq+1), nfe );
               }

               //nje++;
               nje = select( not_done, nje+1, nje );
            }

            // Compute M := 1/(gamma*h) - J
            const ValueType one_hgamma = 1.0 / (h * this->gamma);

            for (int j = 0; j < neq; ++j)
            {
               ValueType *RESTRICT Mcol = &M[(j*neq)];
               ValueType *RESTRICT Jcol = &Jy[(j*neq)];
               for (int i = 0; i < neq; ++i)
                  Mcol[(i)] = -Jcol[(i)];

               Mcol[(j)] += one_hgamma;
            }

            // Factorization M
            this->ludec( neq, M, iwk );
            //nlu++;
            nlu = select(not_done, nlu+1, nlu);
         }

         const BaseIndexType Status_None      = 0,
                             Status_Converged = 1,
                             Status_Diverged  = 2,
                             Status_Finished  = 3;
         IndexType Status = select( not_done, IndexType(Status_None), IndexType(Status_Finished) );

         MaskType  Accepted; // All are intialized inside of stage loop.
         ValueType HScalingFactor;
         ValueType NewtonTheta;

         for (int s = 0; s < this->numStages; s++)
         {
            // Initial the RK stage vectors Z_i and G.
            ValueType *z_s = &z[(s*neq)];
            this->dzero (neq, z_s);
            this->dzero (neq, g);

            // Compute the function at this stage ...
            if (s)
            {
               for (int j = 0; j < s; ++j)
               {
                  // G = \Sum_j Theta_i,j*Z_j = h * \Sum_j A_i,j*F(Z_j)
                  ValueType *z_j = &z[(j*neq)];
                  this->daxpy1 (neq, Theta_(s,j), z_j, g);

                  // Z_i = \Sum_j Alpha(i,j)*Z_j
                  if (InterpolateNewton)
                     this->daxpy1 (neq, Alpha_(s,j), z_j, z_s);
               }
            }

            // Solve the non-linear problem with the Newton-Raphson method.
            Status = select( not_done, IndexType(Status_None), IndexType(Status_Finished) );
            Accepted = false;
            HScalingFactor = 0.8;
            NewtonTheta = NewtonThetaMin;
            ValueType NewtonResidual;

            //for (int NewtonIter = 0; NewtonIter < MaxNewtonIterations; NewtonIter++, nni++)
            for (int NewtonIter = 0; NewtonIter < MaxNewtonIterations; NewtonIter++)
            {
               //nni++;
               nni = select( (Status == Status_None), nni+1, nni );

               // 1) Build the RHS of the root equation: F := G + (h*gamma)*f(y+Z_s) - Z_s
               for (int k = 0; k < neq; ++k)
                  del[(k)] = y[(k)] + z_s[(k)];

               func (neq, t, del, fy);
               //nfe++;
               nfe = select( (Status == Status_None), nfe+1, nfe );

               const ValueType hgamma = h * this->gamma;
               for (int k = 0; k < neq; ++k)
                  del[(k)] = g[(k)] + hgamma * fy[(k)] - z_s[(k)];
                //del[(k)] = g[(k)] - z_s[(k)] + hgamma * fy[(k)];

               // 2) Solve the linear problem: M*delz = (1/hgamma)F ... so scale F first.
               this->dscal (neq, 1.0/hgamma, del);
               this->lusol (neq, M, iwk, del);

               // 3) Check the convergence and convergence rate
               // 3.1) Compute the norm of the correction.
               const ValueType dnorm = this->wnorm ( del, y);

               //std::cout << "miter: " << iter <<" " << s << " "<< NewtonIter << dnorm << NewtonResidual << Status << "\n";

               // 3.2) If not the first iteration, estimate the rate.
               if (NewtonIter != 0)
               {
                  NewtonTheta = dnorm / NewtonResidual; // How much as the residual changed from last iteration.

                  MaskType isDiverging = (NewtonTheta >= NewtonThetaMax) && not_done;

                  ValueType ConvergenceRate = NewtonTheta / (1.0 - NewtonTheta); // Possible FLTEXP here.

                  //std::cout << " " << ConvergenceRate << isDiverging << NewtonTheta << "\n";

                  // If one is diverging, may just give up early all around.
                  // Use the status flag to change or hold h.
                  if ( any( isDiverging ) )
                  {
                     Status = select( isDiverging, Status_Diverged, Status );
                     //std::cout << "any(NewtonTheta >= NewtonThetaMax) " << isDiverging << "\n";
                     break;
                  }

                  // So, nobody is diverging. Test if any are converging too slowly.
                  // We'll shrink h and restart, too.
                  // Predict the error after Max iterations with the current rate:
                  // ... res * Theta^(ItersLeft/(1-Theta))
                  ValueType PredictedError = dnorm * pow( NewtonTheta,
                                       (MaxNewtonIterations-NewtonIter)/(1.0-NewtonTheta));

                  MaskType PredictedErrorTooLarge = ( PredictedError > NewtonTolerance ) && not_done;
                  if ( any( PredictedErrorTooLarge ) )
                  {
                     // Doubtful we'll converge, shrink h and try again.
                     ValueType QNewton = fmin(10.0, PredictedError / NewtonTolerance);
                     HScalingFactor = select( PredictedErrorTooLarge,
                                              0.9 * pow( QNewton, -1.0 / (1.0 + MaxNewtonIterations - NewtonIter)),
                                              1.0 );
                     //std::cout << "PredictedError > NewtonTolerance " << HScalingFactor << PredictedError << PredictedErrorTooLarge << "\n";
                     break;
                  }

                  // So, nobody is diverging and none are predicted to fail. Test for new convergence.
                  Accepted = (ConvergenceRate * dnorm < NewtonTolerance) || not(not_done);
                  //std::cout << "Accepted: " << Accepted << Status << "\n";
               }

               // Save the residual norm for the next iteration unless I'm already converged.
               MaskType UpdateSolution = MaskType(Status == Status_None) && not_done;
               NewtonResidual = select( UpdateSolution, dnorm, NewtonResidual );

               // 4) Update the solution if newly accepted: Z_s <- Z_s + delta
               ValueType OneOrZero = select( UpdateSolution, ValueType(1), ValueType(0) );
               this->daxpy (neq, OneOrZero, del, z_s);

               // Finally, test if everyone's finally finished.
               Status = select( Accepted, Status_Converged, Status );
               if ( all( Accepted || MaskType(Status == Status_Finished) ) )
                  break;
            }

            if ( any(!Accepted) )
            {
               //printf("Any(!Accepted) %d %d.\n", iter, s);
               ComputeJ = 0; // Jacobian is still valid
               ComputeM = 1; // M is invalid since h will change (perhaps) drastically.
               break;
               //return 0;
            }

         } // ... stages

         if ( all( Accepted || MaskType(Status == Status_Finished) ) )
         {
            // Compute the error estimation of the trial solution.
            this->dzero (neq, yerr);
            for (int j = 0; j < this->numStages; ++j)
            {
               ValueType *z_j = &z[(j*neq)];
               if (this->E[j] != 0.0)
                  this->daxpy1 (neq, this->E[j], z_j, yerr);
            }

            ValueType herr = fmax(1.0e-20, this->wnorm ( yerr, y));

            // Is the error acceptable?
            Accepted = ((herr <= 1.0) || (h <= this->h_min)) && not_done;
            if ( any( Accepted ) )
            {
               // If stiffly-accurate, Z_s with s := numStages, is the solution.
               // Else, sum the stage solutions: y_n+1 <- y_n + \Sum_j D_j * Z_j
               for (int j = 0; j < this->numStages; ++j)
               {
                  ValueType *z_j = &z[(j*neq)];
                  if (this->D[j] != 0.0)
                  {
                     // Only update solutions that were accepted.
                     ValueType ZeroOrDj = select( Accepted, ValueType(this->D[j]), 0.0 );
                     this->daxpy ( neq, ZeroOrDj, z_j, y );
                  }
               }

               t = select( Accepted, t + h, t );
               nst = select( Accepted, nst+1, nst );
            }

            not_done = fabs(t - this->t_stop) > ValueType( this->t_round );

            HScalingFactor = select( not_done, 0.9 * pow( 1.0 / herr, (1.0 / this->ELO)), ValueType(1) );

            // Reuse the Jacobian if the Newton Solver is converging fast enough
            // ... for everyone. (any-all)
            ComputeJ = any( NewtonTheta > NewtonThetaMin && not_done );

            // Don't refine if it's not a big step and we could reuse the M matrix.
            bool recycle_M = not(ComputeJ) and all( HScalingFactor >= Qmin && HScalingFactor <= Qmax );
            if (recycle_M)
            {
               //std::cout << "Recycling M: " << "\n";
               ComputeM = 0;
               HScalingFactor = 1.0;
            }
            else
            {
               //if ( not(ComputeJ) and any( HScalingFactor >= Qmin && HScalingFactor <= Qmax ) )
               //   std::cout << "Not recycling M due to any-all vote: " << HScalingFactor << any( HScalingFactor >= Qmin & HScalingFactor <= Qmax ) << "\n";
               ComputeM = 1;
            }
         }

         // Restrict the rate of change in dt
         HScalingFactor = fmax( HScalingFactor, 1.0 / this->adaption_limit);
         HScalingFactor = fmin( HScalingFactor,       this->adaption_limit);

         not_done = fabs(t - this->t_stop) > ValueType( this->t_round );

#if defined(VERBOSE) && (VERBOSE > 0)
         if (iter % VERBOSE == 0)
         {
            std::cout << "iter= " << iter;
            std::cout << " accept= " << Accepted;
            std::cout << " not_done= " << not_done;
            std::cout << " t= " << t;
            std::cout << " h= " << h;
            std::cout << " fact= " << HScalingFactor;
            std::cout << " T= " << y[getTempIndex(neq)] << "\n";
         }
#endif

         ValueType h0 = h;

         // Apply grow/shrink factor for next step.
         h *= HScalingFactor;

         // Limit based on the upper/lower bounds
         h = fmin(h, this->h_max);
         h = fmax(h, this->h_min);

         // Don't overshoot the final time ...
         h = select( not_done & ((t + h) > this->t_stop), this->t_stop - t, h );

         // Stretch the final step if we're really close
         // ... and we didn't just fail
         // ... and we're already re-computing M.
         ComputeM = any( h != h0 );
         if ( ComputeM )
         {
            MaskType Stretch_h = Accepted & ( fabs((t + h) - this->t_stop) < this->h_min );
            h = select( Stretch_h, this->t_stop - t, h );
         }

         ++iter;
         if ( this->max_iters && iter > this->max_iters )
         {
            ierr = ERR_TOO_MUCH_WORK;
            //printf("(iter > max_iters)\n");
            break;
         }
      }

      return ierr;
   }

   // SDIRK solver designed to replicate scalar logical behavior.
   template <class Functor, typename Counters>
   int solve_replicate ( ValueType *tcur, ValueType *hcur, Counters* counters, ValueType y[], Functor& func)
   {
      int ierr = ERR_SUCCESS;

      const int neq = this->neq;

      const int InterpolateNewton  = 1;	// Start at zero (0) or interpolate a starting guess (1)
      const int MaxNewtonIterations= 8;	// Max number of newton iterations
      const double NewtonThetaMin  = 0.005; // Minimum convergence rate for the Newton Iteration (0.001)
      const double NewtonThetaMax  = 0.999; // Maximum residual drop acceptable
      const double NewtonTolerance = 0.03;  // Convergence criteria
      const double Qmax = 1.2;		// Max h-adaption to recycle M
      const double Qmin = 1.;		// Min ""

      ValueType *rwk = m_rwk.getPointer();
      IndexType *iwk = m_iwk.getPointer();

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

      // Zero the counters ...
      nst = 0;
      nfe = 0;
      nlu = 0;
      nje = 0;
      nni = 0;
      iter = 0;

      // Set the work arrays ...
      ValueType *RESTRICT fy   = rwk;
      ValueType *RESTRICT del  = fy + (neq);
      ValueType *RESTRICT Jy   = del + (neq);
      ValueType *RESTRICT M    = Jy + (neq*neq);
      ValueType *RESTRICT z    = M + (neq*neq);
      ValueType *RESTRICT g    = z + (neq*this->numStages);
      ValueType *RESTRICT Mtmp = g + neq;
      ValueType *RESTRICT yerr = del;

      IndexType *RESTRICT ipiv = iwk;
      IndexType *RESTRICT ipivtmp = ipiv + (neq);

      MaskType ComputeJ (true);
      MaskType ComputeM (true);

      MaskType not_done = fabs(t - this->t_stop) > ValueType( this->t_round );

      while ( any(not_done) )
      {
         //std::cout << "Step: " << iter << nst << " " << ComputeM << " " << ComputeJ << not_done << t << "\n";

         const bool any_ComputeJ = any( ComputeJ );
         const bool all_ComputeJ = all( ComputeJ );
         const bool any_ComputeM = any( ComputeM );
         const bool all_ComputeM = all( ComputeM );

         // Construct the Iteration matrix ... if needed.
         if ( any(ComputeM) or any_ComputeJ )
         {
            // Compute the Jacobian matrix or recycle an old one.
            if ( any_ComputeJ )
            {
               //if (jac)
               //   jac (neq, t, y, Jy, user_data);
               //else
               {
                  // Compute the RHS ... the fd algorithm expects it.
                  func (neq, t, y, fy);

                  // Can skip the matrix copy if all need the Jac.
                  this->fdjac ( t, h, y, fy, (all_ComputeJ) ? Jy : Mtmp , func );

                  nfe = select( not_done & ComputeJ, nfe + (neq+1), nfe );
               }

               nje = select( not_done & ComputeJ, nje+1, nje );

               // Masked copy.
               if ( not( all_ComputeJ ) )
                  for (int i = 0; i < (neq*neq); ++i)
                     Jy[i] = select( ComputeJ, Mtmp[i], Jy[i] );
            }

            if ( any_ComputeM )
            {
               ValueType *_M = ( all_ComputeM ) ? M : Mtmp;
               IndexType *_ipiv = ( all_ComputeM ) ? ipiv : ipivtmp;

               // Compute M := 1/(gamma*h) - J
               const ValueType one_hgamma = 1.0 / (h * this->gamma);

               for (int j = 0; j < neq; ++j)
               {
                  ValueType *RESTRICT Mcol = &_M[(j*neq)];
                  ValueType *RESTRICT Jcol = &Jy[(j*neq)];
                  for (int i = 0; i < neq; ++i)
                     Mcol[(i)] = -Jcol[(i)];

                  Mcol[(j)] += one_hgamma;
               }

               // Factorization M
               this->ludec( neq, _M, _ipiv );

               nlu = select( not_done & ComputeM, nlu+1, nlu);

               if ( not( all_ComputeM ) )
               {
                  for (int i = 0; i < (neq*neq); ++i)
                     M[i] = select( ComputeM, Mtmp[i], M[i] );

                  for (int i = 0; i < neq; ++i)
                     ipiv[i] = select( ComputeM, ipivtmp[i], ipiv[i] );
               }
            }
         }

         typedef IndexType StatusType;
         const BaseIndexType Status_Active   = 0,
                             Status_Failed   = (1<<1),
                             Status_Inactive = (1<<2),
                             Status_Converged= (1<<3);

         MaskType isDone = not(not_done);
         StatusType Stage_Status = select( isDone, StatusType(Status_Inactive), StatusType(Status_Active) );

         MaskType  Accepted = false; // All are intialized inside of stage loop.
         ValueType HScalingFactor = select( isDone, ValueType(1.0), ValueType(0.8) );
         ValueType NewtonTheta = NewtonThetaMin;

         for (int s = 0; s < this->numStages && any(Stage_Status == Status_Active); s++)
         {
            //std::cout << "stage " << s << " " << Stage_Status << "\n";

            // Initial the RK stage vectors Z_i and G.
            ValueType *z_s = &z[(s*neq)];
            this->dzero (neq, z_s);
            this->dzero (neq, g);

            // Compute the function at this stage ...
            if (s)
            {
               for (int j = 0; j < s; ++j)
               {
                  // G = \Sum_j Theta_i,j*Z_j = h * \Sum_j A_i,j*F(Z_j)
                  ValueType *z_j = &z[(j*neq)];
                  this->daxpy1 (neq, Theta_(s,j), z_j, g);

                  // Z_i = \Sum_j Alpha(i,j)*Z_j
                  if (InterpolateNewton)
                     this->daxpy1 (neq, Alpha_(s,j), z_j, z_s);
               }
            }

            // Solve the non-linear problem with the Newton-Raphson method.
            StatusType Newton_Status = select( Stage_Status == Status_Active,
                                               StatusType( Status_Active),
                                               StatusType( Status_Inactive) );
            Accepted = false;
            HScalingFactor = select( Stage_Status == Status_Active, 0.8, HScalingFactor );
            NewtonTheta = select( Stage_Status == Status_Active, NewtonThetaMin, NewtonTheta );
            ValueType NewtonResidual;

            for (int NewtonIter = 0;
                     NewtonIter < MaxNewtonIterations && any(Newton_Status == Status_Active);
                     ++NewtonIter)
            {
               nni = select( (Newton_Status == Status_Active), nni+1, nni );

               // 1) Build the RHS of the root equation: F := G + (h*gamma)*f(y+Z_s) - Z_s
               for (int k = 0; k < neq; ++k)
                  del[(k)] = y[(k)] + z_s[(k)];

               func (neq, t, del, fy);
               nfe = select( (Newton_Status == Status_Active), nfe+1, nfe );

               const ValueType hgamma = h * this->gamma;
               for (int k = 0; k < neq; ++k)
                  del[(k)] = g[(k)] + hgamma * fy[(k)] - z_s[(k)];
                //del[(k)] = g[(k)] - z_s[(k)] + hgamma * fy[(k)];

               // 2) Solve the linear problem: M*delz = (1/hgamma)F ... so scale F first.
               this->dscal (neq, 1.0/hgamma, del);
               this->lusol (neq, M, ipiv, del);

               // 3) Check the convergence and convergence rate
               // 3.1) Compute the norm of the correction.
               const ValueType dnorm = this->wnorm ( del, y);

               //std::cout << "miter: " << iter <<" " << s << " "<< NewtonIter << dnorm << NewtonResidual << Newton_Status << "\n";

               // 3.2) If not the first iteration, estimate the rate.
               if (NewtonIter != 0)
               {
                  NewtonTheta = dnorm / NewtonResidual; // How much as the residual changed from last iteration.

                  MaskType Diverging  = (NewtonTheta >= NewtonThetaMax);
                           Diverging &= (Newton_Status == Status_Active);

                  ValueType ConvergenceRate = NewtonTheta / (1.0 - NewtonTheta); // Possible FLTEXP here.

                  Newton_Status = select( Diverging, StatusType(Status_Failed), Newton_Status );

                  //std::cout << " " << ConvergenceRate << Diverging << NewtonTheta << Newton_Status << "\n";

                  // If one is diverging, may just give up early all around.
                  // Use the status flag to change or hold h.
                  if ( all( Newton_Status != Status_Active ) )
                  {
                     //std::cout << "all inactive @ (NewtonTheta >= NewtonThetaMax)\n";
                     break;
                  }

                  // So, not everyone is diverging. Test if any are converging too slowly.
                  // We'll shrink h and restart, too.
                  // Predict the error after Max iterations with the current rate:
                  // ... res * Theta^(ItersLeft/(1-Theta))
                  ValueType PredictedError = dnorm * pow( NewtonTheta,
                                       (MaxNewtonIterations-NewtonIter)/(1.0-NewtonTheta));

                  MaskType PredictedErrorTooLarge  = ( PredictedError > NewtonTolerance );
                           PredictedErrorTooLarge &= ( Newton_Status == Status_Active );

                  if ( any( PredictedErrorTooLarge ) )
                  {
                     //std::cout << "PredictedErrorTooLarge " << PredictedErrorTooLarge << "\n";

                     // Doubtful we'll converge, shrink h and try again.
                     ValueType QNewton = fmin(10.0, PredictedError / NewtonTolerance);
                     ValueType PredictedHScalingFactor =
                                  0.9 * pow( QNewton, -1.0 / (1.0 + MaxNewtonIterations - NewtonIter));

                     HScalingFactor = select( PredictedErrorTooLarge, PredictedHScalingFactor, HScalingFactor );
                     Newton_Status = select( PredictedErrorTooLarge, StatusType(Status_Failed), Newton_Status );

                     // Don't let a finished lane hold up an exit.
                     if ( all( Newton_Status != Status_Active ) )
                     {
                        //std::cout << "all inactive @ PredictedErrorTooLarge\n";
                        break;
                     }
                  }

                  // So, at least someone is making progress. Test for new convergence.
                  MaskType isConverged =   ( ConvergenceRate * dnorm < NewtonTolerance )
                                         & MaskType( Newton_Status == Status_Active );

                  Accepted |= isConverged;
                  //std::cout << "isConverged: " << isConverged << Newton_Status << "\n";
               }

               // Save the residual norm for the next iteration unless I'm already converged.
               MaskType UpdateSolution =  ( Newton_Status == Status_Active );
               NewtonResidual = select( UpdateSolution, dnorm, NewtonResidual );

               // 4) Update the solution if newly accepted: Z_s <- Z_s + delta
               ValueType OneOrZero = select( UpdateSolution, ValueType(1), ValueType(0) );
               this->daxpy (neq, OneOrZero, del, z_s);

               Newton_Status = select( Accepted, StatusType( Status_Inactive ), Newton_Status );
            }

            Stage_Status = select( Newton_Status == Status_Failed,
                                   StatusType( Status_Failed ), Stage_Status );
            //std::cout << "Accepted: " << Accepted << " " << Stage_Status << Newton_Status << "\n";

            //ComputeJ = select( Accepted, ComputeJ, MaskType(false) ); // Jacobian is still valid if not Accepted.
            //ComputeM = select( Accepted, ComputeM, MaskType(true ) ); // M is invalid since h will change (perhaps) drastically.
            ComputeJ =     false; // Jacobian is still valid if not Accepted.
            ComputeM =  Accepted; // M is invalid since h will change (perhaps) drastically.

         } // ... stages

         if ( any( Accepted ) )
         {
            // Compute the error estimation of the trial solution.
            this->dzero (neq, yerr);
            for (int j = 0; j < this->numStages; ++j)
            {
               ValueType *z_j = &z[(j*neq)];
               if (this->E[j] != 0.0)
                  this->daxpy1 (neq, this->E[j], z_j, yerr);
            }

            ValueType herr = fmax(1.0e-20, this->wnorm ( yerr, y));

            // Is the error acceptable?
            MaskType AcceptedSolution = Accepted && (((herr <= 1.0) || (h <= this->h_min)) && not_done);
            if ( any( AcceptedSolution ) )
            {
               // If stiffly-accurate, Z_s with s := numStages, is the solution.
               // Else, sum the stage solutions: y_n+1 <- y_n + \Sum_j D_j * Z_j
               for (int j = 0; j < this->numStages; ++j)
               {
                  ValueType *z_j = &z[(j*neq)];
                  if (this->D[j] != 0.0)
                  {
                     // Only update solutions that were accepted.
                     ValueType ZeroOrDj = select( AcceptedSolution, ValueType(this->D[j]), 0.0 );
                     this->daxpy ( neq, ZeroOrDj, z_j, y );
                  }
               }

               t = select( AcceptedSolution, t + h, t );
               nst = select( AcceptedSolution, nst+1, nst );
            }

            not_done = fabs(t - this->t_stop) > ValueType( this->t_round );

            HScalingFactor = select( Accepted, 0.9 * pow( 1.0 / herr, (1.0 / this->ELO)),
                                                       /*failed*/ HScalingFactor );
            HScalingFactor = select( not_done, HScalingFactor, /*hold fixed*/ ValueType(1) );

            // Reuse the Jacobian if the Newton Solver is converging fast enough.
            // ... no need to recompute if I've failed since nothing has changed.
            ComputeJ = Accepted & ( NewtonTheta > NewtonThetaMin ) & not_done;

            // Don't refine if it's not a big step and we could reuse the M matrix.
            MaskType recycle_M = not ( ComputeJ )
                                 and ( HScalingFactor >= Qmin && HScalingFactor <= Qmax )
                                 and ( Accepted );

            //std::cout << "Recycle M: " << recycle_M << "\n";
            ComputeM       = not( recycle_M );
            HScalingFactor = select( recycle_M, ValueType(1), HScalingFactor );
         }

         // Restrict the rate of change in dt
         HScalingFactor = fmax( HScalingFactor, 1.0 / this->adaption_limit);
         HScalingFactor = fmin( HScalingFactor,       this->adaption_limit);

         not_done = fabs(t - this->t_stop) > ValueType( this->t_round );

#if defined(VERBOSE) && (VERBOSE > 0)
         if (iter % VERBOSE == 0)
         {
            std::cout << "iter= " << iter;
            std::cout << " accept= " << Accepted;
            std::cout << " not_done= " << not_done;
            std::cout << " t= " << t;
            std::cout << " h= " << h;
            std::cout << " fact= " << HScalingFactor;
            std::cout << " T= " << y[getTempIndex(neq)] << "\n";
         }
#endif

         ValueType h0 = h;

         // Apply grow/shrink factor for next step.
         h *= HScalingFactor;

         // Limit based on the upper/lower bounds
         h = fmin(h, this->h_max);
         h = fmax(h, this->h_min);

         // Stretch the final step if we're really close
         // ... and we didn't just fail
         MaskType Stretch_h = Accepted & ( fabs((t + h) - this->t_stop) < this->h_min );
         h = select( Stretch_h, this->t_stop - t, h );

         // Don't overshoot the final time ...
         h = select( not_done & ((t + h) > this->t_stop), this->t_stop - t, h );

         ComputeM = ( h != h0 );

         ++iter;
         if ( this->max_iters && iter > this->max_iters )
         {
            ierr = ERR_TOO_MUCH_WORK;
            //printf("(iter > max_iters)\n");
            break;
         }
      }

      return ierr;
   }

   #undef nst
   #undef nfe
   #undef nje
   #undef nlu
   #undef nni
   #undef iter
   #undef h
   #undef t

   // SDIRK solver designed to replicate scalar logical behavior.
   template <class Functor, typename Counters>
   int solve ( ValueType *tcur, ValueType *hcur, Counters* counters, ValueType y[], Functor& func)
   {
      if ( _replicate )
         return this->solve_replicate ( tcur, hcur, counters, y, func );
      else
         return this->solve_anyall ( tcur, hcur, counters, y, func );
   }

};
   #undef __matrix_index
   #undef A_
   #undef Theta_
   #undef Alpha_

template <typename Functor, typename RHSptr>
void simd_rk_driver ( const int numProblems, const double *u_in, const double t_stop, const Functor& func, const RHSptr rhs_func, const ckdata_t *RESTRICT ck )
{
   const int kk = ck->n_species;
   const int neq = kk+1;
   const double p = func.getPressure();

   typedef typename SIMD_TypeSelector<double,DefaultVectorLength>::value_type SimdType;
   typedef typename SIMD_MaskSelector<SimdType>::mask_type MaskType;
   const int VectorLength = SIMD_Length<SimdType>::length;

   //printf("Instruction Set= %d %s %d %s\n", INSTRSET, typeid(SimdType).name(), VectorLength, typeid( MaskType).name());
   printf("SimdType= %s %d %s\n", typeid(SimdType).name(), VectorLength, typeid( MaskType).name());

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

      if (i % VectorLength == 0)
         printf("%d: final %d %d %d %e %e %f %f\n", i, _nst, _nit, _nfe, u[ getTempIndex(neq) ], T0, (u[getTempIndex(neq)]-T0)/T0, 1000*(t_end-t_begin));
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

   nst = 0, nit = 0, nfe = 0;

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

         std::cout << i0 << ": final " << counters.nst
                         << " "  << counters.nit
                         << " "  << v_u[getTempIndex(neq)]
                         << " "  << T0
                         << " "  << (v_u[getTempIndex(neq)]-T0)/T0 << "\n";
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
      double errmax = 0;
      int ierrmax = -1;
      for (int i = 0; i < numProblems; ++i)
      {
         const double *v_out = vector_out.getPointer() + neq*i;
         const double *s_out = scalar_out.getPointer() + neq*i;
         double diff = std::abs( s_out[ getTempIndex(neq) ]
                               - v_out[ getTempIndex(neq) ] );
         err2 += sqr( diff );
         ref2 += sqr( s_out[ getTempIndex(neq) ] );

         if ( diff > errmax ) {
            errmax = diff;
            ierrmax = i;
         }
      }

      printf("err2= %e %e %e %d %e %d\n", err2, ref2, std::sqrt(err2)/std::sqrt(ref2), numProblems % VectorLength, errmax, ierrmax);
   }
   {
      double err2 = 0, ref2 = 0;
      double errmax = 0;
      int ierrmax = -1;
      for (int k = 0; k < kk; ++k)
         for (int i = 0; i < numProblems; ++i)
         {
            const double *v_out = vector_out.getPointer() + neq*i;
            const double *s_out = scalar_out.getPointer() + neq*i;
            double diff = std::abs( s_out[ getFirstSpeciesIndex(neq)+k ]
                                  - v_out[ getFirstSpeciesIndex(neq)+k ] );
            err2 += sqr( diff );
            ref2 += sqr( s_out[ getFirstSpeciesIndex(neq)+k ] );

            if ( diff > errmax ) {
               errmax = diff;
               ierrmax = i;
            }
         }

      printf("err2= %e %e %e %e %d\n", err2, ref2, std::sqrt(err2)/std::sqrt(ref2), errmax, ierrmax);
   }

   return;
}

template <typename Functor, typename RHSptr>
void simd_ros_driver ( const int numProblems, const double *u_in, const double t_stop, const Functor& func, const RHSptr rhs_func, const ckdata_t *RESTRICT ck )
{
   const int kk = ck->n_species;
   const int neq = kk+1;
   const double p = func.getPressure();

   typedef typename SIMD_TypeSelector<double,DefaultVectorLength>::value_type SimdType;
   typedef typename SIMD_MaskSelector<SimdType>::mask_type MaskType;
   const int VectorLength = SIMD_Length<SimdType>::length;

   printf("SimdType= %s %d %s\n", typeid(SimdType).name(), VectorLength, typeid( MaskType).name());

   VectorType<double,Alignment> scalar_out( neq * numProblems );
   VectorType<double,Alignment> vector_out( neq * numProblems );

   ros_t ros;

   ros_create (&ros, neq, Ros4);

   ros.max_iters = 1000;
   ros.min_iters = 1;

   int lenrwk = ros_lenrwk (&ros);
   int leniwk = ros_leniwk (&ros);
   VectorType<double,Alignment> rwk( VectorLength * lenrwk );
   VectorType<int,Alignment> iwk( VectorLength * leniwk );
   VectorType<double,Alignment> u( VectorLength * neq );

   int nst = 0, nit = 0, nfe = 0, nje = 0;

   auto scalar_solver = [&](const int i, VectorType<double,Alignment>& out, ros_counters_t *counters)
   {
      for (int k = 0; k < neq; ++k)
         u[k] = u_in[ i*neq + k ];

      const double T0 = u_in[ i*neq + getTempIndex(neq) ];

      double t = 0, h = 0;

      ros_init (&ros, t, t_stop);

      double t_begin = WallClock();

      int ierr = ros_solve (&ros, &t, &h, counters, u.getPointer(), iwk.getPointer(), rwk.getPointer(), rhs_func, /*jac_func*/NULL, (void*)&func);
      if (ierr != ROS_SUCCESS)
      {
         fprintf(stderr,"%d: ros_solve error %d %d %d\n", i, ierr, counters->niters, ros.max_iters);
         exit(-1);
      }

      double t_end = WallClock();

      const int _nit = counters->niters;
      const int _nst = counters->nst;
      const int _nfe = counters->nfe;
      const int _nje = counters->nje;

      nit += _nit;
      nst += _nst;
      nfe += _nfe;
      nje += _nje;

      for (int k = 0; k < neq; ++k)
         out[i*neq + k] = u[k];

      if (i % VectorLength == 0)
         printf("%d: final %d %d %d %e %e %f %f\n", i, _nst, _nit, _nfe, u[ getTempIndex(neq) ], T0, (u[ getTempIndex(neq) ]-T0)/T0, 1000*(t_end-t_begin));
   };

   double time_scalar = WallClock();

   for (int i = 0; i < numProblems; ++i)
   {
      ros_counters_t counters;
      scalar_solver(i, scalar_out, &counters);
      //if (i % 10 == 0)
      //   printf("%d: %d %d %f\n", i, counters.nsteps, counters.niters, scalar_out[i*neq+getTempIndex(neq)] );
   }

   ros_destroy(&ros);

   time_scalar = WallClock() - time_scalar;

   printf("Scalar: %f %d %d %d %d\n", 1000*time_scalar, nst, nit, nfe, nje);

   simd_cklib_functor<SimdType> simd_func( ck, p );
   typedef SimdRosSolverType<SimdType> SimdSolverType;
   SimdSolverType simd_solver( neq );

   nst = nit = nfe = nje = 0;

   double time_vector = WallClock();

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

         nit += counters.nit * VectorLength;

         nst += sum( counters.nst );
         nfe += sum( counters.nfe );
         nje += sum( counters.nje );
         //for (int i = 0; i < VectorLength; ++i)
         //{
         //   nst += counters.nst[i];
         //   nfe += counters.nfe[i];
         //   nje += counters.nje[i];
         //}

         std::cout << i0 << ": final " << counters.nst
                         << " " << counters.nit
                         << " " << counters.nfe
                         << " " << v_u[getTempIndex(neq)]
                         << " " << T0
                         << " " << ((v_u[getTempIndex(neq)]-T0)/T0) << "\n";
      }
      else
      {
         for (int i = i0; i < numProblems; ++i)
         {
            ros_counters_t counters;

            scalar_solver( i, vector_out, &counters );

            nit += counters.niters;
            nst += counters.nst;
            nfe += counters.nfe;
            nje += counters.nje;
         }
      }
   }

   time_vector = WallClock() - time_vector;

   printf("SIMD: %f %d %d %d %d\n", 1000*time_vector, nst, nit, nfe, nje);

   printf("SIMD timer: %f %f %.1f\n", 1000.*time_vector, 1000.*time_scalar, time_scalar/time_vector);

   {
      double err2 = 0, ref2 = 0;
      double errmax = 0;
      int ierrmax = -1;
      for (int i = 0; i < numProblems; ++i)
      {
         const double *v_out = vector_out.getPointer() + neq*i;
         const double *s_out = scalar_out.getPointer() + neq*i;
         double diff = std::abs( s_out[ getTempIndex(neq) ]
                               - v_out[ getTempIndex(neq) ] );
         err2 += sqr( diff );
         ref2 += sqr( s_out[ getTempIndex(neq) ] );

         if ( diff > errmax ) {
            errmax = diff;
            ierrmax = i;
         }
      }

      printf("err2= %e %e %e %d %e %d\n", err2, ref2, std::sqrt(err2)/std::sqrt(ref2), numProblems % VectorLength, errmax, ierrmax);
   }
   {
      double err2 = 0, ref2 = 0;
      double errmax = 0;
      int ierrmax = -1;
      for (int k = 0; k < kk; ++k)
         for (int i = 0; i < numProblems; ++i)
         {
            const double *v_out = vector_out.getPointer() + neq*i;
            const double *s_out = scalar_out.getPointer() + neq*i;
            double diff = std::abs( s_out[ getFirstSpeciesIndex(neq)+k ]
                                  - v_out[ getFirstSpeciesIndex(neq)+k ] );
            err2 += sqr( diff );
            ref2 += sqr( s_out[ getFirstSpeciesIndex(neq)+k ] );

            if ( diff > errmax ) {
               errmax = diff;
               ierrmax = i;
            }
         }

      printf("err2= %e %e %e %e %d\n", err2, ref2, std::sqrt(err2)/std::sqrt(ref2), errmax, ierrmax);
   }

   return;
}

template <typename Functor, typename RHSptr>
void simd_sdirk_driver ( const int numProblems, const double *u_in, const double t_stop, const Functor& func, const RHSptr rhs_func, const ckdata_t *RESTRICT ck )
{
   const int kk = ck->n_species;
   const int neq = kk+1;
   const double p = func.getPressure();

   typedef typename SIMD_TypeSelector<double,DefaultVectorLength>::value_type SimdType;
   typedef typename SIMD_MaskSelector<SimdType>::mask_type MaskType;
   const int VectorLength = SIMD_Length<SimdType>::length;

   printf("SimdType = %s %d %s\n", typeid(SimdType).name(), VectorLength, typeid( MaskType).name());

   VectorType<double,Alignment> scalar_out( neq * numProblems );
   VectorType<double,Alignment> vector_out( neq * numProblems );

   sdirk_t sdirk;

   sdirk_create (&sdirk, neq, S4a);

   sdirk.max_iters = 1000;
   sdirk.min_iters = 1;

   int lenrwk = sdirk_lenrwk (&sdirk);
   int leniwk = sdirk_leniwk (&sdirk);
   VectorType<double,Alignment> rwk( VectorLength * lenrwk );
   VectorType<int,Alignment> iwk( VectorLength * leniwk );
   VectorType<double,Alignment> u( VectorLength * neq );

   int nst = 0, nit = 0, nfe = 0, nje = 0, nlu = 0, nni = 0;

   auto scalar_solver = [&](const int i, VectorType<double,Alignment>& out, sdirk_counters_t *counters)
   {
      for (int k = 0; k < neq; ++k)
         u[k] = u_in[ i*neq + k ];

      const double T0 = u_in[ i*neq + getTempIndex(neq) ];

      double t = 0, h = 0;

      sdirk_init (&sdirk, t, t_stop);

      double t_begin = WallClock();

      int ierr = sdirk_solve (&sdirk, &t, &h, counters, u.getPointer(), iwk.getPointer(), rwk.getPointer(), rhs_func, /*jac_func*/NULL, (void*)&func);
      if (ierr != ROS_SUCCESS)
      {
         fprintf(stderr,"%d: sdirk_solve error %d %d %d\n", i, ierr, counters->niters, sdirk.max_iters);
         exit(-1);
      }

      double t_end = WallClock();

      const int _nit = counters->niters;
      const int _nst = counters->nst;
      const int _nfe = counters->nfe;
      const int _nje = counters->nje;
      const int _nlu = counters->nlu;
      const int _nni = counters->nni;

      nit += _nit;
      nst += _nst;
      nfe += _nfe;
      nje += _nje;
      nlu += _nlu;
      nni += _nni;

      for (int k = 0; k < neq; ++k)
         out[i*neq + k] = u[k];

      //if (i % VectorLength == 0)
         printf("%d: final %d %d %d %d %d %e %e %f %f\n", i, _nst, _nit, _nfe, _nlu, _nni, u[ getTempIndex(neq) ], T0, (u[ getTempIndex(neq) ]-T0)/T0, 1000*(t_end-t_begin));
   };

   double time_scalar = WallClock();

   for (int i = 0; i < numProblems; ++i)
   {
      sdirk_counters_t counters;
      scalar_solver(i, scalar_out, &counters);
      //if (i % 10 == 0)
      //   printf("%d: %d %d %f\n", i, counters.nsteps, counters.niters, scalar_out[i*neq+getTempIndex(neq)] );
   }

   sdirk_destroy(&sdirk);

   time_scalar = WallClock() - time_scalar;

   printf("Scalar: %f %d %d %d %d %d %d\n", 1000*time_scalar, nst, nit, nfe, nje, nlu, nni);

   simd_cklib_functor<SimdType> simd_func( ck, p );
   typedef SimdSdirkSolverType<SimdType,1> SimdSolverType;
   SimdSolverType simd_solver( neq );
   //simd_solver.max_iters=100;

   nst = nit = nfe = nlu = nje = nni = 0;

   double time_vector = WallClock();

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

         nit += counters.nit * VectorLength;

         nst += sum( counters.nst );
         nfe += sum( counters.nfe );
         nje += sum( counters.nje );
         nlu += sum( counters.nlu );
         nni += sum( counters.nni );
         //for (int i = 0; i < VectorLength; ++i)
         //{
         //   nst += counters.nst[i];
         //   nfe += counters.nfe[i];
         //   nje += counters.nje[i];
         //   nlu += counters.nlu[i];
         //   nni += counters.nni[i];
         //}

         std::cout << i0 << ": final " << counters.nst
                         << " " << counters.nit
                         << " " << counters.nfe
                         << " " << counters.nlu
                         << " " << counters.nni
                         << " " << v_u[getTempIndex(neq)]
                         << " " << T0
                         << " " << ((v_u[getTempIndex(neq)]-T0)/T0) << "\n";
      }
      else
      {
         for (int i = i0; i < numProblems; ++i)
         {
            sdirk_counters_t counters;

            scalar_solver( i, vector_out, &counters );

            nit += counters.niters;
            nst += counters.nst;
            nfe += counters.nfe;
            nje += counters.nje;
            nlu += counters.nlu;
            nni += counters.nni;
         }
      }
   }

   time_vector = WallClock() - time_vector;

   printf("SIMD: %f %d %d %d %d %d %d %s\n", 1000*time_vector, nst, nit, nfe, nje, nlu, nni, SimdSolverType::_replicate ? "Replicate-Serial" : "Any-All" );

   printf("SIMD speed-up: %f %f %.1f\n", 1000.*time_vector, 1000.*time_scalar, time_scalar/time_vector);

   {
      double err2 = 0, ref2 = 0;
      double errmax = 0;
      int ierrmax = -1;
      for (int i = 0; i < numProblems; ++i)
      {
         const double *v_out = vector_out.getPointer() + neq*i;
         const double *s_out = scalar_out.getPointer() + neq*i;
         double diff = std::abs( s_out[ getTempIndex(neq) ]
                               - v_out[ getTempIndex(neq) ] );
         err2 += sqr( diff );
         ref2 += sqr( s_out[ getTempIndex(neq) ] );

         if ( diff > errmax ) {
            errmax = diff;
            ierrmax = i;
         }
      }

      printf("err2= %e %e %e %d %e %d\n", err2, ref2, std::sqrt(err2)/std::sqrt(ref2), numProblems % VectorLength, errmax, ierrmax);
   }
   {
      double err2 = 0, ref2 = 0;
      double errmax = 0;
      int ierrmax = -1;
      for (int k = 0; k < kk; ++k)
         for (int i = 0; i < numProblems; ++i)
         {
            const double *v_out = vector_out.getPointer() + neq*i;
            const double *s_out = scalar_out.getPointer() + neq*i;
            double diff = std::abs( s_out[ getFirstSpeciesIndex(neq)+k ]
                                  - v_out[ getFirstSpeciesIndex(neq)+k ] );
            err2 += sqr( diff );
            ref2 += sqr( s_out[ getFirstSpeciesIndex(neq)+k ] );

            if ( diff > errmax ) {
               errmax = diff;
               ierrmax = i;
            }
         }

      printf("err2= %e %e %e %e %d\n", err2, ref2, std::sqrt(err2)/std::sqrt(ref2), errmax, ierrmax);
   }

   return;
}

} // namespace

#endif // ifndef
