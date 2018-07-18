#ifndef __simd_vcl_h
#define __simd_vcl_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <cmath>
#include <string>
#include <typeinfo>

#include "simd_dtypes.hpp"

namespace SIMD
{

// Master VCL header
#ifndef  MAX_VECTOR_SIZE
# define MAX_VECTOR_SIZE (512)
#endif

#define VCL_NAMESPACE VCL
#include "vcl/vectorclass.h"
#include "vcl/vectormath_exp.h"

using namespace VCL;

template<> struct SIMD_TypeSelector<double,2> { typedef VCL::Vec2d value_type; };
template<> struct SIMD_TypeSelector<double,4> { typedef VCL::Vec4d value_type; };
template<> struct SIMD_TypeSelector<double,8> { typedef VCL::Vec8d value_type; };
template<> struct SIMD_TypeSelector<int64_t,2> { typedef VCL::Vec2q value_type; };
template<> struct SIMD_TypeSelector<int64_t,4> { typedef VCL::Vec4q value_type; };
template<> struct SIMD_TypeSelector<int64_t,8> { typedef VCL::Vec8q value_type; };

template<> struct SIMD_MaskSelector<VCL::Vec2d> { typedef VCL::Vec2db mask_type; };
template<> struct SIMD_MaskSelector<VCL::Vec4d> { typedef VCL::Vec4db mask_type; };
template<> struct SIMD_MaskSelector<VCL::Vec8d> { typedef VCL::Vec8db mask_type; };
template<> struct SIMD_MaskSelector<VCL::Vec2q> { typedef VCL::Vec2qb mask_type; };
template<> struct SIMD_MaskSelector<VCL::Vec4q> { typedef VCL::Vec4qb mask_type; };
template<> struct SIMD_MaskSelector<VCL::Vec8q> { typedef VCL::Vec8qb mask_type; };

template<> struct SIMD_Length<VCL::Vec2d> { enum { length = 2 }; };
template<> struct SIMD_Length<VCL::Vec4d> { enum { length = 4 }; };
template<> struct SIMD_Length<VCL::Vec8d> { enum { length = 8 }; };

template <> struct SIMD_ScalarType< VCL::Vec2d > { typedef double  scalar_type; };
template <> struct SIMD_ScalarType< VCL::Vec4d > { typedef double  scalar_type; };
template <> struct SIMD_ScalarType< VCL::Vec8d > { typedef double  scalar_type; };
template <> struct SIMD_ScalarType< VCL::Vec2q > { typedef int64_t scalar_type; };
template <> struct SIMD_ScalarType< VCL::Vec4q > { typedef int64_t scalar_type; };
template <> struct SIMD_ScalarType< VCL::Vec8q > { typedef int64_t scalar_type; };

template <typename MaskType>
inline bool all( const MaskType& mask ) { return horizontal_and( mask ); }

template <typename MaskType>
inline bool any( const MaskType& mask ) { return horizontal_or( mask ); }

template <typename SimdType>
inline
   typename std::enable_if<
                  SIMD_isScalar<SimdType>::value == false,
                  std::string
                          >::type
toString (const SimdType& x)
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

//template <typename ValueType>
//inline ValueType fmax( const ValueType& a, const ValueType& b ) { return select( a > b, a, b ); }
//template <typename ValueType>
//inline ValueType fmax( const ValueType& a, const double b ) { return select( a > ValueType(b), a, b); }
//template <typename ValueType>
//inline ValueType fmax( const double a, const ValueType& b ) { return fmax(b,a); }
template <typename SimdType>
inline
typename std::enable_if<
                 (SIMD_isScalar<SimdType>::value == 0),
                  SimdType
                          >::type
   fmax( const SimdType& a, const SimdType& b ) { return select( a > b, a, b ); }

template <typename SimdType>
inline
typename std::enable_if<
                 (SIMD_isScalar<SimdType>::value == 0),
                  SimdType
                          >::type
   fmax( const SimdType& a, const double b ) { return select( a > SimdType(b), a, b); }
template <typename SimdType>
inline
typename std::enable_if<
                 (SIMD_isScalar<SimdType>::value == 0),
                  SimdType
                          >::type
   fmax( const double a, const SimdType& b ) { return SIMD::fmax<SimdType>(b,a); }

template <typename ValueType>
inline ValueType fmin( const ValueType& a, const ValueType& b ) { return select( a < b, a, b ); }
template <typename ValueType>
inline ValueType fmin( const ValueType& a, const double b ) { return select( a < ValueType(b), a, b ); }
template <typename ValueType>
inline ValueType fmin( const double a, const ValueType& b ) { return fmin(b,a); }

template <typename SimdType>
inline
typename std::enable_if<
                 (SIMD_isScalar<SimdType>::value == 0),
                  SimdType
                          >::type
   fabs( const SimdType& v ) { return VCL::abs(v); }

template <typename SimdType>
inline
   typename std::enable_if<
                  SIMD_isScalar<SimdType>::value == 0,
                  typename SIMD_ScalarType<SimdType>::scalar_type
                          >::type
sum (const SimdType& v) { return horizontal_add( v ); }

} // namespace

#endif // ifndef
