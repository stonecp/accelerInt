#ifndef __simd_dtypes_h
#define __simd_dtypes_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <cmath>
#include <string>
#include <typeinfo>

namespace SIMD
{

#if !defined(__ALIGNMENT)
#warning 'Fixing alignment = 64'
const size_t Alignment = 64;
#else
const size_t Alignment = __ALIGNMENT;
#endif

#ifndef __DEFAULT_VECTOR_LENGTH
#define __DEFAULT_VECTOR_LENGTH (4)
#endif

const int DefaultVectorLength = __DEFAULT_VECTOR_LENGTH;

template <typename T, int VL>
struct SIMD_TypeSelector;

// Scalar version
template<> struct SIMD_TypeSelector<double,1>  { typedef double  value_type; };
template<> struct SIMD_TypeSelector<int,1>     { typedef int     value_type; };
template<> struct SIMD_TypeSelector<int64_t,1> { typedef int64_t value_type; };

template <typename V>
struct SIMD_MaskSelector;

template<> struct SIMD_MaskSelector<double>  { typedef bool mask_type; };
template<> struct SIMD_MaskSelector<int>     { typedef bool mask_type; };
template<> struct SIMD_MaskSelector<int64_t> { typedef bool mask_type; };

template <typename V>
struct SIMD_isScalar { enum { value = 0 }; };

template<> struct SIMD_isScalar <double>  { enum { value = 1 }; };
template<> struct SIMD_isScalar <int>     { enum { value = 1 }; };
template<> struct SIMD_isScalar <int64_t> { enum { value = 1 }; };

template <typename V>
struct SIMD_ScalarType {};

template <> struct SIMD_ScalarType<double > { typedef double  scalar_type; };
template <> struct SIMD_ScalarType<int    > { typedef int     scalar_type; };
template <> struct SIMD_ScalarType<int64_t> { typedef int64_t scalar_type; };
template <> struct SIMD_ScalarType<bool   > { typedef bool    scalar_type; };

template <typename V>
struct SIMD_Length;

template<> struct SIMD_Length<double>  { enum { length = 1 }; };
template<> struct SIMD_Length<int>     { enum { length = 1 }; };
template<> struct SIMD_Length<int64_t> { enum { length = 1 }; };

template <typename MaskType>
inline bool all( const MaskType& mask );

template <>
inline bool all( const bool& mask ) { return mask; }

template <typename MaskType>
inline bool any( const MaskType& mask );

template <>
inline bool any( const bool& mask ) { return mask; }

//template <typename SimdType>
//std::string toString (const SimdType& x);

template <typename ScalarType>
inline
   typename std::enable_if<
                  SIMD_isScalar<ScalarType>::value,
                  std::string
                          >::type
toString (const ScalarType& x)
{
   std::ostringstream oss;
   oss << "[" << x << "]";

   return std::string( oss.str() );
}

template <typename ScalarType>
inline
   typename std::enable_if<
                  SIMD_isScalar<ScalarType>::value,
                  typename SIMD_ScalarType<ScalarType>::scalar_type
                          >::type
sum (const ScalarType& v) { return v; }

template <typename ScalarType>
inline
typename std::enable_if<
                  SIMD_isScalar<ScalarType>::value,
                  ScalarType
                          >::type
select ( const bool& mask, const ScalarType& ifTrue, const ScalarType& ifFalse )
{
   return mask ? ifTrue : ifFalse;
}

//template <typename ScalarType>
//inline
//typename std::enable_if<
//                  SIMD_isScalar<ScalarType>::value,
//                  ScalarType
//                          >::type
//   fmax( const ScalarType& a, const ScalarType& b ) { return std::fmax(a,b); }
//
//template <typename ScalarType>
//inline
//typename std::enable_if<
//                  SIMD_isScalar<ScalarType>::value,
//                  ScalarType
//                          >::type
//   fabs( const ScalarType& a ) { return std::fabs(a); }

} // namespace

#if defined(__DEFAULT_VECTOR_LENGTH) && ( __DEFAULT_VECTOR_LENGTH > 1 )
# include "simd_vcl.hpp"
#endif

#endif // ifndef
