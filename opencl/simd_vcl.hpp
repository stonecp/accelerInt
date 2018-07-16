#ifndef __simd_vcl_h
#define __simd_vcl_h

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

// Master VCL header
#ifndef  MAX_VECTOR_SIZE
# define MAX_VECTOR_SIZE (512)
#endif
#define VCL_NAMESPACE VCL
#include "vcl/vectorclass.h"
#include "vcl/vectormath_exp.h"

using namespace VCL;

template <typename T, int VL>
struct SIMD_TypeSelector;

template<> struct SIMD_TypeSelector<double,2> { typedef VCL::Vec2d value_type; };
template<> struct SIMD_TypeSelector<double,4> { typedef VCL::Vec4d value_type; };
template<> struct SIMD_TypeSelector<double,8> { typedef VCL::Vec8d value_type; };
template<> struct SIMD_TypeSelector<int64_t,2> { typedef VCL::Vec2q value_type; };
template<> struct SIMD_TypeSelector<int64_t,4> { typedef VCL::Vec4q value_type; };
template<> struct SIMD_TypeSelector<int64_t,8> { typedef VCL::Vec8q value_type; };

template <typename V>
struct SIMD_MaskSelector;

template<> struct SIMD_MaskSelector<VCL::Vec2d> { typedef VCL::Vec2db mask_type; };
template<> struct SIMD_MaskSelector<VCL::Vec4d> { typedef VCL::Vec4db mask_type; };
template<> struct SIMD_MaskSelector<VCL::Vec8d> { typedef VCL::Vec8db mask_type; };
template<> struct SIMD_MaskSelector<VCL::Vec2q> { typedef VCL::Vec2qb mask_type; };
template<> struct SIMD_MaskSelector<VCL::Vec4q> { typedef VCL::Vec4qb mask_type; };
template<> struct SIMD_MaskSelector<VCL::Vec8q> { typedef VCL::Vec8qb mask_type; };

template <typename V>
struct SIMD_Length;

template<> struct SIMD_Length<VCL::Vec2d> { enum { length = 2 }; };
template<> struct SIMD_Length<VCL::Vec4d> { enum { length = 4 }; };
template<> struct SIMD_Length<VCL::Vec8d> { enum { length = 8 }; };

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

template <typename ValueType>
inline ValueType fmax( const ValueType& a, const ValueType& b ) { return select( a > b, a, b ); }
template <typename ValueType>
inline ValueType fmax( const ValueType& a, const double b ) { return select( a > ValueType(b), a, b); }
template <typename ValueType>
inline ValueType fmax( const double a, const ValueType& b ) { return fmax(b,a); }

template <typename ValueType>
inline ValueType fmin( const ValueType& a, const ValueType& b ) { return select( a < b, a, b ); }
template <typename ValueType>
inline ValueType fmin( const ValueType& a, const double b ) { return select( a < ValueType(b), a, b ); }
template <typename ValueType>
inline ValueType fmin( const double a, const ValueType& b ) { return fmin(b,a); }

} // namespace

#endif // ifndef
