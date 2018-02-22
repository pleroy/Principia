
#pragma once

#include "quantities/wide.hpp"

namespace principia {
namespace quantities {
namespace internal_wide {

template<typename Q>
Quantity128<Q>::Quantity128(Q const x) : wide(ToM128D(x)) {}

template<typename D>
__m128d ToM128D(Quantity<D> const x) {
  return _mm_set1_pd(x.magnitude_);
}

template<typename Q>
__m128d ToM128D(Quantity128<Q> const x) {
  return x.wide;
}

template<typename T, typename>
__m128d ToM128D(T const x) {
  return _mm_set1_pd(static_cast<double>(x));
}

template<typename LScalar, typename RScalar, typename>
typename internal_generators::ProductGenerator<LScalar, RScalar>::Type
operator*(Quantity128<LScalar> const& left, RScalar const& right) {
  static_assert(!std::is_same<LScalar, LScalar>::value,
                "Should not call operator*");
  return
      typename internal_generators::ProductGenerator<LScalar, RScalar>::Type();
}

template<typename LScalar, typename RScalar, typename>
typename internal_generators::ProductGenerator<LScalar, RScalar>::Type
operator*(LScalar const& left, Quantity128<RScalar> const& right) {
  static_assert(!std::is_same<LScalar, LScalar>::value,
                "Should not call operator*");
  return
      typename internal_generators::ProductGenerator<LScalar, RScalar>::Type();
}

template<typename LScalar, typename RScalar, typename>
typename internal_generators::QuotientGenerator<LScalar, RScalar>::Type
operator/(LScalar const& left, Quantity128<RScalar> const& right) {
  static_assert(!std::is_same<LScalar, LScalar>::value,
                "Should not call operator/");
  return
      typename internal_generators::QuotientGenerator<LScalar, RScalar>::Type();
}

}  // namespace internal_wide
}  // namespace quantities
}  // namespace principia
