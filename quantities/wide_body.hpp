
#pragma once

#include "quantities/wide.hpp"

namespace principia {
namespace quantities {
namespace internal_wide {

template<typename T>
Wide<T>::Wide(T const x) : wide_(ToM128D(x)) {}

template<typename D>
__m128d ToM128D(Quantity<D> const x) {
  return _mm_set1_pd(x.magnitude_);
}

template<typename T>
__m128d ToM128D(Wide<T> x) {
  return x.wide_;
}

template<typename T, typename>
inline __m128d ToM128D(T const x) {
  return _mm_set1_pd(static_cast<double>(x));
}

template<typename LScalar, typename RScalar>
typename internal_generators::ProductGenerator<LScalar, RScalar>::Type
operator*(Wide<LScalar> const& left, RScalar const& right) {
  static_assert(!std::is_same<LScalar, LScalar>::value,
                "Should not call operator*");
  return
      typename internal_generators::ProductGenerator<LScalar, RScalar>::Type();
}

template<typename LScalar, typename RScalar>
typename internal_generators::ProductGenerator<LScalar, RScalar>::Type
operator*(LScalar const& left, Wide<RScalar> const& right) {
  static_assert(!std::is_same<LScalar, LScalar>::value,
                "Should not call operator*");
  return
      typename internal_generators::ProductGenerator<LScalar, RScalar>::Type();
}

template<typename LScalar, typename RScalar>
typename internal_generators::ProductGenerator<LScalar, RScalar>::Type
operator/(LScalar const& left, Wide<RScalar> const& right) {
  static_assert(!std::is_same<LScalar, LScalar>::value,
                "Should not call operator/");
  return
      typename internal_generators::ProductGenerator<LScalar, RScalar>::Type();
}

}  // namespace internal_wide
}  // namespace quantities
}  // namespace principia
