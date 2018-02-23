
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

template<typename LScalar, typename RDimensions>
typename internal_generators::ProductGenerator<LScalar,
                                               Quantity<RDimensions>>::Type
operator*(Quantity128<LScalar> const&, Quantity<RDimensions> const&) {
  static_assert(!std::is_same<LScalar, LScalar>::value,
                "Should not call operator*");
  return typename internal_generators::
      ProductGenerator<LScalar, Quantity<RDimensions>>::Type();
}

template<typename LDimensions, typename RScalar>
typename internal_generators::ProductGenerator<Quantity<LDimensions>,
                                               RScalar>::Type
operator*(Quantity<LDimensions> const&, Quantity128<RScalar> const&) {
  static_assert(!std::is_same<RScalar, RScalar>::value,
                "Should not call operator*");
  return typename internal_generators::ProductGenerator<Quantity<LDimensions>,
                                                        RScalar>::Type();
}

template<typename LDimensions, typename RScalar>
typename internal_generators::QuotientGenerator<Quantity<LDimensions>,
                                                RScalar>::Type
operator/(Quantity<LDimensions> const&, Quantity128<RScalar> const&) {
  static_assert(!std::is_same<RScalar, RScalar>::value,
                "Should not call operator/");
  return typename internal_generators::QuotientGenerator<Quantity<LDimensions>,
                                                         RScalar>::Type();
}

template<typename LScalar>
typename internal_generators::ProductGenerator<LScalar, double>::Type
operator*(Quantity128<LScalar> const&, double) {
  static_assert(!std::is_same<LScalar, LScalar>::value,
                "Should not call operator*");
  return
      typename internal_generators::ProductGenerator<LScalar, double>::Type();
}

template<typename RScalar>
typename internal_generators::ProductGenerator<double, RScalar>::Type
operator*(double, Quantity128<RScalar> const&) {
  static_assert(!std::is_same<RScalar, RScalar>::value,
                "Should not call operator*");
  return
      typename internal_generators::ProductGenerator<double, RScalar>::Type();
}

template<typename LDimensions, typename RScalar>
typename internal_generators::QuotientGenerator<double, RScalar>::Type
operator/(double, Quantity128<RScalar> const&) {
  static_assert(!std::is_same<LScalar, LScalar>::value,
                "Should not call operator/");
  return
      typename internal_generators::QuotientGenerator<double, RScalar>::Type();
}

}  // namespace internal_wide
}  // namespace quantities
}  // namespace principia
