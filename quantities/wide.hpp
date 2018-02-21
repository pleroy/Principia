
#pragma once

#include <nmmintrin.h>

#include <type_traits>

#include "base/not_constructible.hpp"
#include "quantities/generators.hpp"
#include "quantities/traits.hpp"

namespace principia {
namespace quantities {

namespace internal_quantities {
template<typename D>
class Quantity;
}  // namespace internal_quantities

namespace internal_wide {

using base::not_constructible;
using internal_quantities::Quantity;

// A wrapper for a quantity copied to all entries of a SIMD vector.
template<typename T>
class Wide final {
  static_assert(is_quantity<T>::value, "Not a quantity type");
 public:
  explicit Wide(T x);

 private:
  __m128d wide_;

  template<typename U>
  friend __m128d ToM128D(Wide<U> x);
};

// Fills both halves of the result.
template<typename D>
__m128d ToM128D(Quantity<D> x);
template<typename T>
__m128d ToM128D(Wide<T> x);
template<typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
__m128d ToM128D(T x);

// These operators are declared here to define Product and Quotient, but they
// are not expected to be called.

template<typename LScalar, typename RScalar>
typename internal_generators::ProductGenerator<LScalar, RScalar>::Type
operator*(Wide<LScalar> const&, RScalar const&);

template<typename LScalar, typename RScalar>
typename internal_generators::ProductGenerator<LScalar, RScalar>::Type
operator*(LScalar const&, Wide<RScalar> const&);

template<typename LScalar, typename RScalar>
typename internal_generators::ProductGenerator<LScalar, RScalar>::Type
operator/(LScalar const&, Wide<RScalar> const&);

}  // namespace internal_wide

using internal_wide::ToM128D;
using internal_wide::Wide;

}  // namespace quantities
}  // namespace principia

// Because of circular dependencies, this file doesn't include wide_body.hpp.
// This will be done by quantities.hpp.
