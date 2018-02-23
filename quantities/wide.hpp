
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

// A wrapper for a quantity (or double) copied to all entries of a SIMD vector.
template<typename Q>
struct Quantity128 final {
  explicit Quantity128(Q x);
  __m128d const wide;
};

template<typename Source, typename Target>
using Wide = std::conditional_t<is_quantity<Target>::value,
                                Source,
                                Quantity128<Source>>;

// Fills both halves of the result.
template<typename D>
__m128d ToM128D(Quantity<D> x);
template<typename Q>
__m128d ToM128D(Quantity128<Q> x);
template<typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
__m128d ToM128D(T x);

// These operators are declared here to define Product and Quotient, but they
// are not expected to be called.

template<typename LScalar, typename RDimensions>
typename internal_generators::ProductGenerator<LScalar,
                                               Quantity<RDimensions>>::Type
operator*(Quantity128<LScalar> const&, Quantity<RDimensions> const&);

template<typename LDimensions, typename RScalar>
typename internal_generators::ProductGenerator<Quantity<LDimensions>,
                                               RScalar>::Type
operator*(Quantity<LDimensions> const&, Quantity128<RScalar> const&);

template<typename LDimensions, typename RScalar>
typename internal_generators::QuotientGenerator<Quantity<LDimensions>,
                                                RScalar>::Type
operator/(Quantity<LDimensions> const&, Quantity128<RScalar> const&);

template<typename LScalar>
typename internal_generators::ProductGenerator<LScalar, double>::Type
operator*(Quantity128<LScalar> const&, double);

template<typename RScalar>
typename internal_generators::ProductGenerator<double, RScalar>::Type
operator*(double, Quantity128<RScalar> const&);

template<typename LDimensions, typename RScalar>
typename internal_generators::QuotientGenerator<double, RScalar>::Type
operator/(double, Quantity128<RScalar> const&);

}  // namespace internal_wide

using internal_wide::ToM128D;
using internal_wide::Wide;

}  // namespace quantities
}  // namespace principia

// Because of circular dependencies, this file doesn't include wide_body.hpp.
// This will be done by quantities.hpp.
