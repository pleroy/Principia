#pragma once

#include <functional>
#include <utility>
#include <vector>

#include "base/array.hpp"
#include "geometry/hilbert.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_hermite3 {

using base::BoundedArray;
using quantities::Derivative;
using quantities::Difference;
using geometry::Hilbert;

// A 3rd degree Hermite polynomial defined by its values and derivatives at the
// bounds of some interval.
// TODO(phl): Invert the two template arguments for consistency with Derivative.
template<typename Argument, typename Value>
class Hermite3 final {
  using NormType = typename Hilbert<Difference<Value>>::NormType;
 public:
  using Derivative1 = Derivative<Value, Argument>;

  Hermite3(std::pair<Argument, Argument> arguments,
           std::pair<Value, Value> const& values,
           std::pair<Derivative1, Derivative1> const& derivatives);

  // NOTE(egg): this does not appear to use Casteljau's algorithm; perhaps it
  // should?
  Value Evaluate(Argument const& argument) const;
  Derivative1 EvaluateDerivative(Argument const& argument) const;

  // The result is sorted.
  BoundedArray<Argument, 2> FindExtrema() const;

  // |samples| must be a container; |get_argument| and |get_value| on the
  // elements of |samples| must return |Argument| and |Value| respectively
  // (possibly by reference or const-reference)
  // Returns the largest error (in the given |norm|) between this polynomial and
  // the given |samples|.
  template<typename Samples>
  NormType LInfinityError(
      Samples const& samples,
      std::function<Argument const&(typename Samples::value_type const&)> const&
          get_argument,
      std::function<Value const&(typename Samples::value_type const&)> const&
          get_value) const;

  // Returns true if the |LInfinityError| is less than |tolerance|.  More
  // efficient than the above function in the case where it returns false.
  template<typename Samples>
  bool LInfinityErrorIsWithin(
      Samples const& samples,
      std::function<Argument const&(typename Samples::value_type const&)> const&
          get_argument,
      std::function<Value const&(typename Samples::value_type const&)> const&
          get_value,
      NormType const& tolerance) const;

 private:
  using Derivative2 = Derivative<Derivative1, Argument>;
  using Derivative3 = Derivative<Derivative2, Argument>;

  std::pair<Argument, Argument> const arguments_;

  // The coefficients are relative to |argements.first|.
  Value a0_;
  Derivative1 a1_;
  Derivative2 a2_;
  Derivative3 a3_;
};

}  // namespace internal_hermite3

using internal_hermite3::Hermite3;

}  // namespace numerics
}  // namespace principia

#include "numerics/hermite3_body.hpp"
