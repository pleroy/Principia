﻿
#pragma once

#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
namespace internal_double_precision {

using base::not_null;
using quantities::Difference;

// A simple container for accumulating a value using compensated summation.  The
// type of the value must be an affine space.  The value constructor is not
// explicit to make it easy to construct an object with no error.
template<typename T>
struct DoublePrecision final {
  constexpr DoublePrecision() = default;
  constexpr DoublePrecision(T const& value);  // NOLINT(runtime/explicit)

  void Increment(Difference<T> const& increment);
  void Increment(DoublePrecision<Difference<T>> const& increment);

  // This class doesn't make it easy to convert e.g., a |Position| to a
  // |Displacement|.
  DoublePrecision<Difference<T>> ToDifference() const;

  // |multiplier| must be a signed power of two.
  // TODO(phl): Only for vector spaces.
  static DoublePrecision Multiply(double const multiplier,
                                  DoublePrecision const& multiplicand);

  void WriteToMessage(not_null<serialization::DoublePrecision*> message) const;
  static DoublePrecision ReadFromMessage(
      serialization::DoublePrecision const& message);

  T value;
  Difference<T> error;
};

template<typename T>
std::ostream& operator<<(std::ostream& os,
                         const DoublePrecision<T>& double_precision);

}  // namespace internal_double_precision

using internal_double_precision::DoublePrecision;

}  // namespace numerics
}  // namespace principia

#include "numerics/double_precision_body.hpp"
