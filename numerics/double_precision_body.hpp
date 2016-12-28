﻿
#pragma once

#include "numerics/double_precision.hpp"

#include "base/macros.hpp"
#include "geometry/serialization.hpp"

namespace principia {
namespace numerics {
namespace internal_double_precision {

using geometry::PointOrMultivectorSerializer;
using geometry::QuantityOrMultivectorSerializer;

template<typename T>
constexpr DoublePrecision<T>::DoublePrecision(T const& value)
    : value(value),
      error() {}

template<typename T>
FORCE_INLINE void DoublePrecision<T>::Increment(
    Difference<T> const& increment) {
  // The naming conventions follow Higham, Accuracy and Stability of Numerical
  // Algorithms, Algorithm 4.2.
  T const temp = value;
  Difference<T> const y = increment + error;
  value = temp + y;
  error = (temp - value) + y;
}

template<typename T>
void DoublePrecision<T>::Increment(
    DoublePrecision<Difference<T>> const& increment) {
  Increment(increment.value);
  Increment(increment.error);
}

template<typename T>
DoublePrecision<Difference<T>> DoublePrecision<T>::ToDifference() const {
  DoublePrecision<Difference<T>> result;
  result.value = value - T();
  result.error = error;
  return result;
}

template<typename T>
DoublePrecision<T> DoublePrecision<T>::Multiply(
    double const multiplier,
    DoublePrecision const & multiplicand) {
#ifdef _DEBUG
  int exponent;
  double const mantissa = std::frexp(multiplier, &exponent);
  CHECK_EQ(0.5, std::fabs(mantissa));
#endif
  DoublePrecision result;
  result.value = multiplicand.value * multiplier;
  result.error = multiplicand.error * multiplier;
  return result;
}

template<typename T>
void DoublePrecision<T>::WriteToMessage(
    not_null<serialization::DoublePrecision*> const message) const {
  using ValueSerializer = PointOrMultivectorSerializer<
                              T, serialization::DoublePrecision::Value>;
  using ErrorSerializer = QuantityOrMultivectorSerializer<
                              Difference<T>,
                              serialization::DoublePrecision::Error>;
  ValueSerializer::WriteToMessage(value, message->mutable_value());
  ErrorSerializer::WriteToMessage(error, message->mutable_error());
}

template<typename T>
DoublePrecision<T> DoublePrecision<T>::ReadFromMessage(
    serialization::DoublePrecision const& message) {
  using ValueSerializer = PointOrMultivectorSerializer<
                              T, serialization::DoublePrecision::Value>;
  using ErrorSerializer = QuantityOrMultivectorSerializer<
                              Difference<T>,
                              serialization::DoublePrecision::Error>;
  DoublePrecision double_precision;
  double_precision.value = ValueSerializer::ReadFromMessage(message.value());
  double_precision.error = ErrorSerializer::ReadFromMessage(message.error());
  return double_precision;
}

template<typename T>
std::ostream& operator<<(std::ostream& os,
                         const DoublePrecision<T>& double_precision) {
  os << double_precision.value << "|" << double_precision.error;
  return os;
}

}  // namespace internal_double_precision
}  // namespace numerics
}  // namespace principia
