
#pragma once

#include <type_traits>

#include "base/not_constructible.hpp"

namespace principia {
namespace quantities {

namespace internal_quantities {
template<typename D>
class Quantity;
}  // namespace internal_quantities

namespace internal_wide {
template<typename Q>
struct Quantity128;
}  // namespace internal_wide;

namespace internal_traits {

using base::not_constructible;
using internal_quantities::Quantity;
using internal_wide::Quantity128;

// A type trait for testing if a type is a quantity.
template<typename T>
struct is_quantity : std::is_arithmetic<T>, not_constructible {};
template<typename D>
struct is_quantity<Quantity<D>> : std::true_type, not_constructible {};
template<typename Q>
struct is_quantity<Quantity128<Q>> : std::true_type, not_constructible {};

// This doesn't quite work in VS2015.  Let's not use it for now.
template<typename T>
constexpr bool is_quantity_v = is_quantity<T>::value;

}  // namespace internal_traits

using internal_traits::is_quantity;

}  // namespace quantities
}  // namespace principia
