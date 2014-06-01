#pragma once

#include<vector>

#include "quantities/quantities.hpp"

using principia::quantities::Product;

namespace principia {
namespace testing_utilities {

template<typename Position, typename Momentum>
inline void ComputeHarmonicOscillatorForce(
    Time const t,
    std::vector<Position> const& q,
    std::vector<Quotient<Momentum, Time>>* result) {
  (*result)[0] = -q[0] * Quotient<Momentum, Product<Position, Time>>::SIUnit();
}

template<typename Position, typename Momentum>
inline void ComputeHarmonicOscillatorVelocity(
    std::vector<Momentum> const& p,
    std::vector<Quotient<Position, Time>>* result) {
  (*result)[0] = p[0] * Quotient<Position, Product<Momentum, Time>>::SIUnit();
}

}  // namespace testing_utilities
}  // namespace principia
