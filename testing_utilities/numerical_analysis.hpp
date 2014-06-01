#pragma once

#include <vector>

#include "quantities/quantities.hpp"

// Right-hand sides for various differential equations frequently used to test
// the properties of integrators.

using principia::quantities::Quotient;
using principia::quantities::Time;

namespace principia {
namespace testing_utilities {

// The one-dimensional unit harmonic oscillator,
// q' = p,  |ComputeHarmonicOscillatorVelocity|,
// p' = -q, |ComputeHarmonicOscillatorForce|.

template<typename Position, typename Momentum>
void ComputeHarmonicOscillatorVelocity(
    std::vector<Momentum> const& p,
    std::vector<Quotient<Position, Time>>* result);

template<typename Position, typename Momentum>
void ComputeHarmonicOscillatorForce(
    Time const t,
    std::vector<Position> const& q,
    std::vector<Quotient<Momentum, Time>>* result);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/numerical_analysis_body.hpp"
