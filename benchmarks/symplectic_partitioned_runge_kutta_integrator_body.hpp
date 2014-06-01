#pragma once

#include <vector>

#undef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/numerical_analysis.hpp"

using principia::integrators::SPRKIntegrator;
using principia::quantities::Time;

namespace principia {
namespace benchmarks {

template<typename Position, typename Momentum>
inline void SolveHarmonicOscillator(
    SPRKIntegrator::Solution<Position, Momentum>* solution) {
  using principia::testing_utilities::ComputeHarmonicOscillatorForce;
  using principia::testing_utilities::ComputeHarmonicOscillatorVelocity;
  SPRKIntegrator integrator;
  SPRKIntegrator::Parameters<Position, Momentum> parameters;

  integrator.Initialize(integrator.Order5Optimal());

  parameters.q0 = {1.0 * Position::SIUnit()};
  parameters.p0 = {0.0 * Momentum::SIUnit()};
  parameters.t0 = 0.0 * Time::SIUnit();
#ifdef _DEBUG
  parameters.tmax = 100.0 * Time::SIUnit();
#else
  parameters.tmax = 1000.0 * Time::SIUnit();
#endif
  parameters.Δt = 1.0E-4 * Time::SIUnit();
  parameters.sampling_period = 1;
  integrator.Solve(&ComputeHarmonicOscillatorForce<Position, Momentum>,
                   &ComputeHarmonicOscillatorVelocity<Position, Momentum>,
                   parameters,
                   solution);
}

}  // namespace benchmarks
}  // namespace principia
