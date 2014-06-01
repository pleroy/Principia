#include "benchmarks/symplectic_partitioned_runge_kutta_integrator.hpp"

// .\Release\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=300                                                    // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2014/06/01-10:56:16
// Benchmark                           Time(ns)    CPU(ns) Iterations
// ------------------------------------------------------------------
// BM_SolveHarmonicOscillator        1158462428 1157407419         52                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1164918290 1163707460         52                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1161977190 1163107456         52                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1161909218 1160707440         52                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator        1163742772 1162207450         52                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator_mean   1162201980 1161427445        260                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator_stddev    2185080    2249814        260                                 1.37019e-013, 1.37057e-013  // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <algorithm>

#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

using principia::integrators::SPRKIntegrator;
using principia::quantities::Abs;
using principia::quantities::AngularFrequency;
using principia::quantities::Cos;
using principia::quantities::Length;
using principia::quantities::Max;
using principia::quantities::Momentum;
using principia::quantities::Sin;

namespace principia {
namespace benchmarks {

template<typename Position, typename Momentum>
void SolveHarmonicOscillatorAndComputeError(benchmark::State* state,
                                            Position* q_error,
                                            Momentum* p_error) {
  SPRKIntegrator::Solution<Position, Momentum> solution;

  SolveHarmonicOscillator<Position, Momentum>(&solution);

  state->PauseTiming();
  *q_error = 0 * Position::SIUnit();
  *p_error = 0 * Momentum::SIUnit();
  for (size_t i = 0; i < solution.time.quantities.size(); ++i) {
    *q_error = Max(
        *q_error,
        Abs(solution.position[0].quantities[i] -
            1 * Length::SIUnit() * Cos(1 * AngularFrequency::SIUnit() *
                                           solution.time.quantities[i])));
    *p_error = Max(
        *p_error,
        Abs(solution.momentum[0].quantities[i] +
            1 * Momentum::SIUnit() * Sin(1 * AngularFrequency::SIUnit() *
                                             solution.time.quantities[i])));
  }
  state->ResumeTiming();
}

static void BM_SolveHarmonicOscillator(
    benchmark::State& state) {  // NOLINT(runtime/references)
  Length   q_error;
  Momentum p_error;
  while (state.KeepRunning()) {
    SolveHarmonicOscillatorAndComputeError<Length, Momentum>(
        &state, &q_error, &p_error);
  }
  std::stringstream ss;
  ss << q_error << ", " << p_error;
  state.SetLabel(ss.str());
}
BENCHMARK(BM_SolveHarmonicOscillator);

}  // namespace benchmarks
}  // namespace principia
