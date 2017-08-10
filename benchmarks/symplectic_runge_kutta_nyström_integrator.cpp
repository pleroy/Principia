﻿
// .\Release\x64\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=5 --benchmark_filter=SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator                                                                                                                 // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <algorithm>
#include <functional>
#include <type_traits>
#include <vector>

#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/physics.pb.h"

namespace principia {

using geometry::Displacement;
using geometry::Frame;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using quantities::Abs;
using quantities::Acceleration;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Mass;
using quantities::Sin;
using quantities::SIUnit;
using quantities::Speed;
using quantities::Stiffness;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

namespace integrators {

namespace {

using World = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST, true>;

// TODO(egg): use the one from testing_utilities/integration again when everyone
// uses |Instant|s.
void ComputeHarmonicOscillatorAcceleration1D(
    Instant const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>& result) {
  result[0] = -q[0] * (SIUnit<Stiffness>() / SIUnit<Mass>());
}

void ComputeHarmonicOscillatorAcceleration3D(
    Instant const& t,
    std::vector<Position<World>> const& q,
    std::vector<Vector<Acceleration, World>>& result) {
  result[0] = (World::origin - q[0]) * (SIUnit<Stiffness>() / SIUnit<Mass>());
}

}  // namespace

template<typename Integrator>
void SolveHarmonicOscillatorAndComputeError1D(benchmark::State& state,
                                              Length& q_error,
                                              Speed& v_error,
                                              Integrator const& integrator) {
  using ODE = SpecialSecondOrderDifferentialEquation<Length>;

  state.PauseTiming();
  Length const q_initial = 1 * Metre;
  Speed const v_initial;
  Instant const t_initial;
#ifdef _DEBUG
  Instant const t_final = t_initial + 100 * Second;
#else
  Instant const t_final = t_initial + 1000 * Second;
#endif
  Time const step = 3.0e-4 * Second;

  std::vector<ODE::SystemState> solution;
  solution.reserve(static_cast<int>((t_final - t_initial) / step));
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration1D, _1, _2, _3);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{q_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.emplace_back(state);
  };

  auto const instance = integrator.NewInstance(problem, append_state, step);

  state.ResumeTiming();
  instance->Solve(t_final);
  state.PauseTiming();

  q_error = Length();
  v_error = Speed();
  for (auto const& state : solution) {
    q_error = std::max(q_error,
                       Abs(state.positions[0].value  -
                           Metre *
                           Cos((state.time.value - t_initial) *
                               (Radian / Second))));
    v_error = std::max(v_error,
                       Abs(state.velocities[0].value +
                           (Metre / Second) *
                           Sin((state.time.value - t_initial) *
                               (Radian / Second))));
  }
  state.ResumeTiming();
}

template<typename Integrator>
void SolveHarmonicOscillatorAndComputeError3D(benchmark::State& state,
                                              Length& q_error,
                                              Speed& v_error,
                                              Integrator const& integrator) {
  using ODE = SpecialSecondOrderDifferentialEquation<Position<World>>;

  state.PauseTiming();
  Displacement<World> const q_initial({1 * Metre, 0 * Metre, 0 * Metre});
  Velocity<World> const v_initial;
  Instant const t_initial;
#ifdef _DEBUG
  Instant const t_final = t_initial + 100 * Second;
#else
  Instant const t_final = t_initial + 1000 * Second;
#endif
  Time const step = 3.0e-4 * Second;

  std::vector<ODE::SystemState> solution;
  solution.reserve(static_cast<int>((t_final - t_initial) / step));
  ODE harmonic_oscillator;
  harmonic_oscillator.compute_acceleration =
      std::bind(ComputeHarmonicOscillatorAcceleration3D, _1, _2, _3);
  IntegrationProblem<ODE> problem;
  problem.equation = harmonic_oscillator;
  problem.initial_state = {{World::origin + q_initial}, {v_initial}, t_initial};
  auto const append_state = [&solution](ODE::SystemState const& state) {
    solution.emplace_back(state);
  };

  auto const instance = integrator.NewInstance(problem, append_state, step);

  state.ResumeTiming();
  instance->Solve(t_final);
  state.PauseTiming();

  q_error = Length();
  v_error = Speed();
  for (auto const& state : solution) {
    q_error = std::max(q_error,
                       ((state.positions[0].value - World::origin) -
                           q_initial *
                           Cos((state.time.value - t_initial) *
                               (Radian / Second))).Norm());
    v_error = std::max(v_error,
                       (state.velocities[0].value +
                           (q_initial / Second) *
                           Sin((state.time.value - t_initial) *
                               (Radian / Second))).Norm());
  }
  state.ResumeTiming();
}

template<typename Integrator, Integrator const& (*integrator)()>
void BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D(
    benchmark::State& state) {
  Length q_error;
  Speed v_error;
  while (state.KeepRunning()) {
    SolveHarmonicOscillatorAndComputeError1D(state, q_error, v_error,
                                             integrator());
  }
  std::stringstream ss;
  ss << q_error << ", " << v_error;
  state.SetLabel(ss.str());
}

template<typename Integrator, Integrator const& (*integrator)()>
void BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D(
    benchmark::State& state) {
  Length q_error;
  Speed v_error;
  while (state.KeepRunning()) {
    SolveHarmonicOscillatorAndComputeError3D(state, q_error, v_error,
                                             integrator());
  }
  std::stringstream ss;
  ss << q_error << ", " << v_error;
  state.SetLabel(ss.str());
}

BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    decltype(McLachlanAtela1992Order4Optimal<Length>()),
    &McLachlanAtela1992Order4Optimal<Length>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    decltype(McLachlan1995SB3A4<Length>()),
    &McLachlan1995SB3A4<Length>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    decltype(McLachlan1995SB3A5<Length>()),
    &McLachlan1995SB3A5<Length>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    decltype(BlanesMoan2002SRKN6B<Length>()),
    &BlanesMoan2002SRKN6B<Length>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    decltype(McLachlanAtela1992Order5Optimal<Length>()),
    &McLachlanAtela1992Order5Optimal<Length>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    decltype(OkunborSkeel1994Order6Method13<Length>()),
    &OkunborSkeel1994Order6Method13<Length>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    decltype(BlanesMoan2002SRKN11B<Length>()),
    &BlanesMoan2002SRKN11B<Length>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator1D,
    decltype(BlanesMoan2002SRKN14A<Length>()),
    &BlanesMoan2002SRKN14A<Length>);

BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    decltype(McLachlanAtela1992Order4Optimal<Position<World>>()),
    &McLachlanAtela1992Order4Optimal<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    decltype(McLachlan1995SB3A4<Position<World>>()),
    &McLachlan1995SB3A4<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    decltype(McLachlan1995SB3A5<Position<World>>()),
    &McLachlan1995SB3A5<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    decltype(BlanesMoan2002SRKN6B<Position<World>>()),
    &BlanesMoan2002SRKN6B<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    decltype(McLachlanAtela1992Order5Optimal<Position<World>>()),
    &McLachlanAtela1992Order5Optimal<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    decltype(OkunborSkeel1994Order6Method13<Position<World>>()),
    &OkunborSkeel1994Order6Method13<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    decltype(BlanesMoan2002SRKN11B<Position<World>>()),
    &BlanesMoan2002SRKN11B<Position<World>>);
BENCHMARK_TEMPLATE2(
    BM_SymplecticRungeKuttaNyströmIntegratorSolveHarmonicOscillator3D,
    decltype(BlanesMoan2002SRKN14A<Position<World>>()),
    &BlanesMoan2002SRKN14A<Position<World>>);

}  // namespace integrators
}  // namespace principia
