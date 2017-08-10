﻿
#include "physics/apsides.hpp"

#include <limits>
#include <map>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "quantities/astronomy.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace physics {
namespace internal_apsides {

using base::not_null;
using geometry::Displacement;
using geometry::Frame;
using geometry::Velocity;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::QuinlanTremaine1990Order12;
using quantities::GravitationalParameter;
using quantities::Pow;
using quantities::Sin;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::astronomy::SolarMass;
using quantities::constants::GravitationalConstant;
using quantities::si::AstronomicalUnit;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using ::testing::Eq;

class ApsidesTest : public ::testing::Test {
 protected:
  using World =
      Frame<serialization::Frame::TestTag, serialization::Frame::TEST1, true>;
};

TEST_F(ApsidesTest, ComputeApsidesDiscreteTrajectory) {
  Instant const t0;
  GravitationalParameter const μ = GravitationalConstant * SolarMass;
  auto const b = new MassiveBody(μ);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<World>> initial_state;
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b));
  initial_state.emplace_back(World::origin, Velocity<World>());

  Ephemeris<World>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0,
          5 * Milli(Metre),
          Ephemeris<World>::FixedStepParameters(
              QuinlanTremaine1990Order12<Position<World>>(),
              10 * Minute));

  Displacement<World> r(
      {1 * AstronomicalUnit, 2 * AstronomicalUnit, 3 * AstronomicalUnit});
  Length const r_norm = r.Norm();
  Velocity<World> v({4 * Kilo(Metre) / Second,
                     5 * Kilo(Metre) / Second,
                     6 * Kilo(Metre) / Second});
  Speed const v_norm = v.Norm();

  Time const T = 2 * π * Sqrt(-(Pow<3>(r_norm) * Pow<2>(μ) /
                                Pow<3>(r_norm * Pow<2>(v_norm) - 2 * μ)));
  Length const a = -r_norm * μ / (r_norm * Pow<2>(v_norm) - 2 * μ);

  DiscreteTrajectory<World> trajectory;
  trajectory.Append(t0, DegreesOfFreedom<World>(World::origin + r, v));

  ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<World>::NoIntrinsicAcceleration,
      t0 + 10 * JulianYear,
      Ephemeris<World>::AdaptiveStepParameters(
          DormandElMikkawyPrince1986RKN434FM<Position<World>>(),
          std::numeric_limits<std::int64_t>::max(),
          1e-3 * Metre,
          1e-3 * Metre / Second),
      Ephemeris<World>::unlimited_max_ephemeris_steps,
      /*last_point_only=*/false);

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(*ephemeris.trajectory(b),
                 trajectory.Begin(),
                 trajectory.End(),
                 apoapsides,
                 periapsides);

  std::experimental::optional<Instant> previous_time;
  std::map<Instant, DegreesOfFreedom<World>> all_apsides;
  for (auto it = apoapsides.Begin(); it != apoapsides.End(); ++it) {
    Instant const time = it.time();
    all_apsides.emplace(time, it.degrees_of_freedom());
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(T, 118, 2079));
    }
    previous_time = time;
  }

  previous_time = std::experimental::nullopt;
  for (auto it = periapsides.Begin(); it != periapsides.End(); ++it) {
    Instant const time = it.time();
    all_apsides.emplace(time, it.degrees_of_freedom());
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(T, 143, 257));
    }
    previous_time = time;
  }

  EXPECT_EQ(6, all_apsides.size());

  previous_time = std::experimental::nullopt;
  std::experimental::optional<Position<World>> previous_position;
  for (auto const pair : all_apsides) {
    Instant const time = pair.first;
    Position<World> const position = pair.second.position();
    if (previous_time) {
      EXPECT_THAT(time - *previous_time,
                  AlmostEquals(0.5 * T, 103, 3567));
      EXPECT_THAT((position - *previous_position).Norm(),
                  AlmostEquals(2.0 * a, 0, 176));
    }
    previous_time = time;
    previous_position = position;
  }
}

TEST_F(ApsidesTest, ComputeNodes) {
  Instant const t0;
  GravitationalParameter const μ = GravitationalConstant * SolarMass;
  auto const b = new MassiveBody(μ);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<World>> initial_state;
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b));
  initial_state.emplace_back(World::origin, Velocity<World>());

  Ephemeris<World>
      ephemeris(
          std::move(bodies),
          initial_state,
          t0,
          5 * Milli(Metre),
          Ephemeris<World>::FixedStepParameters(
              QuinlanTremaine1990Order12<Position<World>>(),
              10 * Minute));

  KeplerianElements<World> elements;
  elements.eccentricity = 0.25;
  elements.semimajor_axis = 1 * AstronomicalUnit;
  elements.inclination = 10 * Degree;
  elements.longitude_of_ascending_node = 42 * Degree;
  elements.argument_of_periapsis = 100 * Degree;
  elements.mean_anomaly = 0 * Degree;
  KeplerOrbit<World> const orbit{
      *ephemeris.bodies()[0], MasslessBody{}, elements, t0};
  elements = orbit.elements_at_epoch();

  DiscreteTrajectory<World> trajectory;
  trajectory.Append(t0, initial_state[0] + orbit.StateVectors(t0));

  ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<World>::NoIntrinsicAcceleration,
      t0 + 10 * JulianYear,
      Ephemeris<World>::AdaptiveStepParameters(
          DormandElMikkawyPrince1986RKN434FM<Position<World>>(),
          std::numeric_limits<std::int64_t>::max(),
          1e-3 * Metre,
          1e-3 * Metre / Second),
      Ephemeris<World>::unlimited_max_ephemeris_steps,
      /*last_point_only=*/false);

  Vector<double, World> const north({0, 0, 1});

  DiscreteTrajectory<World> ascending_nodes;
  DiscreteTrajectory<World> descending_nodes;
  ComputeNodes(trajectory.Begin(),
               trajectory.End(),
               north,
               ascending_nodes,
               descending_nodes);

  std::experimental::optional<Instant> previous_time;
  for (auto it = ascending_nodes.Begin(); it != ascending_nodes.End(); ++it) {
    Instant const time = it.time();
    EXPECT_THAT((it.degrees_of_freedom().position() - World::origin)
                    .coordinates()
                    .ToSpherical()
                    .longitude,
                AlmostEquals(elements.longitude_of_ascending_node, 2, 100));
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(*elements.period, 0, 20));
    }
    previous_time = time;
  }

  previous_time = std::experimental::nullopt;
  for (auto it = descending_nodes.Begin(); it != descending_nodes.End(); ++it) {
    Instant const time = it.time();
    EXPECT_THAT(
        (it.degrees_of_freedom().position() - World::origin)
                .coordinates()
                .ToSpherical()
                .longitude,
        AlmostEquals(elements.longitude_of_ascending_node - π * Radian, 0, 25));
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(*elements.period, 0, 29));
    }
    previous_time = time;
  }

  EXPECT_THAT(ascending_nodes.Size(), Eq(10));
  EXPECT_THAT(descending_nodes.Size(), Eq(10));

  DiscreteTrajectory<World> south_ascending_nodes;
  DiscreteTrajectory<World> south_descending_nodes;
  Vector<double, World> const mostly_south({1, 1, -1});
  ComputeNodes(trajectory.Begin(),
               trajectory.End(),
               mostly_south,
               south_ascending_nodes,
               south_descending_nodes);
  EXPECT_THAT(south_ascending_nodes.Size(), Eq(10));
  EXPECT_THAT(south_descending_nodes.Size(), Eq(10));

  for (auto south_ascending_it  = south_ascending_nodes.Begin(),
            ascending_it        = ascending_nodes.Begin(),
            south_descending_it = south_descending_nodes.Begin(),
            descending_it       = descending_nodes.Begin();
       south_ascending_it != south_ascending_nodes.End();
       ++south_ascending_it,
       ++ascending_it,
       ++south_descending_it,
       ++descending_it) {
    EXPECT_THAT(south_ascending_it.degrees_of_freedom(),
                Eq(descending_it.degrees_of_freedom()));
    EXPECT_THAT(south_ascending_it.time(), Eq(descending_it.time()));
    EXPECT_THAT(south_descending_it.degrees_of_freedom(),
                Eq(ascending_it.degrees_of_freedom()));
    EXPECT_THAT(south_descending_it.time(), Eq(ascending_it.time()));
  }
}

}  // namespace internal_apsides
}  // namespace physics
}  // namespace principia
