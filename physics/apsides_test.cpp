#include "physics/apsides.hpp"

#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/discrete_trajectory_factories.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.

namespace principia {
namespace physics {

using ::testing::Eq;
using ::testing::SizeIs;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::physics::_apsides;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_rotating_body;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_componentwise;
using namespace principia::testing_utilities::_discrete_trajectory_factories;
using namespace principia::testing_utilities::_is_near;

class ApsidesTest : public ::testing::Test {
 protected:
  using World = Frame<struct WorldTag, Inertial>;
};

#if !defined(_DEBUG)

TEST_F(ApsidesTest, ComputeApsidesDiscreteTrajectory) {
  Instant const t0;
  GravitationalParameter const μ = SolarGravitationalParameter;
  auto const b = new MassiveBody(μ);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<World>> initial_state;
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b));
  initial_state.emplace_back(World::origin, World::unmoving);

  Ephemeris<World> ephemeris(
      std::move(bodies),
      initial_state,
      t0,
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Metre,
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<World>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<World>::NewtonianMotionEquation>(),
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
  EXPECT_OK(trajectory.Append(t0,
                              DegreesOfFreedom<World>(World::origin + r, v)));

  EXPECT_OK(ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<World>::NoIntrinsicAcceleration,
      t0 + 10 * JulianYear,
      Ephemeris<World>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<World>::NewtonianMotionEquation>(),
          std::numeric_limits<std::int64_t>::max(),
          1e-3 * Metre,
          1e-3 * Metre / Second),
      Ephemeris<World>::unlimited_max_ephemeris_steps));

  DiscreteTrajectory<World> apoapsides;
  DiscreteTrajectory<World> periapsides;
  ComputeApsides(*ephemeris.trajectory(b),
                 trajectory,
                 trajectory.begin(),
                 trajectory.end(),
                 /*max_points=*/std::numeric_limits<int>::max(),
                 apoapsides,
                 periapsides);

  std::optional<Instant> previous_time;
  std::map<Instant, DegreesOfFreedom<World>> all_apsides;
  for (auto const& [time, degrees_of_freedom] : apoapsides) {
    all_apsides.emplace(time, degrees_of_freedom);
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(T, 118, 2824));
    }
    previous_time = time;
  }

  previous_time = std::nullopt;
  for (auto const& [time, degrees_of_freedom] : periapsides) {
    all_apsides.emplace(time, degrees_of_freedom);
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(T, 134, 257));
    }
    previous_time = time;
  }

  EXPECT_EQ(6, all_apsides.size());

  previous_time = std::nullopt;
  std::optional<Position<World>> previous_position;
  for (auto const& [time, degrees_of_freedom] : all_apsides) {
    Position<World> const position = degrees_of_freedom.position();
    if (previous_time) {
      EXPECT_THAT(time - *previous_time,
                  AlmostEquals(0.5 * T, 103, 5098));
      EXPECT_THAT((position - *previous_position).Norm(),
                  AlmostEquals(2.0 * a, 0, 176));
    }
    previous_time = time;
    previous_position = position;
  }
}

TEST_F(ApsidesTest, ComputeCollision) {
  Instant const t0;
  DiscreteTrajectory<World> reference_trajectory;
  DiscreteTrajectory<World> vessel_trajectory;

  // At |t0| the vessel is inside the celestial, so we expect the collision at a
  // negative time.
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({0 * Metre, -1 * Metre, 0 * Metre}),
              Velocity<World>({1 * Metre / Second,
                               0 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t0,
          /*t1=*/t0 - 10 * Second,
          /*t2=*/t0 + 10 * Second),
      reference_trajectory);
  AppendTrajectoryTimeline(
      NewLinearTrajectoryTimeline(
          DegreesOfFreedom<World>(
              World::origin +
                  Displacement<World>({1 * Metre, -1 * Metre, 0 * Metre}),
              Velocity<World>({0 * Metre / Second,
                               -1 * Metre / Second,
                               0 * Metre / Second})),
          /*Δt=*/1 * Second,
          t0,
          /*t1=*/t0 - 10 * Second,
          /*t2=*/t0 + 1 * Second),
      vessel_trajectory);

  RotatingBody<World> const body(
      1 * Kilogram,
      RotatingBody<World>::Parameters(
          /*min_radius=*/1 * Metre,
          /*mean_radius=*/2 * Metre,
          /*max_radius=*/3 * Metre,
          /*reference_angle=*/0 * Radian,
          /*reference_instant=*/t0,
          /*angular_frequency=*/2 * π * Radian / Second,
          /*right_ascension_of_pole=*/0 * Radian,
          /*declination_of_pole=*/π / 2 * Radian));

  // The celestial is infinite in the z direction and has four lobes in the x-y
  // plane.  Think of a LEGO® axle.
  auto radius = [](Angle const& latitude, Angle const& longitude) {
    return (Cos(4 * longitude) + 2) * Metre;
  };

  auto const maybe_collision = ComputeCollision(body,
                                                reference_trajectory,
                                                vessel_trajectory,
                                                vessel_trajectory.begin(),
                                                vessel_trajectory.end(),
                                                radius);
  ASSERT_TRUE(maybe_collision.has_value());
  auto const& collision = maybe_collision.value();

  // The collision was verified with Mathematica to the given accuracy.
  EXPECT_THAT(collision.time - t0,
              IsNear(-0.5254924180437539_(1) * Second));
  EXPECT_THAT(collision.degrees_of_freedom.position() - World::origin,
              Componentwise(1 * Metre,
                            IsNear(-0.4745075819562462_(1) * Metre),
                            0 * Metre));
  EXPECT_THAT(
      collision.degrees_of_freedom.velocity(),
      AlmostEquals(Velocity<World>({0 * Metre / Second,
                                    -1 * Metre / Second,
                                    0 * Metre / Second}), 0));
}

TEST_F(ApsidesTest, ComputeNodes) {
  Instant const t0;
  GravitationalParameter const μ = SolarGravitationalParameter;
  auto const b = new MassiveBody(μ);

  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<World>> initial_state;
  bodies.emplace_back(std::unique_ptr<MassiveBody const>(b));
  initial_state.emplace_back(World::origin, World::unmoving);

  Ephemeris<World> ephemeris(
      std::move(bodies),
      initial_state,
      t0,
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Metre,
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<World>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<World>::NewtonianMotionEquation>(),
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
  EXPECT_OK(trajectory.Append(t0, initial_state[0] + orbit.StateVectors(t0)));

  EXPECT_OK(ephemeris.FlowWithAdaptiveStep(
      &trajectory,
      Ephemeris<World>::NoIntrinsicAcceleration,
      t0 + 10 * JulianYear,
      Ephemeris<World>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<World>::NewtonianMotionEquation>(),
          std::numeric_limits<std::int64_t>::max(),
          1e-3 * Metre,
          1e-3 * Metre / Second),
      Ephemeris<World>::unlimited_max_ephemeris_steps));

  Vector<double, World> const north({0, 0, 1});

  DiscreteTrajectory<World> ascending_nodes;
  DiscreteTrajectory<World> descending_nodes;
  EXPECT_OK(ComputeNodes(trajectory,
                         trajectory.begin(),
                         trajectory.end(),
                         north,
                         /*max_points=*/std::numeric_limits<int>::max(),
                         ascending_nodes,
                         descending_nodes));

  std::optional<Instant> previous_time;
  for (auto const& [time, degrees_of_freedom] : ascending_nodes) {
    EXPECT_THAT((degrees_of_freedom.position() - World::origin)
                    .coordinates()
                    .ToSpherical()
                    .longitude,
                AlmostEquals(elements.longitude_of_ascending_node, 0, 104));
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(*elements.period, 0, 20));
    }
    previous_time = time;
  }

  previous_time = std::nullopt;
  for (auto const& [time, degrees_of_freedom] : descending_nodes) {
    EXPECT_THAT(
        (degrees_of_freedom.position() - World::origin)
                .coordinates()
                .ToSpherical()
                .longitude,
        AlmostEquals(elements.longitude_of_ascending_node - π * Radian, 0, 29));
    if (previous_time) {
      EXPECT_THAT(time - *previous_time, AlmostEquals(*elements.period, 0, 29));
    }
    previous_time = time;
  }

  EXPECT_THAT(ascending_nodes, SizeIs(10));
  EXPECT_THAT(descending_nodes, SizeIs(10));

  DiscreteTrajectory<World> south_ascending_nodes;
  DiscreteTrajectory<World> south_descending_nodes;
  Vector<double, World> const mostly_south({1, 1, -1});
  EXPECT_OK(ComputeNodes(trajectory,
                         trajectory.begin(),
                         trajectory.end(),
                         mostly_south,
                         /*max_points=*/std::numeric_limits<int>::max(),
                         south_ascending_nodes,
                         south_descending_nodes));
  EXPECT_THAT(south_ascending_nodes, SizeIs(10));
  EXPECT_THAT(south_descending_nodes, SizeIs(10));

  for (auto south_ascending_it  = south_ascending_nodes.begin(),
            ascending_it        = ascending_nodes.begin(),
            south_descending_it = south_descending_nodes.begin(),
            descending_it       = descending_nodes.begin();
       south_ascending_it != south_ascending_nodes.end();
       ++south_ascending_it,
       ++ascending_it,
       ++south_descending_it,
       ++descending_it) {
    EXPECT_THAT(south_ascending_it->degrees_of_freedom,
                Eq(descending_it->degrees_of_freedom));
    EXPECT_THAT(south_ascending_it->time, Eq(descending_it->time));
    EXPECT_THAT(south_descending_it->degrees_of_freedom,
                Eq(ascending_it->degrees_of_freedom));
    EXPECT_THAT(south_descending_it->time, Eq(ascending_it->time));
  }
}

#endif

}  // namespace physics
}  // namespace principia
