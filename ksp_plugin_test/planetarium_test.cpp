﻿
#include "ksp_plugin/planetarium.hpp"

#include <random>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/perspective.hpp"
#include "geometry/rotation.hpp"
#include "gtest/gtest.h"
#include "physics/massive_body.hpp"
#include "physics/mock_continuous_trajectory.hpp"
#include "physics/mock_dynamic_frame.hpp"
#include "physics/mock_ephemeris.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/numbers.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/serialization.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using astronomy::InfinitePast;
using astronomy::InfiniteFuture;
using base::make_not_null_unique;
using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::Displacement;
using geometry::LinearMap;
using geometry::Perspective;
using geometry::RigidTransformation;
using geometry::Rotation;
using geometry::Vector;
using geometry::Velocity;
using physics::MassiveBody;
using physics::MockContinuousTrajectory;
using physics::MockDynamicFrame;
using physics::MockEphemeris;
using physics::RigidMotion;
using physics::RotatingBody;
using quantities::Cos;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Time;
using quantities::si::ArcMinute;
using quantities::si::Degree;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::ReadFromBinaryFile;
using testing_utilities::VanishesBefore;
using ::testing::_;
using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Le;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::SizeIs;

class PlanetariumTest : public ::testing::Test {
 protected:
  PlanetariumTest()
      :  // The camera is located as {0, 20, 0} and is looking along -y.
        perspective_(
            RigidTransformation<Navigation, Camera>(
                Navigation::origin + Displacement<Navigation>(
                                         {0 * Metre, 20 * Metre, 0 * Metre}),
                Camera::origin,
                Rotation<Navigation, Camera>(
                    Vector<double, Navigation>({1, 0, 0}),
                    Vector<double, Navigation>({0, 0, 1}),
                    Bivector<double, Navigation>({0, -1, 0}))
                    .Forget()),
            /*focal=*/5 * Metre),
        // A body of radius 1 m located at the origin.
        body_(MassiveBody::Parameters(1 * Kilogram),
              RotatingBody<Barycentric>::Parameters(
                  /*mean_radius=*/1 * Metre,
                  /*reference_angle=*/0 * Radian,
                  /*reference_instant=*/t0_,
                  /*angular_frequency=*/10 * Radian / Second,
                  /*ascension_of_pole=*/0 * Radian,
                  /*declination_of_pole=*/π / 2 * Radian)),
        bodies_({&body_}) {
    ON_CALL(plotting_frame_, t_min()).WillByDefault(Return(InfinitePast));
    ON_CALL(plotting_frame_, t_max()).WillByDefault(Return(InfiniteFuture));
    EXPECT_CALL(plotting_frame_, ToThisFrameAtTime(_))
        .WillRepeatedly(Return(RigidMotion<Barycentric, Navigation>(
            RigidTransformation<Barycentric, Navigation>::Identity(),
            AngularVelocity<Barycentric>(),
            Velocity<Barycentric>())));
    EXPECT_CALL(ephemeris_, bodies()).WillRepeatedly(ReturnRef(bodies_));
    EXPECT_CALL(ephemeris_, trajectory(_))
        .WillRepeatedly(Return(&continuous_trajectory_));
    EXPECT_CALL(continuous_trajectory_, EvaluatePosition(_))
        .WillRepeatedly(Return(Barycentric::origin));
  }

  not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>>
  NewCircularTrajectory(Time const& period,
                        Time const& step,
                        Time const& last) {
    auto discrete_trajectory =
        make_not_null_unique<DiscreteTrajectory<Barycentric>>();
    for (Time t; t <= last; t += step) {
      Angle const ɑ = 2 * π * t * Radian / period;
      DegreesOfFreedom<Barycentric> const degrees_of_freedom(
          Barycentric::origin +
              Displacement<Barycentric>({10 * Metre * Sin(ɑ),
                                         10 * Metre * Cos(ɑ),
                                         0 * Metre}),
          Velocity<Barycentric>({(20 * π / period) * Metre * Cos(ɑ),
                                 -(20 * π / period) * Metre * Sin(ɑ),
                                 0 * Metre / Second}));
      discrete_trajectory->Append(t0_ + t, degrees_of_freedom);
    }
    return discrete_trajectory;
  }

  Instant const t0_;
  Perspective<Navigation, Camera> const perspective_;
  MockDynamicFrame<Barycentric, Navigation> plotting_frame_;
  RotatingBody<Barycentric> const body_;
  std::vector<not_null<MassiveBody const*>> const bodies_;
  MockContinuousTrajectory<Barycentric> continuous_trajectory_;
  MockEphemeris<Barycentric> ephemeris_;
};

TEST_F(PlanetariumTest, PlotMethod0) {
  auto const discrete_trajectory = NewCircularTrajectory(/*period=*/10 * Second,
                                                         /*step=*/1 * Second,
                                                         /*last=*/10 * Second);

  // No dark area, infinite acuity, wide field of view.
  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1,
      /*angular_resolution=*/0 * Degree,
      /*field_of_view=*/90 * Degree);
  Planetarium planetarium(
      parameters, perspective_, &ephemeris_, &plotting_frame_);
  auto const rp2_lines =
      planetarium.PlotMethod0(discrete_trajectory->Begin(),
                              discrete_trajectory->End(),
                              t0_ + 10 * Second);

  // Because of the way the trajectory was constructed we have two lines which
  // meet in front of the camera and are separated by a hole behind the planet.
  EXPECT_THAT(rp2_lines, SizeIs(2));
  EXPECT_THAT(rp2_lines[0].front().x() - rp2_lines[1].back().x(),
              VanishesBefore(1 * Metre, 6));
  EXPECT_THAT(rp2_lines[0].back().x() - rp2_lines[1].front().x(),
              AlmostEquals(10.0 / Sqrt(399.0) * Metre, 48, 94));

  for (auto const& rp2_line : rp2_lines) {
    for (auto const& rp2_point : rp2_line) {
      // The following limit is obtained by elementary geometry by noticing that
      // the trajectory is viewed from the camera under an angle of π / 6.
      EXPECT_THAT(rp2_point.x(),
                  AllOf(Ge(-5.0 / Sqrt(3.0) * Metre),
                        Le(5.0 / Sqrt(3.0) * Metre)));
      EXPECT_THAT(rp2_point.y(), VanishesBefore(1 * Metre, 5, 13));
    }
  }
}

TEST_F(PlanetariumTest, PlotMethod1) {
  // A quarter of a circular trajectory around the origin, with many small
  // segments.
  auto const discrete_trajectory =
      NewCircularTrajectory(/*period=*/100'000 * Second,
                            /*step=*/1 * Second,
                            /*last=*/25'000 * Second);

  // No dark area, human visual acuity, wide field of view.
  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1,
      /*angular_resolution=*/0.4 * ArcMinute,
      /*field_of_view=*/90 * Degree);
  Planetarium planetarium(
      parameters, perspective_, &ephemeris_, &plotting_frame_);
  auto const rp2_lines =
      planetarium.PlotMethod1(discrete_trajectory->Begin(),
                              discrete_trajectory->End(),
                              t0_ + 10 * Second);

  EXPECT_THAT(rp2_lines, SizeIs(1));
  EXPECT_THAT(rp2_lines[0], SizeIs(4954));
  for (auto const& rp2_point : rp2_lines[0]) {
    EXPECT_THAT(rp2_point.x(),
                AllOf(Ge(0 * Metre),
                      Le(5.0 / Sqrt(3.0) * Metre)));
    EXPECT_THAT(rp2_point.y(), VanishesBefore(1 * Metre, 0, 14));
  }
}

TEST_F(PlanetariumTest, PlotMethod2) {
  // A quarter of a circular trajectory around the origin, with many small
  // segments.
  auto const discrete_trajectory =
      NewCircularTrajectory(/*period=*/100'000 * Second,
                            /*step=*/1 * Second,
                            /*last=*/25'000 * Second);

  // No dark area, human visual acuity, wide field of view.
  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1,
      /*angular_resolution=*/0.4 * ArcMinute,
      /*field_of_view=*/90 * Degree);
  Planetarium planetarium(
      parameters, perspective_, &ephemeris_, &plotting_frame_);
  auto const rp2_lines =
      planetarium.PlotMethod2(discrete_trajectory->Begin(),
                              discrete_trajectory->End(),
                              t0_ + 10 * Second);

  EXPECT_THAT(rp2_lines, SizeIs(1));
  EXPECT_THAT(rp2_lines[0], SizeIs(43));
  for (auto const& rp2_point : rp2_lines[0]) {
    EXPECT_THAT(rp2_point.x(),
                AllOf(Ge(0 * Metre),
                      Le((5.0 / Sqrt(3.0)) * Metre)));
    EXPECT_THAT(rp2_point.y(), VanishesBefore(1 * Metre, 0, 14));
  }
}

#if !defined(_DEBUG)
TEST_F(PlanetariumTest, RealSolarSystem) {
  serialization::DiscreteTrajectory discrete_trajectory_message;
  discrete_trajectory_message.ParseFromString(ReadFromBinaryFile(
      SOLUTION_DIR / "ksp_plugin_test" / "planetarium_trajectory.proto.bin"));
  auto discrete_trajectory =
      DiscreteTrajectory<Barycentric>::ReadFromMessage(
          discrete_trajectory_message, {});

  serialization::Ephemeris ephemeris_message;
  ephemeris_message.ParseFromString(ReadFromBinaryFile(
      SOLUTION_DIR / "ksp_plugin_test" / "planetarium_ephemeris.proto.bin"));
  auto ephemeris = Ephemeris<Barycentric>::ReadFromMessage(ephemeris_message);

  serialization::DynamicFrame plotting_frame_message;
  plotting_frame_message.ParseFromString(
      ReadFromBinaryFile(SOLUTION_DIR / "ksp_plugin_test" /
                         "planetarium_plotting_frame.proto.bin"));
  auto plotting_frame =
      NavigationFrame::ReadFromMessage(plotting_frame_message, ephemeris.get());

  serialization::AffineMap affine_map_message;
  affine_map_message.ParseFromString(ReadFromBinaryFile(
      SOLUTION_DIR / "ksp_plugin_test" / "planetarium_to_camera.proto.bin"));
  auto rigid_transformation =
      RigidTransformation<Navigation, Camera>::ReadFromMessage(
          affine_map_message);

  EXPECT_EQ(23423, discrete_trajectory->Size());

  Planetarium::Parameters parameters(
      /*sphere_radius_multiplier=*/1,
      /*angular_resolution=*/0.4 * ArcMinute,
      /*field_of_view=*/90 * Degree);
  Planetarium planetarium(parameters,
                          Perspective<Navigation, Camera>(rigid_transformation,
                                                          /*focal=*/1 * Metre),
                          ephemeris.get(),
                          plotting_frame.get());
  auto const rp2_lines =
      planetarium.PlotMethod2(discrete_trajectory->Begin(),
                              discrete_trajectory->End(),
                              discrete_trajectory->last().time());

  EXPECT_EQ(2, rp2_lines.size());
  EXPECT_EQ(2, rp2_lines[0].size());
  EXPECT_EQ(9, rp2_lines[1].size());
}
#endif

}  // namespace internal_planetarium
}  // namespace ksp_plugin
}  // namespace principia
