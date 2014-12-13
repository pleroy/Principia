﻿
#include "ksp_plugin/plugin.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>

#include "base/heap_checker.hpp"
#include "geometry/permutation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/mock_n_body_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system.hpp"

using principia::geometry::Bivector;
using principia::geometry::Permutation;
using principia::physics::MockNBodySystem;
using principia::quantities::Abs;
using principia::quantities::ArcTan;
using principia::quantities::Sqrt;
using principia::si::Day;
using principia::si::Hour;
using principia::si::Minute;
using principia::si::Radian;
using principia::si::AstronomicalUnit;
using principia::testing_utilities::AbsoluteError;
using principia::testing_utilities::RelativeError;
using principia::testing_utilities::AlmostEquals;
using principia::testing_utilities::ICRFJ2000Ecliptic;
using principia::testing_utilities::SolarSystem;
using testing::AllOf;
using testing::Eq;
using testing::Ge;
using testing::Gt;
using testing::InSequence;
using testing::Le;
using testing::Lt;
using testing::Ref;
using testing::SizeIs;
using testing::StrictMock;
using testing::_;

namespace principia {
namespace ksp_plugin {

namespace {

int const kNotABody = 1729;

// Appends a |DegreesOfFreedom| equal to the last one at the given |time| to
// each |Trajectory| in the |k|th parameter of the expected call.
// This parameter must be an |NBodySystem<Barycentric>::Trajectories|, |time|
// must be an |Instant|.
ACTION_TEMPLATE(AppendTimeToTrajectories,
                HAS_1_TEMPLATE_PARAMS(int, k),
                AND_1_VALUE_PARAMS(time)) {
  for (auto* trajectory : static_cast<NBodySystem<Barycentric>::Trajectories>(
                              std::tr1::get<k>(args))) {
    trajectory->Append(time, trajectory->last().degrees_of_freedom());
  }
}

}  // namespace

class TestablePlugin : public Plugin {
 public:
  // Takes ownership of |n_body_system|.
  // TODO(phl): We'd like to pass a |unique_ptr<>| but that seems to confuse
  // gMock.
  TestablePlugin(Instant const& initial_time,
                 Index const sun_index,
                 GravitationalParameter const& sun_gravitational_parameter,
                 Angle const& planetarium_rotation,
                 MockNBodySystem<Barycentric>* n_body_system)
      : Plugin(initial_time,
               sun_index,
               sun_gravitational_parameter,
               planetarium_rotation) {
    n_body_system_.reset(n_body_system);
  }

  Time const& Δt() const {
    return Δt_;
  }

  SPRKIntegrator<Length, Speed> const& prolongation_integrator() const {
    return prolongation_integrator_;
  }

  SPRKIntegrator<Length, Speed> const& history_integrator() const {
    return history_integrator_;
  }
};

class PluginTest : public testing::Test {
 protected:
  PluginTest()
      : looking_glass_(Permutation<ICRFJ2000Ecliptic, AliceSun>::XZY),
        solar_system_(SolarSystem::AtСпутник1Launch(
            SolarSystem::Accuracy::kMajorBodiesOnly)),
        bodies_(solar_system_->massive_bodies()),
        initial_time_(solar_system_->trajectories().front()->last().time()),
        sun_gravitational_parameter_(
            bodies_[SolarSystem::kSun]->gravitational_parameter()),
        planetarium_rotation_(1 * Radian) {
    satellite_initial_displacement_ =
        Displacement<AliceSun>({3111.0 * Kilo(Metre),
                                4400.0 * Kilo(Metre),
                                3810.0 * Kilo(Metre)});
    auto const tangent =
        satellite_initial_displacement_ * Bivector<double, AliceSun>({1, 2, 3});
    Vector<double, AliceSun> const unit_tangent = Normalize(tangent);
    EXPECT_THAT(
        InnerProduct(unit_tangent,
                     satellite_initial_displacement_ /
                         satellite_initial_displacement_.Norm()),
        Eq(0));
    // This yields a circular orbit.
    satellite_initial_velocity_ =
        Sqrt(bodies_[SolarSystem::kEarth]->gravitational_parameter() /
                 satellite_initial_displacement_.Norm()) * unit_tangent;

    n_body_system_ = new MockNBodySystem<Barycentric>();
    plugin_ = std::make_unique<StrictMock<TestablePlugin>>(
                  initial_time_,
                  SolarSystem::kSun,
                  sun_gravitational_parameter_,
                  planetarium_rotation_,
                  n_body_system_);
  }

  void InsertAllSolarSystemBodies() {
    for (std::size_t index = SolarSystem::kSun + 1;
         index < bodies_.size();
         ++index) {
      Index const parent_index = SolarSystem::parent(index);
      Displacement<AliceSun> const from_parent_position = looking_glass_(
          solar_system_->trajectories()[index]->
              last().degrees_of_freedom().position -
          solar_system_->trajectories()[parent_index]->
              last().degrees_of_freedom().position);
      Velocity<AliceSun> const from_parent_velocity = looking_glass_(
          solar_system_->trajectories()[index]->
              last().degrees_of_freedom().velocity -
          solar_system_->trajectories()[parent_index]->
              last().degrees_of_freedom().velocity);
      plugin_->InsertCelestial(index,
                               bodies_[index]->gravitational_parameter(),
                               parent_index,
                               from_parent_position,
                               from_parent_velocity);
    }
  }

  // The time of the |step|th synchronized history step of |plugin_|.
  // |HistoryTime(0)| is |initial_time_|.
  Instant HistoryTime(int const step) {
    return initial_time_ + step * plugin_->Δt();
  }

  // Keeps the vessel with the given |guid| during the next call to
  // |AdvanceTime|.  The vessel must be present.
  void KeepVessel(GUID const& guid) {
    bool const inserted = plugin_->InsertOrKeepVessel(guid,
                                                      SolarSystem::kEarth);
    EXPECT_FALSE(inserted) << guid;
  }

  // Inserts a vessel with the given |guid| and makes it a satellite of Earth
  // with relative position |satellite_initial_displacement_| and velocity
  // |satellite_initial_velocity_|.  The vessel must not already be present.
  // Increments the counter |*number_of_new_vessels|.  |number_of_new_vessels|
  // must not be null.
  void InsertVessel(GUID const& guid,
                    std::size_t* const number_of_new_vessels) {
    bool const inserted = plugin_->InsertOrKeepVessel(guid,
                                                      SolarSystem::kEarth);
    EXPECT_TRUE(inserted) << guid;
    plugin_->SetVesselStateOffset(guid,
                                  satellite_initial_displacement_,
                                  satellite_initial_velocity_);
    ++*CHECK_NOTNULL(number_of_new_vessels);
  }

  Permutation<ICRFJ2000Ecliptic, AliceSun> looking_glass_;
  MockNBodySystem<Barycentric>* n_body_system_;  // Not owned.
  std::unique_ptr<SolarSystem> solar_system_;
  SolarSystem::Bodies bodies_;
  Instant initial_time_;
  GravitationalParameter sun_gravitational_parameter_;
  Angle planetarium_rotation_;

  std::unique_ptr<StrictMock<TestablePlugin>> plugin_;

  // These initial conditions will yield a low circular orbit around Earth.
  Displacement<AliceSun> satellite_initial_displacement_;
  Velocity<AliceSun> satellite_initial_velocity_;

 private:
  base::HeapChecker heap_checker_;
};

using PluginDeathTest = PluginTest;

TEST_F(PluginTest, Initialization) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    Index const parent_index = SolarSystem::parent(index);
    EXPECT_THAT(solar_system_->trajectories()[index]->
                    last().degrees_of_freedom().position -
                solar_system_->trajectories()[parent_index]->
                    last().degrees_of_freedom().position,
                AlmostEquals(looking_glass_.Inverse()(
                    plugin_->CelestialDisplacementFromParent(index)),
                    1, 216320));
    EXPECT_THAT(solar_system_->trajectories()[index]->
                    last().degrees_of_freedom().velocity -
                solar_system_->trajectories()[parent_index]->
                    last().degrees_of_freedom().velocity,
                AlmostEquals(looking_glass_.Inverse()(
                    plugin_->CelestialParentRelativeVelocity(index)),
                    0, 936));
  }
}

TEST_F(PluginDeathTest, InsertCelestialError) {
  Displacement<AliceSun> const from_parent_position = looking_glass_(
      solar_system_->trajectories().front()->
          last().degrees_of_freedom().position -
      solar_system_->trajectories().front()->
          last().degrees_of_freedom().position);
  Velocity<AliceSun> const from_parent_velocity = looking_glass_(
      solar_system_->trajectories().front()->
          last().degrees_of_freedom().velocity -
      solar_system_->trajectories().front()->
          last().degrees_of_freedom().velocity);
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->InsertCelestial(42,
                             bodies_.front()->gravitational_parameter(),
                             SolarSystem::kSun,
                             from_parent_position,
                             from_parent_velocity);
  }, "before the end of initialization");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->InsertCelestial(42,
                             bodies_.front()->gravitational_parameter(),
                             kNotABody,
                             from_parent_position,
                             from_parent_velocity);
  }, "No body at index");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->InsertCelestial(SolarSystem::kEarth,
                             bodies_.front()->gravitational_parameter(),
                             SolarSystem::kSun,
                             from_parent_position,
                             from_parent_velocity);
  }, "Body already exists");
}

TEST_F(PluginDeathTest, UpdateCelestialHierarchyError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->UpdateCelestialHierarchy(SolarSystem::kSun, SolarSystem::kPluto);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->UpdateCelestialHierarchy(kNotABody, SolarSystem::kPluto);
  }, "No body at index");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->UpdateCelestialHierarchy(SolarSystem::kSun, kNotABody);
  }, "No body at index");
}

TEST_F(PluginDeathTest, InsertOrKeepVesselError) {
  GUID const guid = "Syrio Forel";
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->InsertOrKeepVessel(guid, SolarSystem::kSun);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->InsertOrKeepVessel(guid, kNotABody);
  }, "No body at index");
}

TEST_F(PluginDeathTest, SetVesselStateOffsetError) {
  GUID const guid = "Test Satellite";
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->SetVesselStateOffset(guid,
                                  satellite_initial_displacement_,
                                  satellite_initial_velocity_);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->SetVesselStateOffset(guid,
                                  satellite_initial_displacement_,
                                  satellite_initial_velocity_);
  }, "No vessel with GUID");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->InsertOrKeepVessel(guid, SolarSystem::kSun);
    plugin_->SetVesselStateOffset(guid,
                                  satellite_initial_displacement_,
                                  satellite_initial_velocity_);
    plugin_->SetVesselStateOffset(guid,
                                  satellite_initial_displacement_,
                                  satellite_initial_velocity_);
  }, "already has a trajectory");
}

TEST_F(PluginDeathTest, AdvanceTimeError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->AdvanceTime(Instant(), Angle());
  }, "Check failed: !initializing");
}

TEST_F(PluginDeathTest, VesselDisplacementFromParentError) {
  GUID const guid = "Test Satellite";
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->VesselDisplacementFromParent(guid);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->VesselDisplacementFromParent(guid);
  }, "No vessel with GUID");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->InsertOrKeepVessel(guid, SolarSystem::kSun);
    plugin_->VesselDisplacementFromParent(guid);
  }, "not given an initial state");
}

TEST_F(PluginDeathTest, VesselParentRelativeVelocityError) {
  GUID const guid = "Test Satellite";
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->VesselParentRelativeVelocity(guid);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->VesselParentRelativeVelocity(guid);
  }, "No vessel with GUID");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->InsertOrKeepVessel(guid, SolarSystem::kSun);
    plugin_->VesselParentRelativeVelocity(guid);
  }, "not given an initial state");
}

TEST_F(PluginDeathTest, CelestialDisplacementFromParentError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->CelestialDisplacementFromParent(SolarSystem::kEarth);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->CelestialDisplacementFromParent(kNotABody);
  }, "No body at index");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->CelestialDisplacementFromParent(SolarSystem::kSun);
  }, "is the sun");
}

TEST_F(PluginDeathTest, CelestialParentRelativeVelocityError) {
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->CelestialParentRelativeVelocity(SolarSystem::kEarth);
  }, "Check failed: !initializing");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->CelestialParentRelativeVelocity(kNotABody);
  }, "No body at index");
  EXPECT_DEATH({
    InsertAllSolarSystemBodies();
    plugin_->EndInitialization();
    plugin_->CelestialParentRelativeVelocity(SolarSystem::kSun);
  }, "is the sun");
}

TEST_F(PluginTest, VesselInsertionAtInitialization) {
  GUID const guid = "Test Satellite";
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  bool const inserted = plugin_->InsertOrKeepVessel(guid,
                                                    SolarSystem::kEarth);
  EXPECT_TRUE(inserted);
  plugin_->SetVesselStateOffset(guid,
                                satellite_initial_displacement_,
                                satellite_initial_velocity_);
  EXPECT_THAT(plugin_->VesselDisplacementFromParent(guid),
              AlmostEquals(satellite_initial_displacement_, 7460));
  EXPECT_THAT(plugin_->VesselParentRelativeVelocity(guid),
              AlmostEquals(satellite_initial_velocity_, 3));
}

// Checks that the plugin correctly uses its 10-second-step history even when
// advanced with smaller timesteps.
TEST_F(PluginTest, AdvanceTimeWithCelestialsOnly) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  Time const δt = 0.02 * Second;
  Angle const planetarium_rotation = 42 * Radian;
  for (int step = 0; step < 10; ++step) {
    for (Instant t = HistoryTime(step) + δt;
         t <= HistoryTime(step + 1);
         t += δt) {
      // Called to compute the prolongations.
      EXPECT_CALL(*n_body_system_,
                  Integrate(Ref(plugin_->prolongation_integrator()), t,
                            plugin_->Δt(), 0, true, SizeIs(bodies_.size())))
          .RetiresOnSaturation();
      plugin_->AdvanceTime(t, planetarium_rotation);
    }
    // Called to advance the synchronized histories.
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->history_integrator()),
                          HistoryTime(step + 1) + δt,
                          plugin_->Δt(), 0, false,
                          SizeIs(bodies_.size())))
        .WillOnce(AppendTimeToTrajectories<5>(HistoryTime(step + 1)))
        .RetiresOnSaturation();
    // Called to compute the prolongations.
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->prolongation_integrator()),
                          HistoryTime(step + 1) + δt, plugin_->Δt(), 0, true,
                          SizeIs(bodies_.size())))
        .RetiresOnSaturation();
    plugin_->AdvanceTime(HistoryTime(step + 1) + δt, planetarium_rotation);
  }
}

// Checks that the plugin correctly advances the history of newly inserted
// vessels with the prolongation integrator (using small steps), then switches
// to the history integrator.  Also tests the removal of vessels.
// The sequence of additions and removals is as follows (using the constants
// defined below; recall that |initial_time_ == HistoryTime(0)|).
// * |HistoryTime(0)|               : Insert |enterprise|.
// * |HistoryTime(0) + a_while|     : Insert |enterprise_d|.
// * |HistoryTime(1) + a_while|     : Insert |stargazer| and |bradbury|.
// * |HistoryTime(1) + half_a_step| : Remove |bradbury| (it is present at
//                                    |HistoryTime(1) + half_a_step|, but not at
//                                    |HistoryTime(1) + half_a_step + δt|).
//                                    It is thus removed while unsynchronized.
// * |HistoryTime(2)|               : Insert |constantinople|.
// * |HistoryTime(3)|               : Remove |enterprise|.
TEST_F(PluginTest, AdvanceTimeWithVessels) {
  Time const δt = 0.02 * Second;
  Time const a_while = 10 * δt;
  Time const half_a_step = 250 * δt;
  Time const ε_δt = 0.1 * δt;
  EXPECT_THAT(half_a_step, Eq(plugin_->Δt() / 2));
  GUID const enterprise = "NCC-1701";
  GUID const enterprise_d = "NCC-1701-D";
  GUID const stargazer = "NCC-2893";
  GUID const bradbury = "NX-72307";
  GUID const constantinople = "NCC-43622";
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  Angle const planetarium_rotation = 42 * Radian;
  std::size_t expected_number_of_new_vessels = 0;
  std::size_t expected_number_of_old_vessels = 0;
  InsertVessel(enterprise, &expected_number_of_new_vessels);
  for (int step = 0; step < 10; ++step) {
    for (Instant t = HistoryTime(step) + δt;
         t <= HistoryTime(step + 1);
         t += δt) {
      // Keep our vessels.  Make sure we're not inserting new ones.
      if (step <= 3) {
        KeepVessel(enterprise);
      }
      if (t > HistoryTime(0) + a_while + ε_δt) {
        KeepVessel(enterprise_d);
      }
      if (t > HistoryTime(1) + a_while + ε_δt) {
        KeepVessel(stargazer);
      }
      if (t > HistoryTime(1) + a_while + ε_δt &&
          t < HistoryTime(1) + half_a_step - ε_δt) {
        KeepVessel(bradbury);
      } else if (AbsoluteError(t - HistoryTime(1), half_a_step) < ε_δt) {
        // We will be removing |bradbury| in this step.
        --expected_number_of_new_vessels;
      }
      if (step > 2) {
        KeepVessel(constantinople);
      }
      // Called to compute the prolongations and advance the unsynchronized
      // histories.
      EXPECT_CALL(*n_body_system_,
                  Integrate(Ref(plugin_->prolongation_integrator()), t,
                            plugin_->Δt(), 0, true,
                            SizeIs(bodies_.size() +
                                       expected_number_of_old_vessels +
                                       expected_number_of_new_vessels)))
          .RetiresOnSaturation();
      plugin_->AdvanceTime(t, planetarium_rotation);
      if (AbsoluteError(t - HistoryTime(0), a_while) < ε_δt) {
        InsertVessel(enterprise_d, &expected_number_of_new_vessels);
      } else if (AbsoluteError(t - HistoryTime(1), a_while) < ε_δt) {
        InsertVessel(stargazer, &expected_number_of_new_vessels);
        InsertVessel(bradbury, &expected_number_of_new_vessels);
      }
    }
    // Keep the vessels for the history-advancing step.
    if (step <= 3) {
      KeepVessel(enterprise);
    }
    KeepVessel(enterprise_d);
    if (step >= 1) {
      KeepVessel(stargazer);
    }
    if (step > 2) {
      KeepVessel(constantinople);
    }
    // Called to advance the synchronized histories.
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->history_integrator()),
                          HistoryTime(step + 1) + δt,
                          plugin_->Δt(), 0, false,
                          SizeIs(bodies_.size() +
                                     expected_number_of_old_vessels)))
        .WillOnce(AppendTimeToTrajectories<5>(HistoryTime(step + 1)))
        .RetiresOnSaturation();
    if (expected_number_of_new_vessels > 0) {
      // Called to synchronize the new histories.
      EXPECT_CALL(*n_body_system_,
                  Integrate(Ref(plugin_->prolongation_integrator()),
                            HistoryTime(step + 1),
                            plugin_->Δt(), 0, true,
                            SizeIs(bodies_.size() +
                                       expected_number_of_new_vessels)))
          .WillOnce(AppendTimeToTrajectories<5>(HistoryTime(step + 1)))
          .RetiresOnSaturation();
    }
    // Called to compute the prolongations.
    EXPECT_CALL(*n_body_system_,
                Integrate(Ref(plugin_->prolongation_integrator()),
                          HistoryTime(step + 1) + δt, plugin_->Δt(), 0, true,
                          SizeIs(bodies_.size() +
                                     expected_number_of_old_vessels +
                                     expected_number_of_new_vessels)))
        .RetiresOnSaturation();
    plugin_->AdvanceTime(HistoryTime(step + 1) + δt, planetarium_rotation);
    expected_number_of_old_vessels += expected_number_of_new_vessels;
    expected_number_of_new_vessels = 0;
    if (step == 2) {
      InsertVessel(constantinople, &expected_number_of_new_vessels);
    } else if (step == 3) {
      // We will be removing |enterprise|.
      --expected_number_of_old_vessels;
    }
  }
}

TEST_F(PluginTest, UpdateCelestialHierarchy) {
  InsertAllSolarSystemBodies();
  plugin_->EndInitialization();
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    plugin_->UpdateCelestialHierarchy(index, SolarSystem::kSun);
  }
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    EXPECT_THAT(
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom().position -
        solar_system_->trajectories()[SolarSystem::kSun]->
            last().degrees_of_freedom().position,
        AlmostEquals(looking_glass_.Inverse()(
            plugin_->CelestialDisplacementFromParent(index)), 1, 5056));
    EXPECT_THAT(
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom().velocity -
        solar_system_->trajectories()[SolarSystem::kSun]->
            last().degrees_of_freedom().velocity,
        AlmostEquals(looking_glass_.Inverse()(
            plugin_->CelestialParentRelativeVelocity(index)), 1, 936));
  }
}

TEST_F(PluginTest, BodyCentredNonrotatingRenderingIntegration) {
  GUID const satellite = "satellite";
  // This is an integration test, so we need a plugin that will actually
  // integrate.
  Plugin plugin(initial_time_,
                SolarSystem::kSun,
                sun_gravitational_parameter_,
                planetarium_rotation_);
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    Index const parent_index = SolarSystem::parent(index);
    Displacement<AliceSun> const from_parent_position = looking_glass_(
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom().position -
        solar_system_->trajectories()[parent_index]->
            last().degrees_of_freedom().position);
    Velocity<AliceSun> const from_parent_velocity = looking_glass_(
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom().velocity -
        solar_system_->trajectories()[parent_index]->
            last().degrees_of_freedom().velocity);
    plugin.InsertCelestial(index,
                           bodies_[index]->gravitational_parameter(),
                           parent_index,
                           from_parent_position,
                           from_parent_velocity);
  }
  plugin.EndInitialization();
  plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  plugin.SetVesselStateOffset(satellite,
                              satellite_initial_displacement_,
                              satellite_initial_velocity_);
  std::unique_ptr<RenderingFrame> const geocentric =
      plugin.NewBodyCentredNonRotatingFrame(SolarSystem::kEarth);
  // We'll check that our orbit is rendered as circular (actually, we only check
  // that it is rendered within a thin spherical shell around the Earth).
  Length perigee = std::numeric_limits<double>::infinity() * Metre;
  Length apogee = -std::numeric_limits<double>::infinity() * Metre;
  Permutation<AliceSun, World> const alice_sun_to_world =
      Permutation<AliceSun, World>(Permutation<AliceSun, World>::XZY);
  Time const δt_long = 10 * Minute;
#if !defined(_DEBUG)
  Time const δt_short = 0.02 * Second;
  Instant t = initial_time_ + δt_short;
  // Exercise #267 by having small time steps at the beginning of the trajectory
  // that are not synchronized with those of the Earth.
  for (; t < initial_time_ + δt_long; t += δt_short) {
    plugin.AdvanceTime(t,
                       1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  }
#else
  Instant t = initial_time_ + δt_long;
  plugin.AdvanceTime(t, 0 * Radian);  // This ensures the history is nonempty.
  plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  t += δt_long;
#endif
  for (; t < initial_time_ + 12 * Hour; t += δt_long) {
    plugin.AdvanceTime(t,
                       1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
    // We give the sun an arbitrary nonzero velocity in |World|.
    Position<World> const sun_world_position =
        World::origin + Velocity<World>(
            { 0.1 * AstronomicalUnit / Hour,
             -1.0 * AstronomicalUnit / Hour,
              0.0 * AstronomicalUnit / Hour}) * (t - initial_time_);
    RenderedTrajectory<World> const rendered_trajectory =
        plugin.RenderedVesselTrajectory(satellite,
                                        *geocentric,
                                        sun_world_position);
    Position<World> const earth_world_position =
        sun_world_position + alice_sun_to_world(
            plugin.CelestialDisplacementFromParent(SolarSystem::kEarth));
    for (auto const segment : rendered_trajectory) {
      Length const l_min =
          std::min((segment.begin - earth_world_position).Norm(),
                   (segment.end - earth_world_position).Norm());
      Length const l_max =
          std::max((segment.begin - earth_world_position).Norm(),
                   (segment.end - earth_world_position).Norm());
      perigee = std::min(perigee, l_min);
      apogee = std::max(apogee, l_max);
    }
    // Check continuity.
    for (std::size_t i = 0; i + 1 < rendered_trajectory.size(); ++i) {
      EXPECT_THAT(rendered_trajectory[i].end,
                  Eq(rendered_trajectory[i + 1].begin));
    }
    EXPECT_THAT(Abs(apogee - perigee), Lt(1.1 * Metre));
  }
}

TEST_F(PluginTest, BarycentricRotatingRenderingIntegration) {
  GUID const satellite = "satellite";
  // This is an integration test, so we need a plugin that will actually
  // integrate.
  Plugin plugin(initial_time_,
                SolarSystem::kSun,
                sun_gravitational_parameter_,
                planetarium_rotation_);
  for (std::size_t index = SolarSystem::kSun + 1;
       index < bodies_.size();
       ++index) {
    Index const parent_index = SolarSystem::parent(index);
    Displacement<AliceSun> const from_parent_position = looking_glass_(
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom().position -
        solar_system_->trajectories()[parent_index]->
            last().degrees_of_freedom().position);
    Velocity<AliceSun> const from_parent_velocity = looking_glass_(
        solar_system_->trajectories()[index]->
            last().degrees_of_freedom().velocity -
        solar_system_->trajectories()[parent_index]->
            last().degrees_of_freedom().velocity);
    plugin.InsertCelestial(index,
                           bodies_[index]->gravitational_parameter(),
                           parent_index,
                           from_parent_position,
                           from_parent_velocity);
  }
  plugin.EndInitialization();
  plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  // A vessel at the Lagrange point L₅.
  Displacement<ICRFJ2000Ecliptic> const from_the_earth_to_the_moon =
      solar_system_->trajectories()[SolarSystem::kMoon]->
          last().degrees_of_freedom().position -
      solar_system_->trajectories()[SolarSystem::kEarth]->
          last().degrees_of_freedom().position;
  Velocity<ICRFJ2000Ecliptic> const moon_velocity_wrt_earth =
      solar_system_->trajectories()[SolarSystem::kMoon]->
          last().degrees_of_freedom().velocity -
      solar_system_->trajectories()[SolarSystem::kEarth]->
          last().degrees_of_freedom().velocity;
  Displacement<ICRFJ2000Ecliptic> const from_the_earth_to_l5 =
      from_the_earth_to_the_moon / 2 -
          Normalize(moon_velocity_wrt_earth) *
              from_the_earth_to_the_moon.Norm() * Sqrt(3) / 2;
  Velocity<ICRFJ2000Ecliptic> const initial_velocity =
      Rotation<ICRFJ2000Ecliptic, ICRFJ2000Ecliptic>(
          π / 3 * Radian,
          Wedge(moon_velocity_wrt_earth, from_the_earth_to_the_moon))(
              moon_velocity_wrt_earth);
  plugin.SetVesselStateOffset(satellite,
                              looking_glass_(from_the_earth_to_l5),
                              looking_glass_(initial_velocity));
  std::unique_ptr<RenderingFrame> const earth_moon_barycentric =
      plugin.NewBarycentricRotatingFrame(SolarSystem::kEarth,
                                         SolarSystem::kMoon);
  Permutation<AliceSun, World> const alice_sun_to_world =
      Permutation<AliceSun, World>(Permutation<AliceSun, World>::XZY);
  Time const δt_long = 1 * Hour;
#if defined(_DEBUG)
  Time const duration = 1 * Day;
  Instant t = initial_time_ + δt_long;
#else
  Time const δt_short = 0.02 * Second;
  Time const duration = 20 * Day;
  Instant t = initial_time_ + δt_short;
  // Exercise #267 by having small time steps at the beginning of the trajectory
  // that are not synchronized with those of the Earth and Moon.
  for (; t < initial_time_ + δt_long; t += δt_short) {
    plugin.AdvanceTime(t,
                       1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  }
#endif
  for (; t < initial_time_ + duration; t += δt_long) {
    plugin.AdvanceTime(t,
                       1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
    plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  }
  plugin.AdvanceTime(t,
                     1 * Radian / Pow<2>(Minute) * Pow<2>(t - initial_time_));
  plugin.InsertOrKeepVessel(satellite, SolarSystem::kEarth);
  // We give the sun an arbitrary nonzero velocity in |World|.
  Position<World> const sun_world_position =
      World::origin + Velocity<World>(
          { 0.1 * AstronomicalUnit / Hour,
           -1.0 * AstronomicalUnit / Hour,
            0.0 * AstronomicalUnit / Hour}) * (t - initial_time_);
  RenderedTrajectory<World> const rendered_trajectory =
      plugin.RenderedVesselTrajectory(satellite,
                                      *earth_moon_barycentric,
                                      sun_world_position);
  Position<World> const earth_world_position =
      sun_world_position + alice_sun_to_world(
          plugin.CelestialDisplacementFromParent(SolarSystem::kEarth));
  Position<World> const moon_world_position =
      earth_world_position + alice_sun_to_world(
          plugin.CelestialDisplacementFromParent(SolarSystem::kMoon));
  Length const earth_moon =
      (moon_world_position - earth_world_position).Norm();
  for (auto const segment : rendered_trajectory) {
    Length const satellite_earth =
        (segment.begin - earth_world_position).Norm();
    Length const satellite_moon =
        (segment.begin - moon_world_position).Norm();
    EXPECT_THAT(RelativeError(earth_moon, satellite_earth), Lt(0.0907));
    EXPECT_THAT(RelativeError(earth_moon, satellite_moon), Lt(0.131));
    EXPECT_THAT(RelativeError(satellite_moon, satellite_earth), Lt(0.148));
  }
  // Check continuity.
  for (std::size_t i = 0; i + 1 < rendered_trajectory.size(); ++i) {
    EXPECT_THAT(rendered_trajectory[i].end,
                Eq(rendered_trajectory[i + 1].begin));
  }
#if !defined(_DEBUG)
  // Check that there are no spikes in the rendered trajectory, i.e., that three
  // consecutive points form a sufficiently flat triangle.  This tests issue
  // #256.
  for (std::size_t i = 0; i + 2 < rendered_trajectory.size(); ++i) {
    EXPECT_THAT(
        (rendered_trajectory[i].begin - rendered_trajectory[i + 1].end).Norm(),
        Gt(((rendered_trajectory[i].begin -
                 rendered_trajectory[i + 1].begin).Norm() +
            (rendered_trajectory[i].end -
                 rendered_trajectory[i + 1].end).Norm()) / 1.5)) << i;
  }
#endif
}

}  // namespace ksp_plugin
}  // namespace principia
