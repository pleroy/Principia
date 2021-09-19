﻿
#include "physics/discrete_trajectory.hpp"

#include <algorithm>
#include <functional>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/numerics_matchers.hpp"
#include "testing_utilities/trajectory_factories.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory {

using geometry::Displacement;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::Instant;
using geometry::Point;
using geometry::Position;
using geometry::R3Element;
using geometry::Vector;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Speed;
using quantities::Sin;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Micro;
using quantities::si::Milli;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteErrorFrom;
using testing_utilities::AppendTrajectory;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using testing_utilities::NewCircularTrajectory;
using testing_utilities::RelativeError;
using testing_utilities::operator""_⑴;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;
using ::testing::Contains;
using ::testing::Each;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Le;
using ::testing::Lt;
using ::testing::Pair;
using ::testing::Ref;

class DiscreteTrajectoryTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  DiscreteTrajectoryTest()
      : q1_(World::origin +
            Displacement<World>({1 * Metre, 2 * Metre, 3 * Metre})),
        q2_(World::origin +
            Displacement<World>({11 * Metre, 12 * Metre, 13 * Metre})),
        q3_(World::origin +
            Displacement<World>({21 * Metre, 22 * Metre, 23 * Metre})),
        q4_(World::origin +
            Displacement<World>({31 * Metre, 32 * Metre, 33 * Metre})),
        p1_(Velocity<World>({4 * Metre / Second,
                               5 * Metre / Second,
                               6 * Metre / Second})),
        p2_(Velocity<World>({14 * Metre / Second,
                               15 * Metre / Second,
                               16 * Metre / Second})),
        p3_(Velocity<World>({24 * Metre / Second,
                               25 * Metre / Second,
                               26 * Metre / Second})),
        p4_(Velocity<World>({34 * Metre / Second,
                               35 * Metre / Second,
                               36 * Metre / Second})),
        d1_(DegreesOfFreedom<World>(q1_, p1_)),
        d2_(DegreesOfFreedom<World>(q2_, p2_)),
        d3_(DegreesOfFreedom<World>(q3_, p3_)),
        d4_(DegreesOfFreedom<World>(q4_, p4_)) {
    t1_ = t0_ + 7 * Second;
    t2_ = t0_ + 17 * Second;
    t3_ = t0_ + 27 * Second;
    t4_ = t0_ + 37 * Second;

    massive_trajectory_ = std::make_unique<DiscreteTrajectory<World>>();
    massless_trajectory_ = std::make_unique<DiscreteTrajectory<World>>();
  }

  std::map<Instant, Position<World>> Positions(
      DiscreteTrajectory<World> const& trajectory) const {
    std::map<Instant, Position<World>> result;
    for (auto const& [time, degrees_of_freedom] : trajectory) {
      result.emplace_hint(result.end(), time, degrees_of_freedom.position());
    }
    return result;
  }

  std::map<Instant, Velocity<World>> Velocities(
      DiscreteTrajectory<World> const& trajectory) const {
    std::map<Instant, Velocity<World>> result;
    for (auto const& [time, degrees_of_freedom] : trajectory) {
      result.emplace_hint(result.end(), time, degrees_of_freedom.velocity());
    }
    return result;
  }

  std::list<Instant> Times(DiscreteTrajectory<World> const& trajectory) const {
    std::list<Instant> result;
    for (auto const& [time, degrees_of_freedom] : trajectory) {
      result.push_back(time);
    }
    return result;
  }

  Position<World> q1_, q2_, q3_, q4_;
  Velocity<World> p1_, p2_, p3_, p4_;
  DegreesOfFreedom<World> d1_, d2_, d3_, d4_;
  Instant const t0_;
  Instant t1_, t2_, t3_, t4_;
  std::unique_ptr<DiscreteTrajectory<World>> massive_trajectory_;
  std::unique_ptr<DiscreteTrajectory<World>> massless_trajectory_;
};

using DiscreteTrajectoryDeathTest = DiscreteTrajectoryTest;

TEST_F(DiscreteTrajectoryDeathTest, NewForkWithCopyError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    massive_trajectory_->Append(t3_, d3_);
    massive_trajectory_->NewForkWithCopy(t2_);
  }, "nonexistent time");
}

TEST_F(DiscreteTrajectoryTest, NewForkWithCopySuccess) {
  EXPECT_TRUE(massive_trajectory_->Empty());
  massive_trajectory_->Append(t1_, d1_);
  EXPECT_EQ(1, massive_trajectory_->Size());
  EXPECT_FALSE(massive_trajectory_->Empty());
  massive_trajectory_->Append(t2_, d2_);
  EXPECT_EQ(2, massive_trajectory_->Size());
  EXPECT_FALSE(massive_trajectory_->Empty());
  massive_trajectory_->Append(t3_, d3_);
  EXPECT_EQ(3, massive_trajectory_->Size());
  EXPECT_FALSE(massive_trajectory_->Empty());
  not_null<DiscreteTrajectory<World>*> const fork =
      massive_trajectory_->NewForkWithCopy(t2_);
  EXPECT_EQ(3, fork->Size());
  EXPECT_FALSE(fork->Empty());
  fork->Append(t4_, d4_);
  EXPECT_EQ(4, fork->Size());
  EXPECT_FALSE(fork->Empty());
  std::map<Instant, Position<World>> positions =
      Positions(*massive_trajectory_);
  std::map<Instant, Velocity<World>> velocities =
      Velocities(*massive_trajectory_);
  std::list<Instant> times = Times(*massive_trajectory_);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  positions = Positions(*fork);
  velocities = Velocities(*fork);
  times = Times(*fork);
  EXPECT_THAT(
      positions,
      ElementsAre(
          Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_), Pair(t4_, q4_)));
  EXPECT_THAT(
      velocities,
      ElementsAre(
          Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_), Pair(t4_, p4_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_, t4_));
}

TEST_F(DiscreteTrajectoryDeathTest, NewForkWithoutCopyError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    massive_trajectory_->Append(t3_, d3_);
    massive_trajectory_->NewForkWithoutCopy(t2_);
  }, "nonexistent time");
}

TEST_F(DiscreteTrajectoryTest, NewForkWithoutCopySuccess) {
  EXPECT_TRUE(massive_trajectory_->Empty());
  massive_trajectory_->Append(t1_, d1_);
  EXPECT_EQ(1, massive_trajectory_->Size());
  EXPECT_FALSE(massive_trajectory_->Empty());
  massive_trajectory_->Append(t2_, d2_);
  EXPECT_EQ(2, massive_trajectory_->Size());
  EXPECT_FALSE(massive_trajectory_->Empty());
  massive_trajectory_->Append(t3_, d3_);
  EXPECT_EQ(3, massive_trajectory_->Size());
  EXPECT_FALSE(massive_trajectory_->Empty());
  not_null<DiscreteTrajectory<World>*> const fork =
      massive_trajectory_->NewForkWithoutCopy(t2_);
  EXPECT_EQ(2, fork->Size());
  EXPECT_FALSE(fork->Empty());
  fork->Append(t4_, d4_);
  EXPECT_EQ(3, fork->Size());
  EXPECT_FALSE(fork->Empty());
  std::map<Instant, Position<World>> positions =
      Positions(*massive_trajectory_);
  std::map<Instant, Velocity<World>> velocities =
      Velocities(*massive_trajectory_);
  std::list<Instant> times = Times(*massive_trajectory_);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  positions = Positions(*fork);
  velocities = Velocities(*fork);
  times = Times(*fork);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t4_, q4_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t4_, p4_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t4_));
}

TEST_F(DiscreteTrajectoryTest, NewForkWithCopyAtLast) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork1 =
      massive_trajectory_->NewForkWithCopy(t3_);
  not_null<DiscreteTrajectory<World>*> const fork2 =
      fork1->NewForkWithCopy(fork1->back().time);
  not_null<DiscreteTrajectory<World>*> const fork3 =
      fork2->NewForkWithCopy(fork1->back().time);
  EXPECT_EQ(t3_, massive_trajectory_->back().time);
  EXPECT_EQ(t3_, fork1->back().time);

  std::map<Instant, Position<World>> positions = Positions(*fork2);
  std::map<Instant, Velocity<World>> velocities = Velocities(*fork2);
  std::list<Instant> times = Times(*fork2);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork2->back().degrees_of_freedom.position());
  EXPECT_EQ(p3_, fork2->back().degrees_of_freedom.velocity());
  EXPECT_EQ(t3_, fork2->back().time);

  std::vector<Instant> after;
  for (DiscreteTrajectory<World>::Iterator it = fork3->Find(t3_);
       it != fork3->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t3_));

  fork2->ForgetAfter(t3_);
  positions = Positions(*fork2);
  velocities = Velocities(*fork2);
  times = Times(*fork2);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork2->back().degrees_of_freedom.position());
  EXPECT_EQ(p3_, fork2->back().degrees_of_freedom.velocity());
  EXPECT_EQ(t3_, fork2->back().time);

  after.clear();
  for (DiscreteTrajectory<World>::Iterator it = fork2->Find(t3_);
       it != fork2->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t3_));

  fork1->Append(t4_, d4_);
  positions = Positions(*fork2);
  velocities = Velocities(*fork2);
  times = Times(*fork2);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork2->back().degrees_of_freedom.position());
  EXPECT_EQ(p3_, fork2->back().degrees_of_freedom.velocity());
  EXPECT_EQ(t3_, fork2->back().time);

  after.clear();
  for (DiscreteTrajectory<World>::Iterator it = fork1->Find(t3_);
       it != fork1->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  positions = Positions(*fork3);
  velocities = Velocities(*fork3);
  times = Times(*fork3);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork3->back().degrees_of_freedom.position());
  EXPECT_EQ(p3_, fork3->back().degrees_of_freedom.velocity());
  EXPECT_EQ(t3_, fork3->back().time);

  fork2->Append(t4_, d4_);
  after.clear();
  for (DiscreteTrajectory<World>::Iterator it = fork2->Find(t3_);
       it != fork2->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  fork3->Append(t4_, d4_);
  after.clear();
  for (DiscreteTrajectory<World>::Iterator it = fork3->Find(t3_);
       it != fork3->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  after.clear();
  for (DiscreteTrajectory<World>::Iterator it = fork3->Find(t2_);
       it != fork3->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t2_, t3_, t4_));
}

TEST_F(DiscreteTrajectoryTest, NewForkAtLast) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork1 =
      massive_trajectory_->NewForkAtLast();
  not_null<DiscreteTrajectory<World>*> const fork2 =
      fork1->NewForkAtLast();
  not_null<DiscreteTrajectory<World>*> const fork3 =
      fork2->NewForkAtLast();
  EXPECT_EQ(t3_, massive_trajectory_->back().time);
  EXPECT_EQ(t3_, fork1->back().time);

  std::map<Instant, Position<World>> positions = Positions(*fork2);
  std::map<Instant, Velocity<World>> velocities = Velocities(*fork2);
  std::list<Instant> times = Times(*fork2);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork2->back().degrees_of_freedom.position());
  EXPECT_EQ(p3_, fork2->back().degrees_of_freedom.velocity());
  EXPECT_EQ(t3_, fork2->back().time);

  std::vector<Instant> after;
  for (DiscreteTrajectory<World>::Iterator it = fork3->Find(t3_);
       it != fork3->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t3_));

  fork2->ForgetAfter(t3_);
  positions = Positions(*fork2);
  velocities = Velocities(*fork2);
  times = Times(*fork2);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork2->back().degrees_of_freedom.position());
  EXPECT_EQ(p3_, fork2->back().degrees_of_freedom.velocity());
  EXPECT_EQ(t3_, fork2->back().time);

  after.clear();
  for (DiscreteTrajectory<World>::Iterator it = fork2->Find(t3_);
       it != fork2->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t3_));

  fork1->Append(t4_, d4_);
  positions = Positions(*fork2);
  velocities = Velocities(*fork2);
  times = Times(*fork2);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork2->back().degrees_of_freedom.position());
  EXPECT_EQ(p3_, fork2->back().degrees_of_freedom.velocity());
  EXPECT_EQ(t3_, fork2->back().time);

  after.clear();
  for (DiscreteTrajectory<World>::Iterator it = fork1->Find(t3_);
       it != fork1->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  positions = Positions(*fork3);
  velocities = Velocities(*fork3);
  times = Times(*fork3);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
  EXPECT_EQ(q3_, fork3->back().degrees_of_freedom.position());
  EXPECT_EQ(p3_, fork3->back().degrees_of_freedom.velocity());
  EXPECT_EQ(t3_, fork3->back().time);

  fork2->Append(t4_, d4_);
  after.clear();
  for (DiscreteTrajectory<World>::Iterator it = fork2->Find(t3_);
       it != fork2->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  fork3->Append(t4_, d4_);
  after.clear();
  for (DiscreteTrajectory<World>::Iterator it = fork3->Find(t3_);
       it != fork3->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t3_, t4_));

  after.clear();
  for (DiscreteTrajectory<World>::Iterator it = fork3->Find(t2_);
       it != fork3->end();
       ++it) {
    after.push_back(it->time);
  }
  EXPECT_THAT(after, ElementsAre(t2_, t3_, t4_));
}

TEST_F(DiscreteTrajectoryTest, AttachFork) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);

  not_null<std::unique_ptr<DiscreteTrajectory<World>>> fork1 =
      make_not_null_unique<DiscreteTrajectory<World>>();
  fork1->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork2 =
      fork1->NewForkWithoutCopy(t3_);
  fork2->Append(t4_, d4_);

  auto const unowned_fork1 = fork1.get();
  massive_trajectory_->AttachFork(std::move(fork1));

  std::map<Instant, Position<World>> positions =
      Positions(*massive_trajectory_);
  std::map<Instant, Velocity<World>> velocities =
      Velocities(*massive_trajectory_);
  std::list<Instant> times = Times(*massive_trajectory_);
  EXPECT_THAT(positions, ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_)));
  EXPECT_THAT(velocities, ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_));

  positions = Positions(*unowned_fork1);
  velocities = Velocities(*unowned_fork1);
  times = Times(*unowned_fork1);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));

  positions = Positions(*fork2);
  velocities = Velocities(*fork2);
  times = Times(*fork2);
  EXPECT_THAT(
      positions,
      ElementsAre(
          Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_), Pair(t4_, q4_)));
  EXPECT_THAT(
      velocities,
      ElementsAre(
          Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_), Pair(t4_, p4_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_, t4_));
}

TEST_F(DiscreteTrajectoryTest, DetachFork) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork1 =
      massive_trajectory_->NewForkWithCopy(t2_);
  not_null<DiscreteTrajectory<World>*> const fork2 =
      fork1->NewForkWithoutCopy(t2_);
  fork1->Append(t4_, d4_);

  auto const detached1 = fork1->DetachFork();
  std::map<Instant, Position<World>> positions = Positions(*detached1);
  std::map<Instant, Velocity<World>> velocities = Velocities(*detached1);
  std::list<Instant> times = Times(*detached1);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t2_, q2_), Pair(t3_, q3_), Pair(t4_, q4_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t2_, p2_), Pair(t3_, p3_), Pair(t4_, p4_)));
  EXPECT_THAT(times, ElementsAre(t2_, t3_, t4_));

  auto const detached2 = fork2->DetachFork();
  positions = Positions(*detached2);
  velocities = Velocities(*detached2);
  times = Times(*detached2);
  EXPECT_THAT(positions, ElementsAre(Pair(t2_, q2_)));
  EXPECT_THAT(velocities, ElementsAre(Pair(t2_, p2_)));
  EXPECT_THAT(times, ElementsAre(t2_));
}

TEST_F(DiscreteTrajectoryDeathTest, AppendError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t2_, d2_);
    massive_trajectory_->Append(t1_, d1_);
  }, "out of order");
}

TEST_F(DiscreteTrajectoryTest, AppendAtExistingTime) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t1_, d1_);
}

TEST_F(DiscreteTrajectoryTest, AppendSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  std::map<Instant, Position<World>> const positions =
      Positions(*massive_trajectory_);
  std::map<Instant, Velocity<World>> const velocities =
      Velocities(*massive_trajectory_);
  std::list<Instant> const times = Times(*massive_trajectory_);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));
}

TEST_F(DiscreteTrajectoryTest, ForgetAfter) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork->Append(t4_, d4_);

  fork->ForgetAfter(t3_ + (t4_ - t3_) / 2);
  std::map<Instant, Position<World>> positions = Positions(*fork);
  std::map<Instant, Velocity<World>> velocities = Velocities(*fork);
  std::list<Instant> times = Times(*fork);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));

  fork->ForgetAfter(t2_);
  positions = Positions(*fork);
  velocities = Velocities(*fork);
  times = Times(*fork);
  EXPECT_THAT(positions, ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_)));
  EXPECT_THAT(velocities, ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_));
  EXPECT_EQ(q2_, fork->back().degrees_of_freedom.position());
  EXPECT_EQ(p2_, fork->back().degrees_of_freedom.velocity());
  EXPECT_EQ(t2_, fork->back().time);

  positions = Positions(*massive_trajectory_);
  velocities = Velocities(*massive_trajectory_);
  times = Times(*massive_trajectory_);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t1_, q1_), Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t1_, p1_), Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t1_, t2_, t3_));

  massive_trajectory_->ForgetAfter(t1_);
  positions = Positions(*massive_trajectory_);
  velocities = Velocities(*massive_trajectory_);
  times = Times(*massive_trajectory_);
  EXPECT_THAT(positions, ElementsAre(Pair(t1_, q1_)));
  EXPECT_THAT(velocities, ElementsAre(Pair(t1_, p1_)));
  EXPECT_THAT(times, ElementsAre(t1_));
  // Don't use fork, it is dangling.
}

TEST_F(DiscreteTrajectoryDeathTest, ForgetBeforeError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    massive_trajectory_->Append(t2_, d2_);
    massive_trajectory_->Append(t3_, d3_);
    not_null<DiscreteTrajectory<World>*> const fork =
        massive_trajectory_->NewForkWithCopy(t2_);
    fork->Append(t4_, d4_);
    massive_trajectory_->ForgetBefore(t3_);
  }, "found 1 fork");
}

TEST_F(DiscreteTrajectoryTest, ForgetBeforeSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork->Append(t4_, d4_);

  massive_trajectory_->ForgetBefore(t1_ + (t2_ - t1_) / 2);
  std::map<Instant, Position<World>> positions =
      Positions(*massive_trajectory_);
  std::map<Instant, Velocity<World>> velocities =
      Velocities(*massive_trajectory_);
  std::list<Instant> times = Times(*massive_trajectory_);
  EXPECT_THAT(positions, ElementsAre(Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t2_, t3_));
  positions = Positions(*fork);
  velocities = Velocities(*fork);
  times = Times(*fork);
  EXPECT_THAT(positions,
              ElementsAre(Pair(t2_, q2_), Pair(t3_, q3_), Pair(t4_, q4_)));
  EXPECT_THAT(velocities,
              ElementsAre(Pair(t2_, p2_), Pair(t3_, p3_), Pair(t4_, p4_)));
  EXPECT_THAT(times, ElementsAre(t2_, t3_, t4_));

  massive_trajectory_->ForgetBefore(t2_);
  positions = Positions(*massive_trajectory_);
  velocities = Velocities(*massive_trajectory_);
  times = Times(*massive_trajectory_);
  EXPECT_THAT(positions, ElementsAre(Pair(t2_, q2_), Pair(t3_, q3_)));
  EXPECT_THAT(velocities, ElementsAre(Pair(t2_, p2_), Pair(t3_, p3_)));
  EXPECT_THAT(times, ElementsAre(t2_, t3_));
}

TEST_F(DiscreteTrajectoryTest, ForgetBeforeEmpty) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->ForgetBefore(t4_);
  EXPECT_TRUE(massive_trajectory_->Empty());
}

TEST_F(DiscreteTrajectoryDeathTest, SerializationError) {
  EXPECT_DEATH({
    massive_trajectory_->Append(t1_, d1_);
    not_null<DiscreteTrajectory<World>*> const fork =
        massive_trajectory_->NewForkWithCopy(t1_);
    serialization::DiscreteTrajectory message;
        fork->WriteToMessage(/*excluded=*/{}, /*tracked=*/{}, /*exact=*/{},
                             &message);
  }, "is_root");
}

TEST_F(DiscreteTrajectoryTest, SerializationSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork0 =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork0->Append(t4_, d4_);
  not_null<DiscreteTrajectory<World>*> const fork1 =
      massive_trajectory_->NewForkWithCopy(t2_);
  not_null<DiscreteTrajectory<World>*> const fork2 =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork2->Append(t4_, d4_);
  not_null<DiscreteTrajectory<World>*> const fork3 =
      massive_trajectory_->NewForkWithCopy(t3_);
  fork3->Append(t4_, d4_);
  serialization::DiscreteTrajectory message;
  serialization::DiscreteTrajectory reference_message;

  massive_trajectory_->WriteToMessage(
      /*excluded=*/{fork0},
      /*tracked=*/{massive_trajectory_.get(), fork1, fork3, fork2},
      /*exact=*/{},
      &message);
  massive_trajectory_->WriteToMessage(
      /*excluded=*/{fork0},
      /*tracked=*/{massive_trajectory_.get(), fork1, fork3, fork2},
      /*exact=*/{},
      &reference_message);

  DiscreteTrajectory<World>* deserialized_root = nullptr;
  DiscreteTrajectory<World>* deserialized_fork1 = nullptr;
  DiscreteTrajectory<World>* deserialized_fork2 = nullptr;
  DiscreteTrajectory<World>* deserialized_fork3 = nullptr;
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const
      deserialized_trajectory =
          DiscreteTrajectory<World>::ReadFromMessage(message,
                                                     {&deserialized_root,
                                                      &deserialized_fork1,
                                                      &deserialized_fork3,
                                                      &deserialized_fork2});
  EXPECT_EQ(deserialized_trajectory.get(), deserialized_root);
  EXPECT_EQ(t2_, deserialized_fork1->Fork()->time);
  EXPECT_EQ(t2_, deserialized_fork2->Fork()->time);
  EXPECT_EQ(t3_, deserialized_fork3->Fork()->time);
  message.Clear();
  deserialized_trajectory->WriteToMessage(/*excluded=*/{},
                                          /*tracked=*/{deserialized_root,
                                                       deserialized_fork1,
                                                       deserialized_fork3,
                                                       deserialized_fork2},
                                          /*exact=*/{},
                                          &message);
  EXPECT_THAT(reference_message, EqualsProto(message));
  EXPECT_THAT(message.children_size(), Eq(2));
  EXPECT_THAT(message.zfp().timeline_size(), Eq(3));
  EXPECT_THAT(message.children(0).trajectories_size(), Eq(2));
  EXPECT_THAT(message.children(0).trajectories(0).children_size(), Eq(0));
  EXPECT_THAT(message.children(0).trajectories(0).zfp().timeline_size(), Eq(1));
  EXPECT_THAT(message.children(0).trajectories(1).children_size(), Eq(0));
  EXPECT_THAT(message.children(0).trajectories(1).zfp().timeline_size(), Eq(2));
  EXPECT_THAT(message.children(1).trajectories_size(), Eq(1));
  EXPECT_THAT(message.children(1).trajectories(0).children_size(), Eq(0));
  EXPECT_THAT(message.children(1).trajectories(0).zfp().timeline_size(), Eq(1));
}

TEST_F(DiscreteTrajectoryTest, SerializationAfterTime) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  not_null<DiscreteTrajectory<World>*> const fork0 =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork0->Append(t4_, d4_);
  not_null<DiscreteTrajectory<World>*> const fork1 =
      massive_trajectory_->NewForkWithCopy(t2_);
  not_null<DiscreteTrajectory<World>*> const fork2 =
      massive_trajectory_->NewForkWithCopy(t2_);
  fork2->Append(t4_, d4_);
  not_null<DiscreteTrajectory<World>*> const fork3 =
      massive_trajectory_->NewForkWithCopy(t3_);
  fork3->Append(t4_, d4_);

  {
    serialization::DiscreteTrajectory message;
    DiscreteTrajectory<World>* deserialized_root = nullptr;
    DiscreteTrajectory<World>* deserialized_fork1 = nullptr;
    DiscreteTrajectory<World>* deserialized_fork2 = nullptr;
    DiscreteTrajectory<World>* deserialized_fork3 = nullptr;
    massive_trajectory_->WriteToMessage(
        /*after_time=*/t1_,
        /*excluded=*/{fork0},
        /*tracked=*/{massive_trajectory_.get(), fork1, fork2, fork3},
        /*exact=*/{},
        &message);

    not_null<std::unique_ptr<DiscreteTrajectory<World>>> const
        deserialized_trajectory =
            DiscreteTrajectory<World>::ReadFromMessage(message,
                                                       {&deserialized_root,
                                                        &deserialized_fork1,
                                                        &deserialized_fork2,
                                                        &deserialized_fork3});
    EXPECT_EQ(t1_, deserialized_trajectory->begin()->time);
    EXPECT_EQ(deserialized_trajectory.get(), deserialized_root);
    EXPECT_EQ(t2_, deserialized_fork1->Fork()->time);
    EXPECT_EQ(t2_, deserialized_fork2->Fork()->time);
    EXPECT_EQ(t3_, deserialized_fork3->Fork()->time);
  }
  {
    serialization::DiscreteTrajectory message;
    DiscreteTrajectory<World>* deserialized_root = nullptr;
    DiscreteTrajectory<World>* deserialized_fork1 = nullptr;
    DiscreteTrajectory<World>* deserialized_fork2 = nullptr;
    DiscreteTrajectory<World>* deserialized_fork3 = nullptr;
    massive_trajectory_->WriteToMessage(
        /*after_time=*/t2_,
        /*excluded=*/{fork0},
        /*tracked=*/{massive_trajectory_.get(), fork1, fork2, fork3},
        /*exact=*/{},
        &message);

    not_null<std::unique_ptr<DiscreteTrajectory<World>>> const
        deserialized_trajectory =
            DiscreteTrajectory<World>::ReadFromMessage(message,
                                                       {&deserialized_root,
                                                        &deserialized_fork1,
                                                        &deserialized_fork2,
                                                        &deserialized_fork3});
    EXPECT_EQ(t2_, deserialized_trajectory->begin()->time);
    EXPECT_EQ(deserialized_trajectory.get(), deserialized_root);
    EXPECT_EQ(t2_, deserialized_fork1->Fork()->time);
    EXPECT_EQ(t2_, deserialized_fork2->Fork()->time);
    EXPECT_EQ(t3_, deserialized_fork3->Fork()->time);
  }
  {
    serialization::DiscreteTrajectory message;
    DiscreteTrajectory<World>* deserialized_root = nullptr;
    DiscreteTrajectory<World>* deserialized_fork3 = nullptr;
    massive_trajectory_->WriteToMessage(
        /*after_time=*/t3_,
        /*excluded=*/{},
        /*tracked=*/{massive_trajectory_.get(), fork3},
        /*exact=*/{},
        &message);

    not_null<std::unique_ptr<DiscreteTrajectory<World>>> const
        deserialized_trajectory =
            DiscreteTrajectory<World>::ReadFromMessage(message,
                                                       {&deserialized_root,
                                                        &deserialized_fork3});
    EXPECT_EQ(t3_, deserialized_trajectory->begin()->time);
    EXPECT_EQ(deserialized_trajectory.get(), deserialized_root);
    EXPECT_EQ(t3_, deserialized_fork3->Fork()->time);
  }
  {
    serialization::DiscreteTrajectory message;
    DiscreteTrajectory<World>* deserialized_root = nullptr;
    massive_trajectory_->WriteToMessage(
        /*after_time=*/t4_,
        /*excluded=*/{},
        /*tracked=*/{massive_trajectory_.get()},
        /*exact=*/{},
        &message);

    not_null<std::unique_ptr<DiscreteTrajectory<World>>> const
        deserialized_trajectory =
            DiscreteTrajectory<World>::ReadFromMessage(message,
                                                       {&deserialized_root});
    EXPECT_TRUE(deserialized_trajectory->Empty());
    EXPECT_EQ(deserialized_trajectory.get(), deserialized_root);
  }
}

TEST_F(DiscreteTrajectoryDeathTest, LastError) {
  EXPECT_DEATH({
    massive_trajectory_->back();
  }, "parent_.*non NULL");
}

TEST_F(DiscreteTrajectoryTest, LastSuccess) {
  massive_trajectory_->Append(t1_, d1_);
  massive_trajectory_->Append(t2_, d2_);
  massive_trajectory_->Append(t3_, d3_);
  EXPECT_EQ(q3_, massive_trajectory_->back().degrees_of_freedom.position());
  EXPECT_EQ(p3_, massive_trajectory_->back().degrees_of_freedom.velocity());
  EXPECT_EQ(t3_, massive_trajectory_->back().time);
}

TEST_F(DiscreteTrajectoryDeathTest, IteratorError) {
  EXPECT_DEATH({
    DiscreteTrajectory<World>::Iterator it = --massive_trajectory_->end();
  }, "parent_.*non NULL");
}

TEST_F(DiscreteTrajectoryTest, IteratorSuccess) {
  DiscreteTrajectory<World>::Iterator it = massive_trajectory_->begin();
  EXPECT_TRUE(it == massive_trajectory_->end());

  massless_trajectory_->Append(t1_, d1_);
  massless_trajectory_->Append(t2_, d2_);
  massless_trajectory_->Append(t3_, d3_);

  it = massless_trajectory_->begin();
  EXPECT_FALSE(it == massless_trajectory_->end());
  EXPECT_EQ(t1_, it->time);
  EXPECT_EQ(d1_, it->degrees_of_freedom);
  ++it;
  EXPECT_EQ(t2_, it->time);
  EXPECT_EQ(d2_, it->degrees_of_freedom);
  ++it;
  EXPECT_EQ(t3_, it->time);
  EXPECT_EQ(d3_, it->degrees_of_freedom);
  ++it;
  EXPECT_TRUE(it == massless_trajectory_->end());

  not_null<DiscreteTrajectory<World>*> const fork =
      massless_trajectory_->NewForkWithCopy(t2_);
  fork->Append(t4_, d4_);

  it = fork->begin();
  EXPECT_FALSE(it == fork->end());
  EXPECT_EQ(t1_, it->time);
  EXPECT_EQ(d1_, it->degrees_of_freedom);
  ++it;
  EXPECT_EQ(t2_, it->time);
  EXPECT_EQ(d2_, it->degrees_of_freedom);
  ++it;
  EXPECT_EQ(t3_, it->time);
  EXPECT_EQ(d3_, it->degrees_of_freedom);
  ++it;
  EXPECT_EQ(t4_, it->time);
  EXPECT_EQ(d4_, it->degrees_of_freedom);
  ++it;
  EXPECT_TRUE(it == fork->end());
}

TEST_F(DiscreteTrajectoryTest, QuadrilateralCircle) {
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Speed const v = ω * r / Radian;
  Time const period = 2 * π * Radian / ω;
  auto const circle = NewCircularTrajectory<World>(
      ω, r, /*Δt=*/period / 4, /*t1=*/t0_, /*t2=*/t0_ + period);

  double max_r_error = 0;
  double max_v_error = 0;
  for (Instant t = t0_; t <= circle->back().time; t += period / 32) {
    auto const degrees_of_freedom_interpolated =
        circle->EvaluateDegreesOfFreedom(t);
    auto const& q_interpolated = degrees_of_freedom_interpolated.position();
    auto const& v_interpolated = degrees_of_freedom_interpolated.velocity();
    EXPECT_THAT(circle->EvaluatePosition(t), Eq(q_interpolated));
    EXPECT_THAT(circle->EvaluateVelocity(t), Eq(v_interpolated));
    max_r_error = std::max(
        max_r_error, RelativeError(r, (q_interpolated - World::origin).Norm()));
    max_v_error =
        std::max(max_v_error, RelativeError(v, v_interpolated.Norm()));
  }
  EXPECT_THAT(max_r_error, IsNear(0.015_⑴));
  EXPECT_THAT(max_v_error, IsNear(0.012_⑴));
}

TEST_F(DiscreteTrajectoryTest, Downsampling) {
  DiscreteTrajectory<World> circle;
  DiscreteTrajectory<World> downsampled_circle;
  downsampled_circle.SetDownsampling({.max_dense_intervals = 50,
                                      .tolerance = 1 * Milli(Metre)});
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Time const Δt = 10 * Milli(Second);
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 10 * Second;
  AppendTrajectory(*NewCircularTrajectory<World>(ω, r, Δt, t1, t2),
                   /*to=*/circle);
  AppendTrajectory(*NewCircularTrajectory<World>(ω, r, Δt, t1, t2),
                   /*to=*/downsampled_circle);

  EXPECT_THAT(circle.Size(), Eq(1001));
  EXPECT_THAT(downsampled_circle.Size(), Eq(77));
  std::vector<Length> errors;
  for (auto const& [time, degrees_of_freedom] : circle) {
    errors.push_back((downsampled_circle.EvaluatePosition(time) -
                      degrees_of_freedom.position()).Norm());
  }
  EXPECT_THAT(errors, Each(Lt(1 * Milli(Metre))));
  EXPECT_THAT(errors, Contains(Gt(9 * Micro(Metre))))
      << *std::max_element(errors.begin(), errors.end());
}

TEST_F(DiscreteTrajectoryTest, DownsamplingForks) {
  // The number of dense intervals does not divide the lengths of the
  // trajectories.
  int const max_dense_intervals = 73;
  Length const tolerance = 1 * Milli(Metre);

  DiscreteTrajectory<World> root;
  root.SetDownsampling({max_dense_intervals, tolerance});
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Time const Δt = 10 * Milli(Second);
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 10 * Second;
  AppendTrajectory(*NewCircularTrajectory<World>(ω, r, Δt, t1, t2),
                   /*to=*/root);
  auto const root_tmax = root.back().time;

  auto fork1 = root.NewForkAtLast();
  fork1->SetDownsampling({max_dense_intervals, tolerance});
  Instant const t3 = t2 + 10 * Second;
  AppendTrajectory(*NewCircularTrajectory<World>(ω, r, Δt, root_tmax + Δt, t3),
                   /*to=*/*fork1);
  auto const fork1_tmax = fork1->back().time;

  // A short fork with no downsampling
  auto fork2 = fork1->NewForkAtLast();
  fork2->SetDownsampling({max_dense_intervals, tolerance});
  Instant const t4 = t3 + 300 * Milli(Second);
  AppendTrajectory(*NewCircularTrajectory<World>(ω, r, Δt, fork1_tmax + Δt, t4),
                   /*to=*/*fork2);
  auto const fork2_tmax = fork2->back().time;

  auto fork3 = fork2->NewForkAtLast();
  fork3->SetDownsampling({max_dense_intervals, tolerance});
  Instant const t5 = t4 + 10 * Second;
  AppendTrajectory(*NewCircularTrajectory<World>(ω, r, Δt, fork2_tmax + Δt, t5),
                   /*to=*/*fork3);

  // Roughly 55 points per downsampled segments.
  EXPECT_THAT(root.Size(), Eq(56));
  EXPECT_THAT(fork1->Size(), Eq(root.Size() + 54));
  EXPECT_THAT(fork2->Size(), Eq(fork1->Size() + 30));
  EXPECT_THAT(fork3->Size(), Eq(fork2->Size() + 55));

  for (Instant t = t0_; t <= fork3->t_max(); t += Δt) {
    EXPECT_THAT((fork3->EvaluatePosition(t) - World::origin).Norm(),
                AbsoluteErrorFrom(r, Lt(0.99 * Milli(Metre))));
  }
}

TEST_F(DiscreteTrajectoryTest, DownsamplingSerialization) {
  DiscreteTrajectory<World> circle;
  circle.SetDownsampling({.max_dense_intervals = 50,
                          .tolerance = 1 * Milli(Metre)});
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Time const Δt = 10 * Milli(Second);
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 5 * Second;
  AppendTrajectory(*NewCircularTrajectory<World>(ω, r, Δt, t1, t2),
                   /*to=*/circle);
  auto const circle_tmax = circle.back().time;

  serialization::DiscreteTrajectory message;
  circle.WriteToMessage(/*excluded=*/{},
                        /*tracked=*/{},
                        /*exact=*/{circle.LowerBound(t0_ + 2 * Second),
                                   circle.LowerBound(t0_ + 3 * Second)},
                        &message);
  auto deserialized_circle =
      DiscreteTrajectory<World>::ReadFromMessage(message, /*forks=*/{});

  // Serialization/deserialization preserves the size, the times, and nudges the
  // positions by less than the tolerance.  It also preserve the degrees of
  // freedom at the "exact" iterators.
  EXPECT_THAT(circle.Size(), Eq(39));
  EXPECT_THAT(deserialized_circle->Size(), circle.Size());
  for (auto it1 = circle.begin(), it2 = deserialized_circle->begin();
       it1 != circle.end();
       ++it1, ++it2) {
    EXPECT_EQ(it1->time, it2->time);
    EXPECT_THAT(it2->degrees_of_freedom.position(),
                AbsoluteErrorFrom(it1->degrees_of_freedom.position(),
                                  Lt(0.22 * Milli(Metre))));
    EXPECT_THAT(it2->degrees_of_freedom.velocity(),
                AbsoluteErrorFrom(it1->degrees_of_freedom.velocity(),
                                  Lt(1.1 * Milli(Metre) / Second)));
  }
  EXPECT_NE(
      deserialized_circle->LowerBound(t0_ + 1 * Second)->degrees_of_freedom,
      circle.LowerBound(t0_ + 1 * Second)->degrees_of_freedom);
  EXPECT_EQ(
      deserialized_circle->LowerBound(t0_ + 2 * Second)->degrees_of_freedom,
      circle.LowerBound(t0_ + 2 * Second)->degrees_of_freedom);
  EXPECT_EQ(
      deserialized_circle->LowerBound(t0_ + 3 * Second)->degrees_of_freedom,
      circle.LowerBound(t0_ + 3 * Second)->degrees_of_freedom);
  EXPECT_NE(
      deserialized_circle->LowerBound(t0_ + 4 * Second)->degrees_of_freedom,
      circle.LowerBound(t0_ + 4 * Second)->degrees_of_freedom);

  // Appending may result in different downsampling because the positions differ
  // a bit.
  Instant const t3 = t0_ + 10 * Second;
  AppendTrajectory(
      *NewCircularTrajectory<World>(ω, r, Δt, circle_tmax + Δt, t3),
      /*to=*/circle);
  AppendTrajectory(
      *NewCircularTrajectory<World>(ω, r, Δt, circle_tmax + Δt, t3),
      /*to=*/*deserialized_circle);

  // Despite the difference in downsampling (and therefore in size) the two
  // trajectories are still within the tolerance.
  EXPECT_THAT(circle.Size(), Eq(77));
  EXPECT_THAT(deserialized_circle->Size(), Eq(78));
  for (Instant t = t0_;
       t <= std::min(circle.back().time, deserialized_circle->back().time);
       t += Δt) {
    EXPECT_THAT(
        deserialized_circle->EvaluatePosition(t),
        AbsoluteErrorFrom(circle.EvaluatePosition(t), Le(0.23 * Milli(Metre))));
    EXPECT_THAT(deserialized_circle->EvaluateVelocity(t),
                AbsoluteErrorFrom(circle.EvaluateVelocity(t),
                                  Le(5.7 * Milli(Metre) / Second)));
  }
}

TEST_F(DiscreteTrajectoryTest, DownsamplingForgetAfter) {
  DiscreteTrajectory<World> circle;
  DiscreteTrajectory<World> forgotten_circle;
  circle.SetDownsampling({.max_dense_intervals = 50,
                          .tolerance = 1 * Milli(Metre)});
  forgotten_circle.SetDownsampling({.max_dense_intervals = 50,
                                    .tolerance = 1 * Milli(Metre)});
  AngularFrequency const ω = 3 * Radian / Second;
  Length const r = 2 * Metre;
  Time const Δt = 1.0 / 128.0 * Second;  // Yields exact times.
  Instant const t1 = t0_;
  Instant const t2 = t0_ + 10 * Second;
  AppendTrajectory(*NewCircularTrajectory<World>(ω, r, Δt, t1, t2),
                   /*to=*/circle);
  AppendTrajectory(*NewCircularTrajectory<World>(ω, r, Δt, t1, t2),
                   /*to=*/forgotten_circle);

  forgotten_circle.ForgetAfter(t2);
  AppendTrajectory(*NewCircularTrajectory<World>(
                       ω, r, Δt, forgotten_circle.back().time + Δt, t2),
                   /*to=*/forgotten_circle);
  EXPECT_THAT(circle.Size(), Eq(92));
  EXPECT_THAT(forgotten_circle.Size(), Eq(circle.Size()));
  std::vector<Length> errors;
  for (auto const [time, degrees_of_freedom] : forgotten_circle) {
    errors.push_back((circle.Find(time)->degrees_of_freedom.position() -
                      degrees_of_freedom.position()).Norm());
  }
  EXPECT_THAT(errors, Each(Eq(0 * Metre)));
}

}  // namespace internal_discrete_trajectory
}  // namespace physics
}  // namespace principia
