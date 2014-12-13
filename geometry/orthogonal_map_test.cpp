#include "geometry/orthogonal_map.hpp"

#include "geometry/grassmann.hpp"
#include "geometry/rotation.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/heap_checked_test.hpp"

namespace principia {
namespace geometry {

using si::Metre;
using testing::Eq;
using testing_utilities::AlmostEquals;

class OrthogonalMapTest : public testing_utilities::HeapCheckedTest {
 protected:
  struct World;
  using Orth = OrthogonalMap<World, World>;
  using Rot = Rotation<World, World>;

  void SetUp() override {
    vector_ = Vector<quantities::Length, World>(
        R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
    bivector_ = Bivector<quantities::Length, World>(
        R3Element<quantities::Length>(1.0 * Metre, 2.0 * Metre, 3.0 * Metre));
    trivector_ = Trivector<quantities::Length, World>(4.0 * Metre);
    orthogonal_a_ = Orth(Sign(-1),
                         Rot(120 * si::Degree,
                             Bivector<double, World>({1, 1, 1})));
    orthogonal_b_ = Orth(Sign(1),
                         Rot(90 * si::Degree,
                             Bivector<double, World>({1, 0, 0})));
    orthogonal_c_ = Orth(Sign(-1),
                         Rot(90 * si::Degree,
                             Bivector<double, World>({1, 0, 0})));
  }

  Vector<quantities::Length, World> vector_;
  Bivector<quantities::Length, World> bivector_;
  Trivector<quantities::Length, World> trivector_;
  Orth orthogonal_a_;
  Orth orthogonal_b_;
  Orth orthogonal_c_;
};

TEST_F(OrthogonalMapTest, Identity) {
  EXPECT_THAT(vector_, Eq(Orth::Identity()(vector_)));
  EXPECT_THAT(bivector_, Eq(Orth::Identity()(bivector_)));
  EXPECT_THAT(trivector_, Eq(Orth::Identity()(trivector_)));
}

TEST_F(OrthogonalMapTest, AppliedToVector) {
  EXPECT_THAT(orthogonal_a_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(-3.0 * Metre,
                                                -1.0 * Metre,
                                                -2.0 * Metre)), 4));
  EXPECT_THAT(orthogonal_b_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                -3.0 * Metre,
                                                2.0 * Metre)), 1));
}

TEST_F(OrthogonalMapTest, AppliedToBivector) {
  EXPECT_THAT(orthogonal_a_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre)), 4));
  EXPECT_THAT(orthogonal_b_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                -3.0 * Metre,
                                                2.0 * Metre)), 1));
}

TEST_F(OrthogonalMapTest, AppliedToTrivector) {
  EXPECT_THAT(orthogonal_a_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  -4.0 * Metre), 0));
  EXPECT_THAT(orthogonal_b_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre), 0));
}

TEST_F(OrthogonalMapTest, Determinant) {
  EXPECT_TRUE(orthogonal_a_.Determinant().Negative());
  EXPECT_TRUE(orthogonal_b_.Determinant().Positive());
  EXPECT_TRUE(orthogonal_c_.Determinant().Negative());
}

TEST_F(OrthogonalMapTest, Inverse) {
  EXPECT_THAT(orthogonal_a_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(-2.0 * Metre,
                                                -3.0 * Metre,
                                                -1.0 * Metre)), 2));
  EXPECT_THAT(orthogonal_b_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                3.0 * Metre,
                                                -2.0 * Metre)), 1));
}

TEST_F(OrthogonalMapTest, Composition) {
  Orth const orthogonal_ac = orthogonal_a_ * orthogonal_c_;
  EXPECT_THAT(orthogonal_ac(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(2.0 * Metre,
                                                1.0 * Metre,
                                                -3.0 * Metre)), 4));
  EXPECT_TRUE((orthogonal_a_ * orthogonal_b_).Determinant().Negative());
  EXPECT_TRUE((orthogonal_a_ * orthogonal_c_).Determinant().Positive());
  EXPECT_TRUE((orthogonal_b_ * orthogonal_c_).Determinant().Negative());
}

}  // namespace geometry
}  // namespace principia
