#include "geometry/rotation.hpp"

#include "base/heap_checker.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace geometry {

using si::Metre;
using testing::Eq;
using testing_utilities::AlmostEquals;

class RotationTest : public testing::Test {
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
    e1_ = Vector<double, World>(R3Element<double>({1, 0, 0}));
    e2_ = Vector<double, World>(R3Element<double>({0, 1, 0}));
    e3_ = Vector<double, World>(R3Element<double>({0, 0, 1}));
    rotation_a_ = Rot(120 * si::Degree,
                      Bivector<double, World>({1, 1, 1}));
    rotation_b_ = Rot(90 * si::Degree,
                      Bivector<double, World>({1, 0, 0}));
    rotation_c_ = Rot(R3x3Matrix({{0.5, 0.5 * sqrt(3), 0},
                                  {-0.5 * sqrt(3), 0.5, 0},
                                  {0, 0, 1}}));
  }

  Vector<quantities::Length, World> vector_;
  Bivector<quantities::Length, World> bivector_;
  Trivector<quantities::Length, World> trivector_;
  Vector<double, World> e1_;
  Vector<double, World> e2_;
  Vector<double, World> e3_;
  Rot rotation_a_;
  Rot rotation_b_;
  Rot rotation_c_;

 private:
  base::HeapChecker heap_checker_;
};

TEST_F(RotationTest, Identity) {
  EXPECT_THAT(vector_, Eq(Rot::Identity()(vector_)));
  EXPECT_THAT(bivector_, Eq(Rot::Identity()(bivector_)));
  EXPECT_THAT(trivector_, Eq(Rot::Identity()(trivector_)));
}

TEST_F(RotationTest, AppliedToVector) {
  EXPECT_THAT(rotation_a_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre)), 4));
  EXPECT_THAT(rotation_b_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                -3.0 * Metre,
                                                2.0 * Metre)), 1));
  EXPECT_THAT(rotation_c_(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>((0.5 + sqrt(3.0)) * Metre,
                                                (1.0 - 0.5 * sqrt(3.0)) * Metre,
                                                3.0 * Metre)), 0));
}

TEST_F(RotationTest, AppliedToBivector) {
  EXPECT_THAT(rotation_a_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre)), 4));
  EXPECT_THAT(rotation_b_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                -3.0 * Metre,
                                                2.0 * Metre)), 1));
  EXPECT_THAT(rotation_c_(bivector_),
              AlmostEquals(Bivector<quantities::Length, World>(
                  R3Element<quantities::Length>((0.5 + sqrt(3.0)) * Metre,
                                                (1.0 - 0.5 * sqrt(3.0)) * Metre,
                                                3.0 * Metre)), 0));
}

TEST_F(RotationTest, AppliedToTrivector) {
  EXPECT_THAT(rotation_a_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre), 0));
  EXPECT_THAT(rotation_b_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre), 0));
  EXPECT_THAT(rotation_c_(trivector_),
              AlmostEquals(Trivector<quantities::Length, World>(
                  4.0 * Metre), 0));
}

TEST_F(RotationTest, Determinant) {
  EXPECT_TRUE(rotation_a_.Determinant().Positive());
  EXPECT_TRUE(rotation_b_.Determinant().Positive());
  EXPECT_TRUE(rotation_c_.Determinant().Positive());
}

TEST_F(RotationTest, Inverse) {
  EXPECT_THAT(rotation_a_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(2.0 * Metre,
                                                3.0 * Metre,
                                                1.0 * Metre)), 2));
  EXPECT_THAT(rotation_b_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(1.0 * Metre,
                                                3.0 * Metre,
                                                -2.0 * Metre)), 1));
  EXPECT_THAT(rotation_c_.Inverse()(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>((0.5 - sqrt(3.0)) * Metre,
                                                (1.0 + 0.5 * sqrt(3.0)) * Metre,
                                                3.0 * Metre)), 0));
}

TEST_F(RotationTest, Composition) {
  Rot const rotation_ab = rotation_a_ * rotation_b_;
  EXPECT_THAT(rotation_ab(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(2.0 * Metre,
                                                1.0 * Metre,
                                                -3.0 * Metre)), 4));
}

TEST_F(RotationTest, Forget) {
  Orth const orthogonal_a = rotation_a_.Forget();
  EXPECT_THAT(orthogonal_a(vector_),
              AlmostEquals(Vector<quantities::Length, World>(
                  R3Element<quantities::Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre)), 4));
}

// These four tests cover all the branches of ToQuaternion.
TEST_F(RotationTest, ToQuaternion1) {
  R3Element<double> const v1 = {2, 5, 6};
  R3Element<double> v2 = {-3, 4, 1};
  v1.Orthogonalize(&v2);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix m = {w1, w2, w3};
  Rot rotation(m.Transpose());
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 6));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 5));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 1));
}

TEST_F(RotationTest, ToQuaternion2) {
  R3Element<double> const v1 = {-2, -5, -6};
  R3Element<double> v2 = {-3, 4, 1};
  v1.Orthogonalize(&v2);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix m = {w1, w2, w3};
  Rot rotation(m.Transpose());
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 6));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 5));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 1));
}

TEST_F(RotationTest, ToQuaternion3) {
  R3Element<double> const v1 = {-2, -5, -6};
  R3Element<double> v2 = {-3, -4, 1};
  v1.Orthogonalize(&v2);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix m = {w1, w2, w3};
  Rot rotation(m.Transpose());
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 2));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 1));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 12));
}

TEST_F(RotationTest, ToQuaternion4) {
  R3Element<double> const v1 = {-2, -5, -6};
  R3Element<double> v2 = {-3, -4, -1};
  v1.Orthogonalize(&v2);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix m = {w1, w2, w3};
  Rot rotation(m.Transpose());
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 6));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 1));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 2));
}

}  // namespace geometry
}  // namespace principia
