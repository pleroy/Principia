
#include "numerics/quadrature.hpp"

#include <limits>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace numerics {
namespace quadrature {

using geometry::Displacement;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using quantities::Angle;
using quantities::Cos;
using quantities::Pow;
using quantities::Sin;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::operator""_⑴;
using ::testing::AnyOf;
using ::testing::Eq;

class QuadratureTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;
};

TEST_F(QuadratureTest, Sin) {
  int evaluations = 0;
  auto const f = [&evaluations](Angle const& x) {
    ++evaluations;
    return Sin(x);
  };
  auto const ʃf = (Cos(2.0 * Radian) - Cos(5.0 * Radian)) * Radian;
  EXPECT_THAT(GaussLegendre<1>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(10_⑴)));
  EXPECT_THAT(GaussLegendre<2>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(3.3_⑴)));
  EXPECT_THAT(GaussLegendre<3>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(4.0e-1_⑴)));
  EXPECT_THAT(GaussLegendre<4>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(2.4e-2_⑴)));
  EXPECT_THAT(GaussLegendre<5>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(8.6e-4_⑴)));
  EXPECT_THAT(GaussLegendre<6>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(2.1e-5_⑴)));
  EXPECT_THAT(GaussLegendre<7>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(3.6e-7_⑴)));
  EXPECT_THAT(GaussLegendre<8>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(4.7e-9_⑴)));
  EXPECT_THAT(GaussLegendre<9>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(4.8e-11_⑴)));
  EXPECT_THAT(GaussLegendre<10>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 2495, 2498));
  EXPECT_THAT(GaussLegendre<11>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 20));
  EXPECT_THAT(GaussLegendre<12>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 6, 7));
  EXPECT_THAT(GaussLegendre<13>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 0, 1));
  EXPECT_THAT(GaussLegendre<14>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 1, 2));
  EXPECT_THAT(GaussLegendre<15>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 3, 4));

  evaluations = 0;
  EXPECT_THAT(AutomaticClenshawCurtis(
                  f,
                  -2.0 * Radian,
                  5.0 * Radian,
                  /*max_relative_error=*/std::numeric_limits<double>::epsilon(),
                  /*max_points=*/std::nullopt),
              AlmostEquals(ʃf, 3, 4));
  EXPECT_THAT(evaluations, Eq(65));
}

TEST_F(QuadratureTest, Sin2) {
  auto const f = [](Angle const& x) { return Sin(2 * x); };
  auto const ʃf = (Cos(4 * Radian) - Cos(10 * Radian)) / 2 * Radian;
  EXPECT_THAT(GaussLegendre<1>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(9.7_⑴)));
  EXPECT_THAT(GaussLegendre<2>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(7.6_⑴)));
  EXPECT_THAT(GaussLegendre<3>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(7.6_⑴)));
  EXPECT_THAT(GaussLegendre<4>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(2.4_⑴)));
  EXPECT_THAT(GaussLegendre<5>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(4.2e-1_⑴)));
  EXPECT_THAT(GaussLegendre<6>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(4.6e-2_⑴)));
  EXPECT_THAT(GaussLegendre<7>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(3.5e-3_⑴)));
  EXPECT_THAT(GaussLegendre<8>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(2.0e-4_⑴)));
  EXPECT_THAT(GaussLegendre<9>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(8.5e-6_⑴)));
  EXPECT_THAT(GaussLegendre<10>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(2.9e-7_⑴)));;
  EXPECT_THAT(GaussLegendre<11>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(8.1e-9_⑴)));
  EXPECT_THAT(GaussLegendre<12>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(1.9e-10_⑴)));
  EXPECT_THAT(GaussLegendre<13>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(3.7e-12_⑴)));
  EXPECT_THAT(GaussLegendre<14>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 422, 425));
  EXPECT_THAT(GaussLegendre<15>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 36, 50));
  EXPECT_THAT(GaussLegendre<16>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 152, 159));
}

TEST_F(QuadratureTest, Sin10) {
  int evaluations = 0;
  auto const f = [&evaluations](Angle const& x) {
    ++evaluations;
    return Sin(10 * x);
  };
  auto const ʃf = (Cos(20 * Radian) - Cos(50 * Radian)) / 10 * Radian;
  EXPECT_THAT(AutomaticClenshawCurtis(
                  f,
                  -2.0 * Radian,
                  5.0 * Radian,
                  /*max_relative_error=*/std::numeric_limits<double>::epsilon(),
                  /*max_points=*/std::nullopt),
              AlmostEquals(ʃf, 2, 20));
  EXPECT_THAT(evaluations, AnyOf(Eq(32769), Eq(65537), Eq(262145), Eq(524289)));
}

#if 0
TEST_F(QuadratureTest, NonTerminating) {
  int evaluations = 0;
  auto const left =
      PolynomialInMonomialBasis<Displacement<World>, Time, 3, EstrinEvaluator>(
          {Displacement<World>({-2.64575131106458450e+00 * Metre,
                                +0.00000000000000000e+00 * Metre,
                                +0.00000000000000000e+00 * Metre}),
           Displacement<World>({+4.02334445113228672e-06 * Metre,
                                +0.00000000000000000e+00 * Metre,
                                +0.00000000000000000e+00 * Metre}) /
               Second,
           Displacement<World>({-1.27463011048138714e-12 * Metre,
                                +0.00000000000000000e+00 * Metre,
                                +0.00000000000000000e+00 * Metre}) /
               Pow<2>(Second),
           Displacement<World>({+1.07683673837640856e-19 * Metre,
                                +0.00000000000000000e+00 * Metre,
                                +0.00000000000000000e+00 * Metre}) /
               Pow<3>(Second)});
  auto const right =
      PolynomialInMonomialBasis<Displacement<World>, Time, 2, EstrinEvaluator>(
          {Displacement<World>({-3.40821468615081086e-02 * Metre,
                                +0.00000000000000000e+00 * Metre,
                                +0.00000000000000000e+00 * Metre}),
           Displacement<World>({+2.05092765462352891e-08 * Metre,
                                +0.00000000000000000e+00 * Metre,
                                +0.00000000000000000e+00 * Metre}) /
               Second,
           Displacement<World>({-2.50948167739304229e-15 * Metre,
                                +0.00000000000000000e+00 * Metre,
                                +0.00000000000000000e+00 * Metre}) /
               Pow<2>(Second)});

  auto const f = [&evaluations, &left, &right](Time const& t) {
    ++evaluations;
    return InnerProduct(left(t), right(t));
  };
  EXPECT_THAT(AutomaticClenshawCurtis(
                  f,
                  0 * Second,
                  +7.89120000000000000e+06 * Second,
                  /*max_relative_error=*/0x1p-32,
                  /*max_points=*/std::nullopt),
              AlmostEquals(0 * Pow<2>(Metre) * Second, 0));
  EXPECT_THAT(evaluations, AnyOf(Eq(666)));
}
#endif

}  // namespace quadrature
}  // namespace numerics
}  // namespace principia
