
#include "numerics/polynomial.hpp"

#include <tuple>

#include "base/macros.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/numerics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/matchers.hpp"

#define PRINCIPIA_USE_IACA 0
#if PRINCIPIA_USE_IACA
#include "intel/iacaMarks.h"
#endif

namespace principia {

using geometry::Frame;
using geometry::Displacement;
using geometry::Handedness;
using geometry::Inertial;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using quantities::Acceleration;
using quantities::Energy;
using quantities::Entropy;
using quantities::Length;
using quantities::Product;
using quantities::Quotient;
using quantities::Current;
using quantities::Temperature;
using quantities::Time;
using quantities::constants::BoltzmannConstant;
using quantities::constants::SpeedOfLight;
using quantities::si::Ampere;
using quantities::si::Joule;
using quantities::si::Kelvin;
using quantities::si::Metre;
using quantities::si::Second;
using quantities::si::Watt;
using testing_utilities::AlmostEquals;
using testing_utilities::EqualsProto;
using ::testing::Eq;

namespace numerics {

class PolynomialTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using P2V = PolynomialInMonomialBasis<Displacement<World>, Time, 2,
                                        HornerEvaluator>;
  using P2A = PolynomialInMonomialBasis<Displacement<World>, Instant, 2,
                                        HornerEvaluator>;
  using P2P = PolynomialInMonomialBasis<Position<World>, Instant, 2,
                                        HornerEvaluator>;
  using P17 = PolynomialInMonomialBasis<Displacement<World>, Time, 17,
                                        EstrinEvaluator>;

  PolynomialTest()
      : coefficients_({
            Displacement<World>({0 * Metre,
                                 0 * Metre,
                                 1 * Metre}),
            Velocity<World>({0 * Metre / Second,
                             1 * Metre / Second,
                             0 * Metre / Second}),
            Vector<Acceleration, World>({1 * Metre / Second / Second,
                                         0 * Metre / Second / Second,
                                         0 * Metre / Second / Second})}) {}

  P2V::Coefficients const coefficients_;
};

#if PRINCIPIA_USE_IACA
// A convenient skeleton for analysing code with IACA.
TEST_F(PolynomialTest, DISABLED_IACA) {
  constexpr int degree = 17;
  using E = EstrinEvaluator<Displacement<World>, Time, degree>;
  using P = PolynomialInMonomialBasis<Displacement<World>,
                                      Time,
                                      degree,
                                      EstrinEvaluator>;
  P::Coefficients const coefficients;

  auto iaca = [](P::Coefficients const& c, Time const& t) {
    IACA_VC64_START;
    auto const result = E::Evaluate(c, t);
    IACA_VC64_END;
    return result;
  };
  CHECK_EQ(iaca(coefficients, 2 * Second), iaca(coefficients, 2 * Second));
}
#endif

// Check that coefficients can be accessed and have the right type.
TEST_F(PolynomialTest, Coefficients) {
  Displacement<World> const d = std::get<0>(coefficients_);
  Velocity<World> const v = std::get<1>(coefficients_);
  EXPECT_EQ(1 * Metre, d.coordinates().z);
  EXPECT_EQ(1 * Metre / Second, v.coordinates().y);
}

// Check that a polynomial can be constructed and evaluated.
TEST_F(PolynomialTest, Evaluate2V) {
  P2V const p(coefficients_);
  EXPECT_EQ(2, p.degree());
  Displacement<World> const d = p(0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0.25 * Metre,
                                                   0.5 * Metre,
                                                   1 * Metre}), 0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));
}

// Check that a polynomial can be for an affine argument.
TEST_F(PolynomialTest, Evaluate2A) {
  Instant const t0 = Instant() + 0.3 * Second;
  P2A const p(coefficients_, t0);
  EXPECT_EQ(2, p.degree());
  Displacement<World> const d = p(t0 + 0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(t0 + 0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0.25 * Metre,
                                                   0.5 * Metre,
                                                   1 * Metre}), 0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));

  // This compiles.
  p.Primitive();
}

// Check that a polynomial can return an affine value.
TEST_F(PolynomialTest, Evaluate2P) {
  Instant const t0 = Instant() + 0.3 * Second;
  P2P const p({World::origin + std::get<0>(coefficients_),
               std::get<1>(coefficients_),
               std::get<2>(coefficients_)},
              t0);
  EXPECT_EQ(2, p.degree());
  Position<World> const d = p(t0 + 0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(t0 + 0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(World::origin + Displacement<World>({0.25 * Metre,
                                                                   0.5 * Metre,
                                                                   1 * Metre}),
                              0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));

  // This doesn't compile (and rightly so).
#if 0
  p.Primitive();
#endif
}

// Check that a polynomial of high order may be declared.
TEST_F(PolynomialTest, Evaluate17) {
  P17::Coefficients const coefficients;
  P17 const p(coefficients);
  EXPECT_EQ(17, p.degree());
  Displacement<World> const d = p(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}), 0));
}

// Check that a conversion to increase the degree works.
TEST_F(PolynomialTest, Conversion) {
  P2V const p2v(coefficients_);
  P17 const p17 = P17(p2v);
  Displacement<World> const d = p17(0.5 * Second);
  Velocity<World> const v = p17.EvaluateDerivative(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0.25 * Metre,
                                                   0.5 * Metre,
                                                   1 * Metre}), 0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));
}

TEST_F(PolynomialTest, VectorSpace) {
  P2V const p2v(coefficients_);
  {
    auto const p = p2v + p2v;
    auto const actual = p(0 * Second);
    auto const expected =
        Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = p2v - p2v;
    auto const actual = p(0 * Second);
    auto const expected =
        Displacement<World>({0 * Metre, 0 * Metre, 0 * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = 3.0 * Joule * p2v;
    auto const actual = p(0 * Second);
    auto const expected = Vector<Product<Energy, Length>, World>(
                  {0 * Joule * Metre, 0 * Joule * Metre, 3 * Joule * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = p2v * (3.0 * Joule);
    auto const actual = p(0 * Second);
    auto const expected = Vector<Product<Length, Energy>, World>(
                  {0 * Joule * Metre, 0 * Joule * Metre, 3 * Joule * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = p2v / (4.0 * Joule);
    auto const actual = p(0 * Second);
    auto const expected = Vector<Quotient<Length, Energy>, World>(
                  {0 * Metre / Joule, 0 * Metre / Joule, 0.25 * Metre / Joule});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
}

TEST_F(PolynomialTest, Ring) {
  using P2 = PolynomialInMonomialBasis<Temperature, Time, 2, HornerEvaluator>;
  using P3 = PolynomialInMonomialBasis<Current, Time, 3, HornerEvaluator>;
  P2 const p2({1 * Kelvin, 3 * Kelvin / Second, -8 * Kelvin / Second / Second});
  P3 const p3({2 * Ampere,
               -4 * Ampere / Second,
               3 * Ampere / Second / Second,
               1 * Ampere / Second / Second / Second});
  auto const p = p2 * p3;
  {
    auto const actual = p(0 * Second);
    EXPECT_THAT(actual, AlmostEquals(2 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p(1 * Second);
    EXPECT_THAT(actual, AlmostEquals(-8 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p(-1 * Second);
    EXPECT_THAT(actual, AlmostEquals(-80 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p(2 * Second);
    EXPECT_THAT(actual, AlmostEquals(-350 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p(-2 * Second);
    EXPECT_THAT(actual, AlmostEquals(-518 * Ampere * Kelvin, 0));
  }
}

// Compose contains a fold expression which fails to compile in Clang because of
// https://bugs.llvm.org/show_bug.cgi?id=30590.  That bug will be fixed post-
// 11.0.0.  Since we don't use Compose as of this writing, and working around
// the bug would be hard, we ifdef out the test.
#if PRINCIPIA_COMPILER_MSVC
TEST_F(PolynomialTest, Monoid) {
  using P0 =
      PolynomialInMonomialBasis<Current, Temperature, 0, HornerEvaluator>;
  using P2A =
      PolynomialInMonomialBasis<Temperature, Instant, 2, HornerEvaluator>;
  using P2V =
      PolynomialInMonomialBasis<Temperature, Time, 2, HornerEvaluator>;
  using P3 =
      PolynomialInMonomialBasis<Current, Temperature, 3, HornerEvaluator>;
  Instant const t0;
  P0 const p0(std::tuple{9 * Ampere});
  P2A const p2a({1 * Kelvin,
                 3 * Kelvin / Second,
                 -8 * Kelvin / Second / Second}, t0);
  P2V const p2v({1 * Kelvin,
                 3 * Kelvin / Second,
                 -8 * Kelvin / Second / Second});
  P3 const p3({2 * Ampere,
               -4 * Ampere / Kelvin,
               3 * Ampere / Kelvin / Kelvin,
               1 * Ampere / Kelvin / Kelvin / Kelvin});
  auto const pa = Compose(p3, p2a);
  auto const pv = Compose(p3, p2v);
  {
    auto const actual_a = pa(t0 + 0 * Second);
    auto const actual_v = pv(0 * Second);
    EXPECT_THAT(actual_a, AlmostEquals(2 * Ampere, 0));
    EXPECT_THAT(actual_v, AlmostEquals(2 * Ampere, 0));
  }
  {
    auto const actual_a = pa(t0 + 1 * Second);
    auto const actual_v = pv(1 * Second);
    EXPECT_THAT(actual_a, AlmostEquals(2 * Ampere, 0));
    EXPECT_THAT(actual_v, AlmostEquals(2 * Ampere, 0));
  }
  {
    auto const actual_a = pa(t0 - 1 * Second);
    auto const actual_v = pv(-1 * Second);
    EXPECT_THAT(actual_a, AlmostEquals(-658 * Ampere, 0));
    EXPECT_THAT(actual_v, AlmostEquals(-658 * Ampere, 0));
  }
  {
    auto const actual_a = pa(t0 + 2 * Second);
    auto const actual_v = pv(2 * Second);
    EXPECT_THAT(actual_a, AlmostEquals(-13648 * Ampere, 0));
    EXPECT_THAT(actual_v, AlmostEquals(-13648 * Ampere, 0));
  }
  {
    auto const actual_a = pa(t0 - 2 * Second);
    auto const actual_v = pv(-2 * Second);
    EXPECT_THAT(actual_a, AlmostEquals(-46396 * Ampere, 0));
    EXPECT_THAT(actual_v, AlmostEquals(-46396 * Ampere, 0));
  }
  {
    auto const actual = Compose(p0, p2a)(t0);
    EXPECT_THAT(actual, AlmostEquals(9 * Ampere, 0));
  }
}
#endif

TEST_F(PolynomialTest, PointwiseInnerProduct) {
  P2V::Coefficients const coefficients({
      Displacement<World>({0 * Metre,
                           2 * Metre,
                           3 * Metre}),
      Velocity<World>({-1 * Metre / Second,
                       1 * Metre / Second,
                       0 * Metre / Second}),
      Vector<Acceleration, World>({1 * Metre / Second / Second,
                                   1 * Metre / Second / Second,
                                   -2 * Metre / Second / Second})});
  P2V const p2va(coefficients_);
  P2V const p2vb(coefficients);

  auto const p = PointwiseInnerProduct(p2va, p2vb);
  {
    auto const actual = p(0 * Second);
    EXPECT_THAT(actual, AlmostEquals(3 * Metre * Metre, 0));
  }
  {
    auto const actual = p(1 * Second);
    EXPECT_THAT(actual, AlmostEquals(5 * Metre * Metre, 0));
  }
  {
    auto const actual = p(-1 * Second);
    EXPECT_THAT(actual, AlmostEquals(1 * Metre * Metre, 0));
  }
  {
    auto const actual = p(2 * Second);
    EXPECT_THAT(actual, AlmostEquals(19 * Metre * Metre, 0));
  }
  {
    auto const actual = p(-2 * Second);
    EXPECT_THAT(actual, AlmostEquals(11 * Metre * Metre, 0));
  }
}

TEST_F(PolynomialTest, AtOrigin) {
  Instant const t0 = Instant() + 3 * Second;
  P2A const p(coefficients_, t0);
  P2A const q = p.AtOrigin(Instant() - 2 * Second);
  for (Instant t = Instant() - 10 * Second;
       t < Instant() + 10 * Second;
       t += 0.3 * Second) {
    EXPECT_THAT(q(t), AlmostEquals(p(t), 0, 942));
  }
}

TEST_F(PolynomialTest, Derivative) {
  using P2 = PolynomialInMonomialBasis<Temperature, Time, 2, HornerEvaluator>;
  using P3 = PolynomialInMonomialBasis<Current, Time, 3, HornerEvaluator>;
  P2 const p2({1 * Kelvin, 3 * Kelvin / Second, -8 * Kelvin / Second / Second});
  P3 const p3({2 * Ampere,
               -4 * Ampere / Second,
               3 * Ampere / Second / Second,
               1 * Ampere / Second / Second / Second});

  EXPECT_EQ(3 * Kelvin / Second,
            p2.Derivative<1>()(0 * Second));
  EXPECT_EQ(-16 * Kelvin / Second / Second,
            p2.Derivative<2>()(0 * Second));

  EXPECT_EQ(-4 * Ampere / Second,
            p3.Derivative<1>()(0 * Second));
  EXPECT_EQ(6 * Ampere / Second / Second,
            p3.Derivative<2>()(0 * Second));
  EXPECT_EQ(6 * Ampere / Second / Second / Second,
            p3.Derivative<3>()(0 * Second));
}

TEST_F(PolynomialTest, PrimitiveIntegrate) {
  using P2 = PolynomialInMonomialBasis<Temperature, Time, 2, HornerEvaluator>;
  P2 const p2({1 * Kelvin, 3 * Kelvin / Second, -8 * Kelvin / Second / Second});

  EXPECT_THAT(p2.Primitive()(0 * Second),
              AlmostEquals(0 * Kelvin * Second, 0));
  EXPECT_THAT(p2.Primitive()(1 * Second),
              AlmostEquals(-1.0 / 6.0 * Kelvin * Second, 5));
  EXPECT_THAT(p2.Primitive()(-1 * Second),
              AlmostEquals(19.0 / 6.0 * Kelvin * Second, 1));
  EXPECT_THAT(p2.Primitive()(2 * Second),
              AlmostEquals(-40.0 / 3.0 * Kelvin * Second, 1));

  EXPECT_THAT(p2.Integrate(-1 * Second, 2 * Second),
              AlmostEquals(-99.0 / 6.0 * Kelvin * Second, 3));
}

TEST_F(PolynomialTest, EvaluateConstant) {
  PolynomialInMonomialBasis<Entropy, Time, 0, HornerEvaluator> const
      horner_boltzmann(std::make_tuple(BoltzmannConstant));
  PolynomialInMonomialBasis<Entropy, Time, 0, EstrinEvaluator> const
      estrin_boltzmann(std::make_tuple(BoltzmannConstant));
  EXPECT_THAT(horner_boltzmann(1729 * Second), Eq(BoltzmannConstant));
  EXPECT_THAT(estrin_boltzmann(1729 * Second), Eq(BoltzmannConstant));
  EXPECT_THAT(horner_boltzmann.EvaluateDerivative(1729 * Second),
              Eq(0 * Watt / Kelvin));
  EXPECT_THAT(estrin_boltzmann.EvaluateDerivative(1729 * Second),
              Eq(0 * Watt / Kelvin));
}

TEST_F(PolynomialTest, EvaluateLinear) {
  PolynomialInMonomialBasis<Length, Time, 1, HornerEvaluator> const
      horner_light({0 * Metre, SpeedOfLight});
  PolynomialInMonomialBasis<Length, Time, 1, EstrinEvaluator> const
      estrin_light({0 * Metre, SpeedOfLight});
  constexpr Length light_second = Second * SpeedOfLight;
  EXPECT_THAT(horner_light(1729 * Second), Eq(1729 * light_second));
  EXPECT_THAT(estrin_light(1729 * Second), Eq(1729 * light_second));
  EXPECT_THAT(horner_light.EvaluateDerivative(1729 * Second), Eq(SpeedOfLight));
  EXPECT_THAT(estrin_light.EvaluateDerivative(1729 * Second), Eq(SpeedOfLight));
}

// Check that polynomials may be serialized.
TEST_F(PolynomialTest, Serialization) {
  {
    P2V p2v(coefficients_);
    serialization::Polynomial message;
    p2v.WriteToMessage(&message);
    EXPECT_EQ(2, message.degree());
    EXPECT_TRUE(message.HasExtension(
        serialization::PolynomialInMonomialBasis::extension));
    auto const& extension = message.GetExtension(
        serialization::PolynomialInMonomialBasis::extension);
    EXPECT_EQ(3, extension.coefficient_size());
    for (auto const& coefficient : extension.coefficient()) {
      EXPECT_TRUE(coefficient.has_multivector());
    }
    EXPECT_FALSE(extension.has_origin());

    auto const polynomial_read =
        Polynomial<Displacement<World>, Time>::ReadFromMessage<HornerEvaluator>(
            message);
    EXPECT_EQ(2, polynomial_read->degree());
    EXPECT_THAT(
        (*polynomial_read)(0.5 * Second),
        AlmostEquals(
            Displacement<World>({0.25 * Metre, 0.5 * Metre, 1 * Metre}), 0));
    serialization::Polynomial message2;
    polynomial_read->WriteToMessage(&message2);
    EXPECT_THAT(message2, EqualsProto(message));
  }
  {
    P2A p2a(coefficients_, Instant());
    serialization::Polynomial message;
    p2a.WriteToMessage(&message);
    EXPECT_EQ(2, message.degree());
    EXPECT_TRUE(message.HasExtension(
        serialization::PolynomialInMonomialBasis::extension));
    auto const& extension = message.GetExtension(
        serialization::PolynomialInMonomialBasis::extension);
    EXPECT_EQ(3, extension.coefficient_size());
    for (auto const& coefficient : extension.coefficient()) {
      EXPECT_TRUE(coefficient.has_multivector());
    }
    EXPECT_TRUE(extension.has_origin());
    EXPECT_TRUE(extension.origin().has_scalar());

    auto const polynomial_read =
        Polynomial<Displacement<World>,
                   Instant>::ReadFromMessage<HornerEvaluator>(message);
    EXPECT_EQ(2, polynomial_read->degree());
    EXPECT_THAT(
        (*polynomial_read)(Instant() + 0.5 * Second),
        AlmostEquals(
            Displacement<World>({0.25 * Metre, 0.5 * Metre, 1 * Metre}), 0));
    serialization::Polynomial message2;
    polynomial_read->WriteToMessage(&message2);
    EXPECT_THAT(message2, EqualsProto(message));
  }
  {
    P17::Coefficients const coefficients;
    P17 p17(coefficients);
    serialization::Polynomial message;
    p17.WriteToMessage(&message);
    EXPECT_EQ(17, message.degree());
    EXPECT_TRUE(message.HasExtension(
        serialization::PolynomialInMonomialBasis::extension));
    auto const& extension = message.GetExtension(
        serialization::PolynomialInMonomialBasis::extension);
    EXPECT_EQ(18, extension.coefficient_size());
    for (auto const& coefficient : extension.coefficient()) {
      EXPECT_TRUE(coefficient.has_multivector());
    }
    EXPECT_FALSE(extension.has_origin());

    auto const polynomial_read =
        Polynomial<Displacement<World>, Time>::ReadFromMessage<HornerEvaluator>(
            message);
    EXPECT_EQ(17, polynomial_read->degree());
    EXPECT_THAT((*polynomial_read)(0.5 * Second),
                AlmostEquals(
                    Displacement<World>({0 * Metre, 0 * Metre, 0 * Metre}), 0));
    serialization::Polynomial message2;
    polynomial_read->WriteToMessage(&message2);
    EXPECT_THAT(message2, EqualsProto(message));
  }
}

TEST_F(PolynomialTest, Output) {
  P2V p2v(coefficients_);
  P2A p2a(coefficients_, Instant());
  P17::Coefficients const coefficients;
  P17 p17(coefficients);
  LOG(ERROR) << p2v;
  LOG(ERROR) << p2a;
  LOG(ERROR) << p17;
}

}  // namespace numerics
}  // namespace principia

#undef PRINCIPIA_USE_IACA
