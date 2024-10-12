#include "numerics/sin_cos.hpp"

#include <algorithm>
#include <limits>
#include <random>

#include "boost/multiprecision/cpp_int.hpp"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "numerics/next.hpp"
#include "quantities/numbers.hpp"
#include "testing_utilities/almost_equals.hpp"

// This test lives in `functions` to avoid pulling `boost` into `numerics`.
namespace principia {
namespace numerics {
namespace _sin_cos {

using namespace boost::multiprecision;
using namespace principia::numerics::_next;
using namespace principia::testing_utilities::_almost_equals;

class SinCosTest : public ::testing::Test {};

TEST_F(SinCosTest, Random) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(-2 * π, 2 * π);

  cpp_bin_float_50 max_sin_ulps_error = 0;
  cpp_bin_float_50 max_cos_ulps_error = 0;
  double worst_sin_argument = 0;
  double worst_cos_argument = 0;
  std::int64_t incorrectly_rounded_sin = 0;
  std::int64_t incorrectly_rounded_cos = 0;

#if _DEBUG
  static constexpr std::int64_t iterations = 100;
#else
  static constexpr std::int64_t iterations = 4'000'000;
#endif

  for (std::int64_t i = 0; i < iterations; ++i) {
    double const principia_argument = uniformly_at(random);
    auto const boost_argument = cpp_rational(principia_argument);
    {
      auto const boost_sin =
          functions::_multiprecision::Sin(boost_argument);
      double const principia_sin = Sin(principia_argument);
      auto const sin_error =
          abs(boost_sin - static_cast<cpp_bin_float_50>(principia_sin));
      auto const ulp = NextUp(principia_sin) - principia_sin;
      auto const sin_ulps_error = sin_error / ulp;
      if (sin_ulps_error > max_sin_ulps_error) {
        max_sin_ulps_error = sin_ulps_error;
        worst_sin_argument = principia_argument;
      }
      if (sin_ulps_error > 0.5) {
        ++incorrectly_rounded_sin;
        LOG(ERROR) << "Sin: " << std::setprecision(25) << principia_argument;
      }
    }
    {
      auto const boost_cos =
          functions::_multiprecision::Cos(boost_argument);
      double const principia_cos = Cos(principia_argument);
      auto const cos_error =
          abs(boost_cos - static_cast<cpp_bin_float_50>(principia_cos));
      auto const ulp = NextUp(principia_cos) - principia_cos;
      auto const cos_ulps_error = cos_error / ulp;
      if (cos_ulps_error > max_cos_ulps_error) {
        max_cos_ulps_error = cos_ulps_error;
        worst_cos_argument = principia_argument;
      }
      if (cos_ulps_error > 0.5) {
        ++incorrectly_rounded_cos;
        LOG(ERROR) << "Cos: " << std::setprecision(25) << principia_argument;
      }
    }
  }

  // This implementation is not quite correctly rounded, but not far from it.
  EXPECT_LE(max_sin_ulps_error, 0.500002);
  EXPECT_LE(max_cos_ulps_error, 0.500002);
  EXPECT_LE(incorrectly_rounded_sin, 1);
  EXPECT_LE(incorrectly_rounded_cos, 1);

  LOG(ERROR) << "Sin error: " << max_sin_ulps_error << std::setprecision(25)
             << " ulps for argument: " << worst_sin_argument
             << " value: " << Sin(worst_sin_argument)
             << "; incorrectly rounded: " << std::setprecision(3)
             << incorrectly_rounded_sin / static_cast<double>(iterations);
  LOG(ERROR) << "Cos error: " << max_cos_ulps_error << std::setprecision(25)
             << " ulps for argument: " << worst_cos_argument
             << " value: " << Cos(worst_cos_argument)
             << "; incorrectly rounded: " << std::setprecision(3)
             << incorrectly_rounded_cos / static_cast<double>(iterations);
}

// Values for which the base algorithm gives an error of 1 ULP.
TEST_F(SinCosTest, HardRounding) {
  EXPECT_THAT(Sin(1.777288458404935767021016),
              AlmostEquals(0.9787561457198967196367773, 1));
  EXPECT_THAT(Cos(3.912491942337291916942377),
              AlmostEquals(-0.7172843528140595004137653, 1));
  EXPECT_THAT(Sin(5.528810471911395296729097),
              AlmostEquals(-0.6848332450871304488693046, 1));
  EXPECT_THAT(Sin(2.670333644894535396474566),
              AlmostEquals(0.4540084183741445456039384, 1));
  EXPECT_THAT(Cos(1.486604973422413600303571),
              AlmostEquals(0.0840919279825555407437241, 1));
  EXPECT_THAT(Sin(-2.496680544289484160458414),
              AlmostEquals(-0.6011282027544306294509797, 1));
  EXPECT_THAT(Sin(3.348980952786005715893225),
              AlmostEquals(-0.2059048676737040700634683, 1));
  EXPECT_THAT(Sin(3.523452575387961971387085),
              AlmostEquals(-0.3726470704519433130297035, 1));
  EXPECT_THAT(Cos(-6.265702600230396157598989),
              AlmostEquals(0.99984718137127853720984932, 1));
  EXPECT_THAT(Sin(1.881458091523454001503524),
              AlmostEquals(0.9521314843257784876761001, 1));
  EXPECT_THAT(Sin(-1.763163156774038675678185),
              AlmostEquals(-0.9815544881044536151825223, 1));
  EXPECT_THAT(Cos(-3.885819786017697730073905),
              AlmostEquals(-0.7356116652133562472394118, 1));
  EXPECT_THAT(Sin(-2.58105062034143273308473),
              AlmostEquals(-0.5316453603071467637339815, 1));
  EXPECT_THAT(Sin(1.657419885978818285821035),
              AlmostEquals(0.99625052493662308306103561, 1));
  EXPECT_THAT(Sin(5.094301519947547873812255),
              AlmostEquals(-0.9279535374988051033005616, 1));
  EXPECT_THAT(Sin(5.262137362438826571064965),
              AlmostEquals(-0.8526560125576488347044409, 1));
  EXPECT_THAT(Cos(-5.026994177012682030181168),
              AlmostEquals(0.3094410694753661206223057, 1));
  EXPECT_THAT(Cos(0.2388111698570396512764091),
              AlmostEquals(0.9716198764286143041422587, 1));
}

}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia