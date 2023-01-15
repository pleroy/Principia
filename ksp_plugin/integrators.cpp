#include "ksp_plugin/integrators.hpp"

#include <limits>

#include "geometry/named_quantities.hpp"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_integrators {

using geometry::Position;
using integrators::EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using integrators::methods::BlanesMoan2002SRKN14A;
using integrators::methods::Fine1987RKNG34;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using integrators::methods::Quinlan1999Order8A;
using quantities::si::Minute;
using quantities::si::Second;

DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters
DefaultDownsamplingParameters() {
  return DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters{
      .max_dense_intervals = 10'000,
      .tolerance = 10 * Metre,
  };
}

Ephemeris<Barycentric>::AccuracyParameters
DefaultEphemerisAccuracyParameters() {
  return Ephemeris<Barycentric>::AccuracyParameters(
      /*fitting_tolerance=*/1 * Milli(Metre),
      /*geopotential_tolerance*/ 0x1.0p-24);
}

Ephemeris<Barycentric>::FixedStepParameters
DefaultEphemerisFixedStepParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
      SymplecticRungeKuttaNyströmIntegrator<
          BlanesMoan2002SRKN14A,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*step=*/35 * Minute);
}

Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
DefaultBurnParameters() {
  return Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters(
      EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
          Fine1987RKNG34,
          Ephemeris<Barycentric>::GeneralizedNewtonianMotionEquation>(),
      /*max_steps=*/1000,
      /*length_integration_tolerance=*/1 * Metre,
      /*speed_integration_tolerance=*/1 * Metre / Second);
}

Ephemeris<Barycentric>::FixedStepParameters DefaultHistoryParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
      SymmetricLinearMultistepIntegrator<
          Quinlan1999Order8A,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*step=*/10 * Second);
}

Ephemeris<Barycentric>::AdaptiveStepParameters
DefaultPsychohistoryParameters() {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          DormandالمكاوىPrince1986RKN434FM,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*length_integration_tolerance=*/1 * Milli(Metre),
      /*speed_integration_tolerance=*/1 * Milli(Metre) / Second);
}

Ephemeris<Barycentric>::AdaptiveStepParameters DefaultPredictionParameters() {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          DormandالمكاوىPrince1986RKN434FM,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*max_steps=*/1000,
      /*length_integration_tolerance=*/1 * Metre,
      /*speed_integration_tolerance=*/1 * Metre / Second);
}

}  // namespace internal_integrators
}  // namespace ksp_plugin
}  // namespace principia
