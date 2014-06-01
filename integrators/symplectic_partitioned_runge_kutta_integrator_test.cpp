
#define TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/dimensionless.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/numerical_analysis.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

using principia::quantities::Abs;
using principia::quantities::AngularFrequency;
using principia::quantities::Cos;
using principia::quantities::Dimensionless;
using principia::quantities::Force;
using principia::quantities::Length;
using principia::quantities::Log10;
using principia::quantities::Max;
using principia::quantities::Momentum;
using principia::quantities::Sin;
using principia::quantities::Speed;
using principia::testing_utilities::AbsoluteError;
using principia::testing_utilities::BidimensionalDatasetMathematicaInput;
using principia::testing_utilities::ComputeHarmonicOscillatorForce;
using principia::testing_utilities::ComputeHarmonicOscillatorVelocity;
using principia::testing_utilities::PearsonProductMomentCorrelationCoefficient;
using principia::testing_utilities::Slope;
using testing::AllOf;
using testing::Gt;
using testing::Lt;

namespace principia {
namespace integrators {

class SPRKTest : public testing::Test {
 public:
  static void SetUpTestCase() {
    google::LogToStderr();
  }

 protected:
  void SetUp() override {
    integrator_.Initialize(integrator_.Order5Optimal());
  }

  SPRKIntegrator                               integrator_;
  SPRKIntegrator::Parameters<Length, Momentum> parameters_;
  SPRKIntegrator::Solution<Length, Momentum>   solution_;
};

TEST_F(SPRKTest, HarmonicOscillator) {
  parameters_.q0 = {1.0 * Length::SIUnit()};
  parameters_.p0 = {0.0 * Momentum::SIUnit()};
  parameters_.t0 = 0.0 * Time::SIUnit();
#ifdef _DEBUG
  parameters_.tmax = 100.0 * Time::SIUnit();
#else
  parameters_.tmax = 1000.0 * Time::SIUnit();
#endif
  parameters_.Δt = 1.0E-4 * Time::SIUnit();
  parameters_.sampling_period = 1;
  integrator_.Solve(&ComputeHarmonicOscillatorForce<Length, Momentum>,
                    &ComputeHarmonicOscillatorVelocity<Length, Momentum>,
                    parameters_, &solution_);
  Length   q_error = 0 * Length::SIUnit();
  Momentum p_error = 0 * Momentum::SIUnit();
  for (size_t i = 0; i < solution_.time.quantities.size(); ++i) {
    q_error = Max(
        q_error,
        Abs(solution_.position[0].quantities[i] -
            Length::SIUnit() * Cos(AngularFrequency::SIUnit() *
                                       solution_.time.quantities[i])));
    p_error = Max(
        p_error,
        Abs(solution_.momentum[0].quantities[i] +
            Momentum::SIUnit() * Sin(AngularFrequency::SIUnit() *
                                         solution_.time.quantities[i])));
  }
  LOG(INFO) << "q_error = " << q_error;
  LOG(INFO) << "p_error = " << p_error;
  EXPECT_THAT(q_error, Lt(2E-16 * Speed::SIUnit() * parameters_.tmax));
  EXPECT_THAT(p_error, Lt(2E-16 * Force::SIUnit() * parameters_.tmax));
}

TEST_F(SPRKTest, Convergence) {
  parameters_.q0 = {1.0 * Length::SIUnit()};
  parameters_.p0 = {0.0 * Momentum::SIUnit()};
  parameters_.t0 = 0.0 * Time::SIUnit();
  parameters_.tmax = 100 * Time::SIUnit();
  parameters_.sampling_period = 0;
  // For 0.2 * 1.1⁻²¹ < |Δt| < 0.2 , the correlation between step size and error
  // is very strong. It the step is small enough to converge and large enough to
  // stay clear of floating point inaccuracy.
  parameters_.Δt = 0.2 * Time::SIUnit();
  int const step_sizes = 22;
  double const step_reduction = 1.1;
  std::vector<Dimensionless> log_step_sizes(step_sizes);
  std::vector<Dimensionless> log_q_errors(step_sizes);
  std::vector<Dimensionless> log_p_errors(step_sizes);
  for (int i = 0; i < step_sizes; ++i, parameters_.Δt /= step_reduction) {
    integrator_.Solve(&ComputeHarmonicOscillatorForce<Length, Momentum>,
                      &ComputeHarmonicOscillatorVelocity<Length, Momentum>,
                      parameters_, &solution_);
    log_step_sizes[i] = Log10(parameters_.Δt / (1 * Time::SIUnit()));
    log_q_errors[i] = Log10(
        Abs((solution_.position[0].quantities[0] -
             Length::SIUnit() * Cos(AngularFrequency::SIUnit() *
                                        solution_.time.quantities[0])) /
                 Length::SIUnit()));
    log_p_errors[i] = Log10(
        Abs((solution_.momentum[0].quantities[0] +
             Momentum::SIUnit() * Sin(AngularFrequency::SIUnit() *
                                          solution_.time.quantities[0])) /
                 Momentum::SIUnit()));
  }
  Dimensionless const q_convergence_order = Slope(log_step_sizes, log_q_errors);
  Dimensionless const q_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_q_errors);
  LOG(INFO) << "Convergence order in q : " << q_convergence_order;
  LOG(INFO) << "Correlation            : " << q_correlation;
  LOG(INFO) << "Convergence data for q :\n" <<
      BidimensionalDatasetMathematicaInput(log_step_sizes, log_q_errors);
  EXPECT_THAT(q_convergence_order, AllOf(Gt(4.9), Lt(5.1)));
  EXPECT_THAT(q_correlation, AllOf(Gt(0.999), Lt(1.01)));
  Dimensionless const p_convergence_order = Slope(log_step_sizes, log_p_errors);
  Dimensionless const p_correlation =
      PearsonProductMomentCorrelationCoefficient(log_step_sizes, log_p_errors);
  LOG(INFO) << "Convergence order in p : " << p_convergence_order;
  LOG(INFO) << "Correlation            : " << p_correlation;
  LOG(INFO) << "Convergence data for p :\n" <<
      BidimensionalDatasetMathematicaInput(log_step_sizes, log_q_errors);
  EXPECT_THAT(p_convergence_order, AllOf(Gt(5.9), Lt(6.1)));
  EXPECT_THAT(p_correlation, AllOf(Gt(0.999), Lt(1.01)));
}

}  // namespace integrators
}  // namespace principia
