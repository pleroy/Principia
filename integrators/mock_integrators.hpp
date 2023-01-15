#pragma once

#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "integrators/integrators.hpp"

namespace principia {
namespace integrators {
namespace internal_integrators {

using base::make_not_null_unique;

template<typename DifferentialEquation>
class MockFixedStepSizeIntegrator
    : public FixedStepSizeIntegrator<DifferentialEquation> {
 public:
  using ODE = DifferentialEquation;
  using typename Integrator<ODE>::AppendState;

  class MockInstance : public Integrator<ODE>::Instance {
   public:
    MockInstance() : Integrator<ODE>::Instance() {}

    MOCK_METHOD(absl::Status, Solve, (Instant const& t_final), (override));
    MOCK_METHOD(not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>,
                Clone,
                (),
                (const, override));
    MOCK_METHOD(FixedStepSizeIntegrator<ODE> const&,
                integrator,
                (),
                (const, override));
  };

  MOCK_METHOD(not_null<std::unique_ptr<typename Integrator<ODE>::Instance>>,
              NewInstance,
              (InitialValueProblem<ODE> const& problem,
               typename Integrator<ODE>::AppendState const& append_state,
               Time const& step),
              (const, override));
  MOCK_METHOD(void,
              WriteToMessage,
              (not_null<serialization::FixedStepSizeIntegrator*> message),
              (const, override));

  static MockFixedStepSizeIntegrator const& Get() {
    static MockFixedStepSizeIntegrator const integrator;
    return integrator;
  }

 private:
  MockFixedStepSizeIntegrator() : FixedStepSizeIntegrator<ODE>() {}
};

}  // namespace internal_integrators

using internal_integrators::MockFixedStepSizeIntegrator;

}  // namespace integrators
}  // namespace principia
