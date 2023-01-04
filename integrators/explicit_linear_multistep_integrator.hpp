// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_EXPLICIT_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#define PRINCIPIA_INTEGRATORS_EXPLICIT_LINEAR_MULTISTEP_INTEGRATOR_HPP_

#include <functional>
#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "numerics/fixed_arrays.hpp"
#include "integrators/methods.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/integrators.pb.h"

// TODO(phl): Unify the style of the various integrators, they have spurious
// differences at the moment.
namespace principia {
namespace integrators {
namespace internal_explicit_linear_multistep_integrator {

using base::not_null;
using numerics::FixedStrictlyLowerTriangularMatrix;
using numerics::FixedVector;
using quantities::Variation;

//TODO(phl):comment
// This class solves ordinary differential equations of the form q″ = f(q, t)
// using an embedded Runge-Kutta method.  We follow the standard conventions for
// the coefficients, i.e.,
//   c for the nodes;
//   a for the Runge-Kutta matrix;
//   b̂ for the weights of the high-order method;
//   b for the weights of the low-order method;
// See [DP86] for an example.

// In the implementation, we follow [DP807] in calling the results of the
// right-hand-side evaluations fᵢ (this quantity is not named in [DP86]).
// The order of the template parameters follow the notation of [DP86], whose
// RKq(p)s[F]X has higher order q, lower order p, comprises s stages, and has
// the first-same-as-last property.

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
class ExplicitLinearMultistepIntegrator
    : public FixedStepSizeIntegrator<
          ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable,
                                                         StateElements...>> {
 public:
  using ODE =
      ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable,
                                                     StateElements...>;
  using typename Integrator<ODE>::AppendState;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    absl::Status Solve(IndependentVariable const& s_final) override;
    ExplicitLinearMultistepIntegrator const& integrator() const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;
#if 0
    template<typename... S = StateElements...,
             typename = std::enable_if_t<base::is_serializable_v<S...>>>
    static not_null<std::unique_ptr<Instance>> ReadFromMessage(
        serialization::
            ExplicitLinearMultistepIntegratorInstance const& extension,
        IntegrationProblem<ODE> const& problem,
        AppendState const& append_state,
        IndependentVariableDifference const& step,
        ExplicitLinearMultistepIntegrator const& integrator);
#endif

   private:
    using typename ODE::IndependentVariableDifference;
    using typename ODE::State;
    using typename ODE::StateVariation;
    using typename ODE::SystemState;

    //TODO(phl)should be common?
    struct Step final {
      std::vector<DoublePrecision<State>> states;
      std::vector<StateVariation> state_variations;
      DoublePrecision<IndependentVariable> s;

      void WriteToMessage(
          not_null<serialization::ExplicitLinearMultistepIntegratorInstance::
                       Step*> message) const;
      template<typename P = Position,
               typename = std::enable_if_t<base::is_serializable_v<P>>>
      static Step ReadFromMessage(
          serialization::ExplicitLinearMultistepIntegratorInstance::Step const&
              message);
    };

    Instance(IntegrationProblem<ODE> const& problem,
             AppendState const& append_state,
             IndependentVariableDifference const& step,
             ExplicitLinearMultistepIntegrator const& integrator);

    // Performs the startup integration, i.e., computes enough states to either
    // reach |s_final| or to reach a point where |instance.previous_steps_| has
    // |order - 1| elements.  During startup |instance.current_state_| is
    // updated more frequently than once every |instance.step_|.
    void StartupSolve(IndependentVariable const& s_final);

    static void FillStepFromSystemState(ODE const& equation,
                                        SystemState const& state,
                                        Step& step);

    int startup_step_index_ = 0;
    std::list<Step> previous_steps_;  // At most |order_| elements.
    ExplicitLinearMultistepIntegrator const& integrator_;
    friend class ExplicitLinearMultistepIntegrator;
  };

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state) const override;

  explicit ExplicitLinearMultistepIntegrator(
      FixedStepSizeIntegrator<ODE> const& startup_integrator);

  ExplicitLinearMultistepIntegrator(
      ExplicitLinearMultistepIntegrator const&) = delete;
  ExplicitLinearMultistepIntegrator(
      ExplicitLinearMultistepIntegrator&&) = delete;
  ExplicitLinearMultistepIntegrator& operator=(
      ExplicitLinearMultistepIntegrator const&) = delete;
  ExplicitLinearMultistepIntegrator& operator=(
      ExplicitLinearMultistepIntegrator&&) = delete;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message)
      const override;

 private:
  static constexpr auto order_ = Method::order;
  static constexpr auto b_numerator_ = Method::b_numerator;
  static constexpr auto b_denominator_ = Method::b_denominator;

  FixedStepSizeIntegrator<ODE> const& startup_integrator_;
};

}  // namespace internal_explicit_linear_multistep_integrator

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
internal_explicit_linear_multistep_integrator::
    ExplicitLinearMultistepIntegrator<Method,
                                         IndependentVariable,
                                         StateElements...> const&
ExplicitLinearMultistepIntegrator();

}  // namespace integrators
}  // namespace principia

#include "integrators/explicit_linear_multistep_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_EXPLICIT_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
