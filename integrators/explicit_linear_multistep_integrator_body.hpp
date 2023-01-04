#pragma once

#include "integrators/explicit_linear_multistep_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <optional>
#include <vector>

#include "base/for_all_of.hpp"
#include "base/jthread.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_explicit_linear_multistep_integrator {

using base::for_all_of;
using base::make_not_null_unique;
using geometry::Sign;
using numerics::DoublePrecision;
using quantities::DebugString;
using quantities::Difference;
using quantities::Quotient;

int const startup_step_divisor = 16;

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
absl::Status
ExplicitLinearMultistepIntegrator<Method,
                                  IndependentVariable,
                                  StateElements...>::Instance::
Solve(IndependentVariable const& s_final) {
  using IndependentVariableDifference =
      typename ODE::IndependentVariableDifference;
  using State = typename ODE::State;
  using StateDifference = typename ODE::StateDifference;
  using StateVariation = typename ODE::StateVariation;

  auto const& b_numerator = integrator_.b_numerator_;
  auto const& b_denominator = integrator_.b_denominator_;

  auto& current_state = this->current_state_;
  auto& append_state = this->append_state_;
  auto const& step = this->step_;
  auto const& equation = this->equation_;

  if (previous_steps_.size() < order) {
    StartupSolve(t_final);

    // If |s_final| is not large enough, we may not have generated enough
    // points.  Bail out, we'll continue the next time |Solve| is called.
    if (previous_steps_.size() < order) {
      return absl::OkStatus();
    }
  }
  CHECK_EQ(previous_steps_.size(), order);

  // Argument checks.
  int const dimension = previous_steps_.back().displacements.size();

  // Independent variable step.
  CHECK_LT(IndependentVariable(), step);
  IndependentVariable const& h = step;
  // Current independent variable.
  DoublePrecision<IndependentVariable> s = previous_steps_.back().s;
  // Order.
  int const k = order;

  absl::Status status;
  std::vector<Position> positions(dimension);

  DoubleDisplacements Σⱼ_minus_αⱼ_qⱼ(dimension);
  std::vector<Acceleration> Σⱼ_βⱼ_numerator_aⱼ(dimension);
  while (h <= (t_final - t.value) - t.error) {
    // We take advantage of the symmetry to iterate on the list of previous
    // steps from both ends.
    auto front_it = previous_steps_.begin();
    auto back_it = previous_steps_.rbegin();

    // This block corresponds to j = 0.  We must not pair it with j = k.
    {
      DoubleDisplacements const& qⱼ = front_it->displacements;
      std::vector<Acceleration> const& aⱼ = front_it->accelerations;
      double const αⱼ = α[0];
      double const βⱼ_numerator = β_numerator[0];
      for (int d = 0; d < dimension; ++d) {
        Σⱼ_minus_αⱼ_qⱼ[d] = Scale(-αⱼ, qⱼ[d]);
        Σⱼ_βⱼ_numerator_aⱼ[d] = βⱼ_numerator * aⱼ[d];
      }
      ++front_it;
    }
    // The generic value of j, paired with k - j.
    for (int j = 1; j < k / 2; ++j) {
      DoubleDisplacements const& qⱼ = front_it->displacements;
      DoubleDisplacements const& qₖ₋ⱼ = back_it->displacements;
      std::vector<Acceleration> const& aⱼ = front_it->accelerations;
      std::vector<Acceleration> const& aₖ₋ⱼ = back_it->accelerations;
      double const αⱼ = α[j];
      double const βⱼ_numerator = β_numerator[j];
      for (int d = 0; d < dimension; ++d) {
        Σⱼ_minus_αⱼ_qⱼ[d] -= Scale(αⱼ, qⱼ[d]);
        Σⱼ_minus_αⱼ_qⱼ[d] -= Scale(αⱼ, qₖ₋ⱼ[d]);
        Σⱼ_βⱼ_numerator_aⱼ[d] += βⱼ_numerator * (aⱼ[d] + aₖ₋ⱼ[d]);
      }
      ++front_it;
      ++back_it;
    }
    // This block corresponds to j = k / 2.  We must not pair it with j = k / 2.
    {
      DoubleDisplacements const& qⱼ = front_it->displacements;
      std::vector<Acceleration> const& aⱼ = front_it->accelerations;
      double const αⱼ = α[k / 2];
      double const βⱼ_numerator = β_numerator[k / 2];
      for (int d = 0; d < dimension; ++d) {
        Σⱼ_minus_αⱼ_qⱼ[d] -= Scale(αⱼ, qⱼ[d]);
        Σⱼ_βⱼ_numerator_aⱼ[d] += βⱼ_numerator * aⱼ[d];
      }
    }

    // Create a new step in the instance.
    t.Increment(h);
    previous_steps_.emplace_back();
    Step& current_step = previous_steps_.back();
    current_step.time = t;
    current_step.accelerations.resize(dimension);
    current_step.displacements.reserve(dimension);

    // Fill the new step.  We skip the division by αₖ as it is equal to 1.0.
    double const αₖ = α[0];
    DCHECK_EQ(αₖ, 1.0);
    for (int d = 0; d < dimension; ++d) {
      DoubleDisplacement& current_displacement = Σⱼ_minus_αⱼ_qⱼ[d];
      current_displacement.Increment(h * h *
                                     Σⱼ_βⱼ_numerator_aⱼ[d] / β_denominator);
      current_step.displacements.push_back(current_displacement);
      DoublePosition const current_position =
          DoublePosition() + current_displacement;
      positions[d] = current_position.value;
      current_state.positions[d] = current_position;
    }
    termination_condition::UpdateWithAbort(
        equation.compute_acceleration(t.value,
                                      positions,
                                      current_step.accelerations),
        status);
    previous_steps_.pop_front();

    ComputeVelocityUsingCohenHubbardOesterwinter();

    // Inform the caller of the new state.
    RETURN_IF_STOPPED;
    current_state.time = t;
    append_state(current_state);
    if (absl::IsAborted(status)) {
      return status;
    }
  }

  return status;
}

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
ExplicitLinearMultistepIntegrator<Method,
                                  IndependentVariable,
                                  StateElements...> const&
ExplicitLinearMultistepIntegrator<Method,
                                     IndependentVariable,
                                     StateElements...>::Instance::
integrator() const {
  return integrator_;
}

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
not_null<std::unique_ptr<
    typename Integrator<ExplicitFirstOrderOrdinaryDifferentialEquation<
        IndependentVariable, StateElements...>>::Instance>>
ExplicitLinearMultistepIntegrator<Method,
                                  IndependentVariable,
                                  StateElements...>::Instance::
Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
void ExplicitLinearMultistepIntegrator<Method,
                                      IndependentVariable,
                                      StateElements...>::Instance::
WriteToMessage(not_null<serialization::IntegratorInstance*> message) const {
  FixedStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  [[maybe_unused]] auto* const extension =
      message
          ->MutableExtension(
              serialization::FixedStepSizeIntegratorInstance::extension)
          ->MutableExtension(
              serialization::
                  ExplicitLinearMultistepIntegratorInstance::extension);
}

#if 0
template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
template<typename, typename>
not_null<std::unique_ptr<
    typename ExplicitLinearMultistepIntegrator<Method,
                                                  IndependentVariable,
                                                  StateElements...>::Instance>>
ExplicitLinearMultistepIntegrator<Method,
                                  IndependentVariable,
                                  StateElements...>::Instance::
ReadFromMessage(serialization::
                    ExplicitLinearMultistepIntegratorInstance const&
                        extension,
                IntegrationProblem<ODE> const& problem,
                AppendState const& append_state,
                Time const& time_step,
                bool const first_use,
                ExplicitLinearMultistepIntegrator const& integrator) {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(new Instance(problem,
                                                append_state,
                                                time_step,
                                                first_use,
                                                integrator));
}
#endif

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
ExplicitLinearMultistepIntegrator<Method,
                                  IndependentVariable,
                                  StateElements...>::Instance::
Instance(IntegrationProblem<ODE> const& problem,
         AppendState const& append_state,
         IndependentVariableDifference const& step,
         ExplicitLinearMultistepIntegrator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem,
                                             append_state,
                                             step),
      integrator_(integrator) {}

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
void ExplicitLinearMultistepIntegrator<Method,
                                       IndependentVariable,
                                       StateElements...>::Instance::
Instance::StartupSolve(IndependentVariable const& s_final) {
  auto& current_state = this->current_state_;
  auto const& step = this->step_;
  auto const& equation = this->equation_;

  IndependentVariable const startup_step = step / startup_step_divisor;

  CHECK(!previous_steps_.empty());
  CHECK_LT(previous_steps_.size(), order);

  auto const startup_append_state =
      [this](SystemState const& state) {
        // Stop changing anything once we're done with the startup.  We may be
        // called one more time by the |startup_integrator_|.
        if (previous_steps_.size() < order) {
          this->current_state_ = state;
          // The startup integrator has a smaller step.  We do not record all
          // the states it computes, but only those that are a multiple of the
          // main integrator step.
          if (++startup_step_index_ % startup_step_divisor == 0) {
            CHECK_LT(previous_steps_.size(), order);
            previous_steps_.emplace_back();
            FillStepFromSystemState(this->equation_,
                                    this->current_state_,
                                    previous_steps_.back());
            // This call must happen last for a subtle reason: the callback may
            // want to |Clone| this instance (see |Ephemeris::Checkpoint|) in
            // which cases it is necessary that all the member variables be
            // filled for restartability to work.
            this->append_state_(state);
          }
        }
      };

  auto const startup_instance =
      integrator_.startup_integrator_.NewInstance({equation, current_state},
                                                  startup_append_state,
                                                  startup_step);

  startup_instance->Solve(
      std::min(current_state.time.value +
                   (order - previous_steps_.size()) * step + step / 2.0,
               t_final)).IgnoreError();

  CHECK_LE(previous_steps_.size(), order);
}

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
void ExplicitLinearMultistepIntegrator<Method,
                                       IndependentVariable,
                                       StateElements...>::
Instance::FillStepFromSystemState(ODE const& equation,
                                  SystemState const& state,
                                  Step& step) {
  std::vector<typename ODE::Position> positions;
  step.time = state.time;
  for (auto const& position : state.positions) {
    step.displacements.push_back(position - DoublePrecision<Position>());
    positions.push_back(position.value);
  }
  step.accelerations.resize(step.displacements.size());
  // Ignore the status here.  We are merely computing the acceleration to store
  // it, not to advance an integrator.
  equation.compute_acceleration(step.time.value,
                                positions,
                                step.accelerations).IgnoreError();
}

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
not_null<std::unique_ptr<
    typename Integrator<ExplicitFirstOrderOrdinaryDifferentialEquation<
        IndependentVariable, StateElements...>>::Instance>>
ExplicitLinearMultistepIntegrator<Method,
                                  IndependentVariable,
                                  StateElements...>::
NewInstance(IntegrationProblem<ODE> const& problem,
            AppendState const& append_state) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem,
                   append_state,
                   /*step=*/parameters.first_step,
                   *this));
}

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
ExplicitLinearMultistepIntegrator<Method,
                                  IndependentVariable,
                                  StateElements...>::
ExplicitLinearMultistepIntegrator(
    FixedStepSizeIntegrator<ODE> const& startup_integrator)
    : startup_integrator_(startup_integrator) {
}

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
void ExplicitLinearMultistepIntegrator<Method,
                                       IndependentVariable,
                                       StateElements...>::
WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  message->set_kind(Method::kind);
}

}  // namespace internal_explicit_linear_multistep_integrator

template<typename Method,
         typename IndependentVariable,
         typename... StateElements>
internal_explicit_linear_multistep_integrator::
    ExplicitLinearMultistepIntegrator<Method,
                                      IndependentVariable,
                                      StateElements...> const&
ExplicitLinearMultistepIntegrator() {
  static_assert(
      std::is_base_of<methods::ExplicitLinearMultistep,
                      Method>::value,
      "Method must be derived from ExplicitLinearMultistep");
  static internal_explicit_linear_multistep_integrator::
      ExplicitLinearMultistepIntegrator<Method,
                                        IndependentVariable,
                                        StateElements...> const integrator;
  return integrator;
}

}  // namespace integrators
}  // namespace principia
