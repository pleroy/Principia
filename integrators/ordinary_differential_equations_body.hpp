#pragma once

#include "integrators/ordinary_differential_equations.hpp"

#include <vector>

#include "base/for_all_of.hpp"

namespace principia {
namespace integrators {
namespace termination_condition {

inline void UpdateWithAbort(absl::Status const& updater,
                            absl::Status& updated) {
  if (absl::IsAborted(updater)) {
    updated = updater;
  } else {
    updated.Update(updater);
  }
}

}  // namespace termination_condition

namespace internal_ordinary_differential_equations {

using base::for_all_of;

template<typename IndependentVariable_, typename... DependentVariable>
ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable_,
                                               DependentVariable...>::
State::State(IndependentVariable const& s, DependentVariables const& y)
    : s(s) {
  for_all_of(y, this->y).loop([](auto const& y, auto& this_y) {
    for (auto const& y_i : y) {
      this_y.emplace_back(y_i);
    }
  });
}

template<typename IndependentVariable_, typename... DependentVariable>
void ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable_,
                                                    DependentVariable...>::
State::WriteToMessage(not_null<serialization::State*> message) const {
  // Writing the tuple would be tricky.
  LOG(FATAL) << "NYI";
}

template<typename IndependentVariable_, typename... DependentVariable>
typename ExplicitFirstOrderOrdinaryDifferentialEquation<
    IndependentVariable_,
    DependentVariable...>::State
ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable_,
                                               DependentVariable...>::
State::ReadFromMessage(serialization::State const& message) {
  // Reading the tuple would be tricky.
  LOG(FATAL) << "NYI";
}

template<typename... DependentVariable>
DecomposableFirstOrderDifferentialEquation<DependentVariable...>::
State::State(IndependentVariable const& t, DependentVariables const& y)
    : time(t), y(y) {}

template<typename DependentVariable>
ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>::
State::State(IndependentVariable const& t,
             DependentVariables const& q,
             DependentVariableDerivatives const& v)
    : time(t) {
  for (int i = 0; i < q.size(); ++i) {
    positions.emplace_back(q[i]);
    velocities.emplace_back(v[i]);
  }
}

template<typename DependentVariable>
void ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>::
State::WriteToMessage(not_null<serialization::State*> const message) const {
  for (auto const& position : positions) {
    position.WriteToMessage(message->add_position());
  }
  for (auto const& velocity : velocities) {
    velocity.WriteToMessage(message->add_velocity());
  }
  time.WriteToMessage(message->mutable_time());
}

template<typename DependentVariable>
typename ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>::
    State
ExplicitSecondOrderOrdinaryDifferentialEquation<DependentVariable>::
State::ReadFromMessage(serialization::State const& message) {
  State state;
  for (auto const& p : message.position()) {
    state.positions.push_back(
        DoublePrecision<DependentVariable>::ReadFromMessage(p));
  }
  for (auto const& v : message.velocity()) {
    state.velocities.push_back(
        DoublePrecision<DependentVariableDerivative>::ReadFromMessage(v));
  }
  state.time =
      DoublePrecision<IndependentVariable>::ReadFromMessage(message.time());
  return state;
}

}  // namespace internal_ordinary_differential_equations
}  // namespace integrators
}  // namespace principia
