﻿
#include "ksp_plugin/flight_plan.hpp"

#include <experimental/optional>
#include <vector>

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {

using integrators::DormandElMikkawyPrince1986RKN434FM;

namespace ksp_plugin {

FlightPlan::FlightPlan(
    Mass const& initial_mass,
    Instant const& initial_time,
    DegreesOfFreedom<Barycentric> const& initial_degrees_of_freedom,
    Instant const& final_time,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters)
    : initial_mass_(initial_mass),
      initial_time_(initial_time),
      initial_degrees_of_freedom_(initial_degrees_of_freedom),
      final_time_(final_time),
      root_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      ephemeris_(ephemeris),
      adaptive_step_parameters_(adaptive_step_parameters) {
  CHECK(final_time_ >= initial_time_);

  // Set the (single) point of the root.
  root_->Append(initial_time_, initial_degrees_of_freedom_);

  // Create a fork for the first coasting trajectory.
  segments_.emplace_back(root_->NewForkWithoutCopy(initial_time_));
  CoastLastSegment(final_time_);
}

Instant FlightPlan::initial_time() const {
  return initial_time_;
}

Instant FlightPlan::final_time() const {
  return final_time_;
}

int FlightPlan::number_of_manœuvres() const {
  return manœuvres_.size();
}

NavigationManœuvre const& FlightPlan::GetManœuvre(int const index) const {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_manœuvres());
  return manœuvres_[index];
}

bool FlightPlan::Append(Burn burn) {
  auto manœuvre =
      MakeNavigationManœuvre(
          std::move(burn),
          manœuvres_.empty() ? initial_mass_ : manœuvres_.back().final_mass());
  if (manœuvre.FitsBetween(start_of_last_coast(), final_time_) &&
      !manœuvre.IsSingular()) {
    DiscreteTrajectory<Barycentric>* recomputed_last_coast =
        CoastIfReachesManœuvreInitialTime(last_coast(), manœuvre);
    if (recomputed_last_coast != nullptr) {
      ReplaceLastSegment(recomputed_last_coast);
      Append(std::move(manœuvre));
      return true;
    }
  }
  return false;
}

void FlightPlan::ForgetBefore(Instant const& time,
                              std::function<void()> const& on_empty) {
  // Find the first segment to keep.  Note that incrementing by 2 ensures that
  // we only look at coasts.
  std::experimental::optional<int> first_to_keep;
  for (int i = 0; i < segments_.size(); i += 2) {
    if (time <= segments_[i]->last().time()) {
      first_to_keep = i;
      break;
    }
  }
  if (!first_to_keep) {
    // The entire flight plan needs to go away.
    on_empty();
    return;
  }

  // Detach the first coast to keep, truncate its beginning, and reattach it
  // to a new root.
  std::unique_ptr<DiscreteTrajectory<Barycentric>> new_first_coast =
      segments_[*first_to_keep]->DetachFork();
  new_first_coast->ForgetBefore(time);
  root_ = make_not_null_unique<DiscreteTrajectory<Barycentric>>();
  root_->AttachFork(std::move(new_first_coast));

  // Remove from the vectors the trajectories and manœuvres that we don't want
  // to keep.
  segments_.erase(segments_.cbegin(),
                  segments_.cbegin() + *first_to_keep);
  // For some reason manœuvres_.erase() doesn't work because it wants to
  // copy, hence this dance.
  std::vector<NavigationManœuvre> m;
  std::move(manœuvres_.begin() + *first_to_keep / 2,
            manœuvres_.end(),
            std::back_inserter(m));
  manœuvres_.swap(m);

  auto const root_begin = root_->Begin();
  initial_time_ = root_begin.time();
  initial_degrees_of_freedom_ = root_begin.degrees_of_freedom();
}

void FlightPlan::RemoveLast() {
  CHECK(!manœuvres_.empty());
  manœuvres_.pop_back();
  PopLastSegment();  // Last coast.
  PopLastSegment();  // Last burn.
  ResetLastSegment();
  CoastLastSegment(final_time_);
}

bool FlightPlan::ReplaceLast(Burn burn) {
  CHECK(!manœuvres_.empty());
  auto manœuvre = MakeNavigationManœuvre(std::move(burn),
                                         manœuvres_.back().initial_mass());
  if (manœuvre.FitsBetween(start_of_penultimate_coast(), final_time_) &&
      !manœuvre.IsSingular()) {
    DiscreteTrajectory<Barycentric>* recomputed_penultimate_coast =
        CoastIfReachesManœuvreInitialTime(penultimate_coast(), manœuvre);
    if (recomputed_penultimate_coast != nullptr) {
      manœuvres_.pop_back();
      PopLastSegment();  // Last coast.
      PopLastSegment();  // Last burn.
      ReplaceLastSegment(recomputed_penultimate_coast);
      Append(std::move(manœuvre));
      return true;
    }
  }
  return false;
}

bool FlightPlan::SetFinalTime(Instant const& final_time) {
  if (start_of_last_coast() > final_time) {
    return false;
  } else {
    final_time_ = final_time;
    ResetLastSegment();
    CoastLastSegment(final_time_);
    return true;
  }
}

Ephemeris<Barycentric>::AdaptiveStepParameters const&
FlightPlan::adaptive_step_parameters() const {
  return adaptive_step_parameters_;
}

bool FlightPlan::SetAdaptiveStepParameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters) {
  auto const original_adaptive_step_parameters = adaptive_step_parameters_;
  adaptive_step_parameters_ = adaptive_step_parameters;
  if (RecomputeSegments()) {
    return true;
  } else {
    // If the recomputation fails, leave this place as clean as we found it.
    adaptive_step_parameters_ = original_adaptive_step_parameters;
    CHECK(RecomputeSegments());
    return false;
  }
}

int FlightPlan::number_of_segments() const {
  return segments_.size();
}

void FlightPlan::GetSegment(
    int const index,
    not_null<DiscreteTrajectory<Barycentric>::Iterator*> begin,
    not_null<DiscreteTrajectory<Barycentric>::Iterator*> end) const {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_segments());
  *begin = segments_[index]->Fork();
  *end = segments_[index]->End();
}

void FlightPlan::GetAllSegments(
    not_null<DiscreteTrajectory<Barycentric>::Iterator*> begin,
    not_null<DiscreteTrajectory<Barycentric>::Iterator*> end) const {
  *begin = segments_.back()->Find(segments_.front()->Fork().time());
  *end = segments_.back()->End();
  CHECK(*begin != *end);
}

void FlightPlan::WriteToMessage(
    not_null<serialization::FlightPlan*> const message) const {
  initial_mass_.WriteToMessage(message->mutable_initial_mass());
  initial_time_.WriteToMessage(message->mutable_initial_time());
  initial_degrees_of_freedom_.WriteToMessage(
      message->mutable_initial_degrees_of_freedom());
  final_time_.WriteToMessage(message->mutable_final_time());
  adaptive_step_parameters_.WriteToMessage(
      message->mutable_adaptive_step_parameters());
  for (auto const& manœuvre : manœuvres_) {
    manœuvre.WriteToMessage(message->add_manoeuvre());
  }
}

std::unique_ptr<FlightPlan> FlightPlan::ReadFromMessage(
    serialization::FlightPlan const& message,
    not_null<DiscreteTrajectory<Barycentric>*> const root,
    not_null<Ephemeris<Barycentric>*> const ephemeris) {
  bool const is_pre_буняковский = message.segment_size() > 0;

  std::unique_ptr<Ephemeris<Barycentric>::AdaptiveStepParameters>
      adaptive_step_parameters;
  Instant initial_time = Instant::ReadFromMessage(message.initial_time());
  std::unique_ptr<DegreesOfFreedom<Barycentric>> initial_degrees_of_freedom;
  if (is_pre_буняковский) {
    adaptive_step_parameters =
        std::make_unique<Ephemeris<Barycentric>::AdaptiveStepParameters>(
            AdaptiveStepSizeIntegrator<
                Ephemeris<Barycentric>::NewtonianMotionEquation>::
                    ReadFromMessage(message.integrator()),
            /*max_steps=*/1000,
            Length::ReadFromMessage(message.length_integration_tolerance()),
            Speed::ReadFromMessage(message.speed_integration_tolerance()));
    auto it = root->LowerBound(initial_time);
    if (it.time() != initial_time) {
      --it;
      initial_time = it.time();
    }
    initial_degrees_of_freedom =
        std::make_unique<DegreesOfFreedom<Barycentric>>(
            it.degrees_of_freedom());
  } else {
    CHECK(message.has_adaptive_step_parameters());
    adaptive_step_parameters =
        std::make_unique<Ephemeris<Barycentric>::AdaptiveStepParameters>(
            Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
                message.adaptive_step_parameters()));
    if (message.has_initial_degrees_of_freedom()) {
      initial_degrees_of_freedom =
          std::make_unique<DegreesOfFreedom<Barycentric>>(
              DegreesOfFreedom<Barycentric>::ReadFromMessage(
                  message.initial_degrees_of_freedom()));
    } else {
      auto it = root->LowerBound(initial_time);
      if (it.time() != initial_time) {
        --it;
        initial_time = it.time();
      }
      initial_degrees_of_freedom =
          std::make_unique<DegreesOfFreedom<Barycentric>>(
              it.degrees_of_freedom());
    }
  }

  auto flight_plan = std::make_unique<FlightPlan>(
      Mass::ReadFromMessage(message.initial_mass()),
      initial_time,
      *initial_degrees_of_freedom,
      Instant::ReadFromMessage(message.final_time()),
      ephemeris,
      *adaptive_step_parameters);

  if (is_pre_буняковский) {
    // The constructor has forked a segment.  Remove it.
    flight_plan->PopLastSegment();
    for (auto const& segment : message.segment()) {
      flight_plan->segments_.emplace_back(
          DiscreteTrajectory<Barycentric>::ReadPointerFromMessage(
              segment, root));
    }
    for (int i = 0; i < message.manoeuvre_size(); ++i) {
      auto const& manoeuvre = message.manoeuvre(i);
      flight_plan->manœuvres_.push_back(
          NavigationManœuvre::ReadFromMessage(manoeuvre, ephemeris));
      flight_plan->manœuvres_[i].set_coasting_trajectory(
          flight_plan->segments_[2 * i]);
    }

    // We may end up here with a flight plan that has too many anomalous
    // segments because of past bugs.  The best we can do is to ignore it.
    if (!flight_plan->RecomputeSegments()) {
      flight_plan.reset();
    }
  } else {
    for (int i = 0; i < message.manoeuvre_size(); ++i) {
      auto const& manoeuvre = message.manoeuvre(i);
      flight_plan->manœuvres_.push_back(
          NavigationManœuvre::ReadFromMessage(manoeuvre, ephemeris));
    }
    // We need to forcefully prolong, otherwise we might exceed the ephemeris
    // step limit while recomputing the segments and fail the check.
    flight_plan->ephemeris_->Prolong(flight_plan->start_of_last_coast());
    CHECK(flight_plan->RecomputeSegments()) << message.DebugString();
  }

  return std::move(flight_plan);
}

FlightPlan::FlightPlan()
    : initial_degrees_of_freedom_(Barycentric::origin, Velocity<Barycentric>()),
      root_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()),
      adaptive_step_parameters_(
          DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
          /*max_steps=*/1,
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second) {}

void FlightPlan::Append(NavigationManœuvre manœuvre) {
  manœuvres_.emplace_back(std::move(manœuvre));
  {
    // Hide the moved-from |manœuvre|.
    NavigationManœuvre& manœuvre = manœuvres_.back();
    CHECK_EQ(manœuvre.initial_time(), segments_.back()->last().time());
    manœuvre.set_coasting_trajectory(segments_.back());
    AddSegment();
    BurnLastSegment(manœuvre);
    AddSegment();
    CoastLastSegment(final_time_);
  }
}

bool FlightPlan::RecomputeSegments() {
  // It is important that the segments be destroyed in (reverse chronological)
  // order of the forks.
  while (segments_.size() > 1) {
    PopLastSegment();
  }
  ResetLastSegment();
  for (auto& manœuvre : manœuvres_) {
    CoastLastSegment(manœuvre.initial_time());
    manœuvre.set_coasting_trajectory(segments_.back());
    AddSegment();
    BurnLastSegment(manœuvre);
    AddSegment();
  }
  CoastLastSegment(final_time_);
  return anomalous_segments_ <= 2;
}

void FlightPlan::BurnLastSegment(NavigationManœuvre const& manœuvre) {
  if (anomalous_segments_ > 0) {
    return;
  } else if (manœuvre.initial_time() < manœuvre.final_time()) {
    bool const reached_final_time =
        ephemeris_->FlowWithAdaptiveStep(segments_.back(),
                                         manœuvre.IntrinsicAcceleration(),
                                         manœuvre.final_time(),
                                         adaptive_step_parameters_,
                                         max_ephemeris_steps_per_frame);
    if (!reached_final_time) {
      anomalous_segments_ = 1;
    }
  }
}

void FlightPlan::CoastLastSegment(Instant const& final_time) {
  if (anomalous_segments_ > 0) {
    return;
  } else {
    bool const reached_final_time =
        ephemeris_->FlowWithAdaptiveStep(
                        segments_.back(),
                        Ephemeris<Barycentric>::kNoIntrinsicAcceleration,
                        final_time,
                        adaptive_step_parameters_,
                        max_ephemeris_steps_per_frame);
    if (!reached_final_time) {
      anomalous_segments_ = 1;
    }
  }
}

void FlightPlan::ReplaceLastSegment(
    not_null<DiscreteTrajectory<Barycentric>*> const segment) {
  CHECK_EQ(segment->parent(), segments_.back()->parent());
  CHECK_EQ(segment->Fork().time(), segments_.back()->Fork().time());
  PopLastSegment();
  // |segment| must not be anomalous, so it cannot not follow an anomalous
  // segment.
  CHECK_EQ(0, anomalous_segments_);
  segments_.emplace_back(segment);
}

void FlightPlan::AddSegment() {
  segments_.emplace_back(segments_.back()->NewForkAtLast());
  if (anomalous_segments_ > 0) {
    ++anomalous_segments_;
  }
}

void FlightPlan::ResetLastSegment() {
  segments_.back()->ForgetAfter(segments_.back()->Fork().time());
  if (anomalous_segments_ == 1) {
    // If there was one anomalous segment, it was the last one, which was
    // anomalous because it ended early.  It is no longer anomalous.
    anomalous_segments_ = 0;
  }
}

void FlightPlan::PopLastSegment() {
  DiscreteTrajectory<Barycentric>* trajectory = segments_.back();
  CHECK(!trajectory->is_root());
  trajectory->parent()->DeleteFork(&trajectory);
  segments_.pop_back();
  if (anomalous_segments_ > 0) {
    --anomalous_segments_;
  }
}

DiscreteTrajectory<Barycentric>* FlightPlan::CoastIfReachesManœuvreInitialTime(
    DiscreteTrajectory<Barycentric>& coast,
    NavigationManœuvre const& manœuvre) {
  DiscreteTrajectory<Barycentric>* recomputed_coast =
      coast.parent()->NewForkWithoutCopy(coast.Fork().time());
  bool const reached_manœuvre_initial_time =
      ephemeris_->FlowWithAdaptiveStep(
          recomputed_coast,
          Ephemeris<Barycentric>::kNoIntrinsicAcceleration,
          manœuvre.initial_time(),
          adaptive_step_parameters_,
          max_ephemeris_steps_per_frame);
  if (!reached_manœuvre_initial_time) {
    recomputed_coast->parent()->DeleteFork(&recomputed_coast);
  }
  return recomputed_coast;
}

Instant FlightPlan::start_of_last_coast() const {
  return manœuvres_.empty() ? initial_time_ : manœuvres_.back().final_time();
}

Instant FlightPlan::start_of_penultimate_coast() const {
  return manœuvres_.size() == 1
             ? initial_time_
             : manœuvres_[manœuvres_.size() - 2].final_time();
}

DiscreteTrajectory<Barycentric>& FlightPlan::last_coast() {
  return *segments_.back();
}

DiscreteTrajectory<Barycentric>& FlightPlan::penultimate_coast() {
  // The penultimate coast is the antepenultimate segment.
  return *segments_[segments_.size() - 3];
}

}  // namespace ksp_plugin
}  // namespace principia
