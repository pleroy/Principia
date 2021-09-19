﻿
#pragma once

#include "physics/discrete_trajectory.hpp"

#include <algorithm>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "astronomy/epoch.hpp"
#include "base/flags.hpp"
#include "base/not_null.hpp"
#include "base/zfp_compressor.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "numerics/fit_hermite_spline.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_forkable {

using geometry::Instant;

template<typename Frame>
Instant const& DiscreteTrajectoryTraits<Frame>::time(
    TimelineConstIterator const it) {
  return it->first;
}

template<typename Frame>
DiscreteTrajectoryIterator<Frame>::reference::reference(
    typename DiscreteTrajectoryTraits<Frame>::TimelineConstIterator it)
    : time(it->first), degrees_of_freedom(it->second) {}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::reference
DiscreteTrajectoryIterator<Frame>::operator*() const {
  auto const& it = this->current();
  return reference(it);
}

template<typename Frame>
std::optional<typename DiscreteTrajectoryIterator<Frame>::reference>
    DiscreteTrajectoryIterator<Frame>::operator->() const {
  auto const& it = this->current();
  return std::make_optional<reference>(it);
}

template<typename Frame>
not_null<DiscreteTrajectoryIterator<Frame>*>
DiscreteTrajectoryIterator<Frame>::that() {
  return this;
}

template<typename Frame>
not_null<DiscreteTrajectoryIterator<Frame> const*>
DiscreteTrajectoryIterator<Frame>::that() const {
  return this;
}

}  // namespace internal_forkable

namespace internal_discrete_trajectory {

using astronomy::InfiniteFuture;
using astronomy::InfinitePast;
using base::Flags;
using base::make_not_null_unique;
using base::ZfpCompressor;
using geometry::Displacement;
using numerics::FitHermiteSpline;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Second;

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkWithCopy(Instant const& time) {
  // May be at |timeline_end()| if |time| is the fork time of this object.
  auto timeline_it = timeline_.find(time);
  CHECK(timeline_it != timeline_end() ||
        (!this->is_root() && time == this->Fork()->time))
      << "NewForkWithCopy at nonexistent time " << time;

  auto const fork = this->NewFork(timeline_it);

  // Copy the tail of the trajectory in the child object.
  if (timeline_it != timeline_.end()) {
    fork->timeline_.insert(++timeline_it, timeline_.end());
  }
  return fork;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkWithoutCopy(Instant const& time) {
  // May be at |timeline_end()| if |time| is the fork time of this object.
  auto timeline_it = timeline_.find(time);
  CHECK(timeline_it != timeline_end() ||
        (!this->is_root() && time == this->Fork()->time))
      << "NewForkWithoutCopy at nonexistent time " << time;

  return this->NewFork(timeline_it);
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkAtLast() {
  auto end = timeline_.end();
  if (timeline_.empty()) {
    return this->NewFork(end);
  } else {
    return this->NewFork(--end);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::AttachFork(
    not_null<std::unique_ptr<DiscreteTrajectory<Frame>>> fork) {
  CHECK(fork->is_root());
  CHECK(!this->Empty());
  // It is easy to mess up forks and to end up with a trajectory which is its
  // own parent.  This check catches that problem.  It's not foolproof (there
  // are more complicated anomalies that can arise) but it's still useful.
  CHECK(fork.get() != this);

  auto& fork_timeline = fork->timeline_;
  auto const this_last = --this->end();

  // Determine if |fork| already has a point matching the end of this
  // trajectory.
  bool must_prepend;
  if (fork_timeline.empty()) {
    must_prepend = true;
  } else {
    CHECK_LE(this_last->time, fork_timeline.begin()->first);
    auto const it = fork_timeline.find(this_last->time);
    if (it == fork_timeline.end()) {
      must_prepend = true;
    } else {
      CHECK(it == fork_timeline.begin())
          << it->first << " " << this_last->time;
      must_prepend = false;
    }
  }

  // If needed, prepend to |fork| a copy of the last point of this trajectory.
  // This ensures that |fork| and this trajectory start and end, respectively,
  // with points at the same time (but possibly distinct degrees of freedom).
  if (must_prepend) {
    fork_timeline.emplace_hint(fork_timeline.begin(),
                               this_last->time,
                               this_last->degrees_of_freedom);
  }

  // Attach |fork| to this trajectory.
  this->AttachForkToCopiedBegin(std::move(fork));

  // Remove the first point of |fork| now that it is properly attached to its
  // parent: that point is either redundant (if it was prepended above) or wrong
  // (because we "trust" this trajectory more than |fork|).  The children that
  // might have been forked at the deleted point were relocated by
  // AttachForkToCopiedBegin.
  fork_timeline.erase(fork_timeline.begin());
}

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>>
DiscreteTrajectory<Frame>::DetachFork() {
  CHECK(!this->is_root());

  // Insert a new point in the timeline for the fork time.  It should go at the
  // beginning of the timeline.
  auto const fork_it = this->Fork();
  auto const begin_it = timeline_.emplace_hint(
      timeline_.begin(), fork_it->time, fork_it->degrees_of_freedom);
  CHECK(begin_it == timeline_.begin());

  // Detach this trajectory and tell the caller that it owns the pieces.
  return this->DetachForkWithCopiedBegin();
}

template<typename Frame>
absl::Status DiscreteTrajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  CHECK(this->is_root() || time > this->Fork()->time)
       << "Append at " << time << " which is before fork time "
       << this->Fork()->time;

  if (!timeline_.empty() && timeline_.cbegin()->first == time) {
    LOG(WARNING) << "Append at existing time " << time
                 << ", time range = [" << this->front().time << ", "
                 << this->back().time << "]";
    return absl::OkStatus();
  }
  auto const it = timeline_.emplace_hint(timeline_.end(),
                                         time,
                                         degrees_of_freedom);
  // Decrementing |end()| is much faster than incrementing |it|.  Don't ask.
  CHECK(--timeline_.end() == it)
      << "Append out of order at " << time << ", last time is "
      << (--timeline_.end())->first;
  if (downsampling_.has_value()) {
    return UpdateDownsampling(it);
  } else {
    return absl::OkStatus();
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetAfter(Instant const& time) {
  this->DeleteAllForksAfter(time);
  if (downsampling_.has_value()) {
    downsampling_->ForgetAfter(time);
  }

  // Get an iterator denoting the first entry with time > |time|.  Remove that
  // entry and all the entries that follow it.  This preserves any entry with
  // time == |time|.
  auto const first_removed_in_timeline = timeline_.upper_bound(time);
  timeline_.erase(first_removed_in_timeline, timeline_.end());

  if (!timeline_.empty() &&
      downsampling_.has_value() &&
      downsampling_->empty()) {
    // Further points will be appended to the last remaining point, so this is
    // where the dense timeline will begin.
    downsampling_->Append(--timeline_.cend());
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetBefore(Instant const& time) {
  CHECK(this->is_root());
  this->CheckNoForksBefore(time);
  if (downsampling_.has_value()) {
    downsampling_->ForgetBefore(time);
  }

  // Get an iterator denoting the first entry with time >= |time|.  Remove all
  // the entries that precede it.  This preserves any entry with time == |time|.
  auto const first_kept_in_timeline = timeline_.lower_bound(time);
  timeline_.erase(timeline_.begin(), first_kept_in_timeline);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::SetDownsampling(
    DownsamplingParameters const& downsampling_parameters) {
  CHECK(!downsampling_.has_value());
  downsampling_.emplace(downsampling_parameters,
                        [this](Instant const& t) {
                          auto const it = this->Find(t);
                          CHECK(it != this->end());
                          return it.current();
                        });

  // For a fork, the fork point is taken into account for downsampling: it is
  // always preserved, but the second point preserved after downsampling may be
  // farther down the timeline of this trajectory.
  if (!this->is_root()) {
    downsampling_->Append(this->Fork().current());
  }
  for (auto it = timeline_.cbegin(); it != timeline_.cend(); ++it) {
    // When reading pre-陈景润 saves we may call this function on long
    // trajectories, in which case we need to downsample as we go.
    UpdateDownsampling(it);
  }
}
template<typename Frame>
void DiscreteTrajectory<Frame>::ClearDownsampling() {
  downsampling_.reset();
}

template<typename Frame>
Instant DiscreteTrajectory<Frame>::t_min() const {
  return this->Empty() ? InfiniteFuture : this->front().time;
}

template<typename Frame>
Instant DiscreteTrajectory<Frame>::t_max() const {
  return this->Empty() ? InfinitePast : this->back().time;
}

template<typename Frame>
Position<Frame> DiscreteTrajectory<Frame>::EvaluatePosition(
    Instant const& time) const {
  auto const iter = this->LowerBound(time);
  if (iter->time == time) {
    return iter->degrees_of_freedom.position();
  }
  CHECK_LT(t_min(), time);
  CHECK_GT(t_max(), time);
  return GetInterpolation(iter).Evaluate(time);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectory<Frame>::EvaluateVelocity(
    Instant const& time) const {
  auto const iter = this->LowerBound(time);
  if (iter->time == time) {
    return iter->degrees_of_freedom.velocity();
  }
  CHECK_LT(t_min(), time);
  CHECK_GT(t_max(), time);
  return GetInterpolation(iter).EvaluateDerivative(time);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectory<Frame>::EvaluateDegreesOfFreedom(
    Instant const& time) const {
  auto const iter = this->LowerBound(time);
  if (iter->time == time) {
    return iter->degrees_of_freedom;
  }
  CHECK_LT(t_min(), time);
  CHECK_GT(t_max(), time);
  auto const interpolation = GetInterpolation(iter);
  return {interpolation.Evaluate(time), interpolation.EvaluateDerivative(time)};
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteToMessage(
    std::set<DiscreteTrajectory const*> const& excluded,
    std::vector<DiscreteTrajectory const*> const& tracked,
    std::vector<Iterator> const& exact,
    not_null<serialization::DiscreteTrajectory*> const message) const {
  WriteToMessage(InfinitePast, excluded, tracked, exact, message);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteToMessage(
    Instant const& after_time,
    std::set<DiscreteTrajectory const*> const& excluded,
    std::vector<DiscreteTrajectory const*> const& tracked,
    std::vector<Iterator> const& exact,
    not_null<serialization::DiscreteTrajectory*> const message) const {
  CHECK(this->is_root());

  std::set<DiscreteTrajectory<Frame> const*> mutable_excluded = excluded;
  std::vector<DiscreteTrajectory<Frame> const*> mutable_tracked = tracked;
  WriteSubTreeToMessage(after_time, mutable_excluded, mutable_tracked, message);
  CHECK(std::all_of(mutable_excluded.begin(),
                    mutable_excluded.end(),
                    [](DiscreteTrajectory<Frame> const* const fork) {
                      return fork == nullptr;
                    }));
  CHECK(std::all_of(mutable_tracked.begin(),
                    mutable_tracked.end(),
                    [](DiscreteTrajectory<Frame> const* const fork) {
                      return fork == nullptr;
                    }));

  for (auto const& it : exact) {
    auto* const serialized_exact = message->add_exact();
    it->time.WriteToMessage(serialized_exact->mutable_instant());
    it->degrees_of_freedom.WriteToMessage(
        serialized_exact->mutable_degrees_of_freedom());
  }
}

template<typename Frame>
template<typename, typename>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>>
DiscreteTrajectory<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<DiscreteTrajectory<Frame>**> const& forks) {
  // We use a Timeline as a convenient container for the exact points.
  Timeline exact;
  for (auto const& instantaneous_degrees_of_freedom : message.exact()) {
    auto const t =
        Instant::ReadFromMessage(instantaneous_degrees_of_freedom.instant());
    auto const degrees_of_freedom = DegreesOfFreedom<Frame>::ReadFromMessage(
        instantaneous_degrees_of_freedom.degrees_of_freedom());
    exact.emplace_hint(exact.end(), t, degrees_of_freedom);
  }

  auto trajectory = make_not_null_unique<DiscreteTrajectory>();
  CHECK(std::all_of(forks.begin(),
                    forks.end(),
                    [](DiscreteTrajectory<Frame>** const fork) {
                      return fork != nullptr && *fork == nullptr;
                    }));
  trajectory->FillSubTreeFromMessage(message, forks, exact);
  return trajectory;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*> DiscreteTrajectory<Frame>::that() {
  return this;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame> const*>
DiscreteTrajectory<Frame>::that() const {
  return this;
}

template<typename Frame>
std::pair<typename DiscreteTrajectory<Frame>::TimelineConstIterator, bool>
DiscreteTrajectory<Frame>::timeline_insert(
    const typename TimelineConstIterator::value_type& value) {
  return timeline_.insert(value);
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_begin() const {
  return timeline_.begin();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_end() const {
  return timeline_.end();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_find(Instant const& time) const {
  return timeline_.find(time);
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_lower_bound(Instant const& time) const {
  return timeline_.lower_bound(time);
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::timeline_empty() const {
  return timeline_.empty();
}

template<typename Frame>
std::int64_t DiscreteTrajectory<Frame>::timeline_size() const {
  return timeline_.size();
}

template<typename Frame>
DiscreteTrajectory<Frame>::Downsampling::Downsampling(
    DownsamplingParameters const& downsampling_parameters,
    std::function<TimelineConstIterator(Instant const&)> iterator_for_time)
    : downsampling_parameters_(downsampling_parameters),
      iterator_for_time_(std::move(iterator_for_time)) {
  // This contains points, hence one more than intervals.
  dense_iterators_.reserve(max_dense_intervals() + 1);
}

template<typename Frame>
std::int64_t DiscreteTrajectory<Frame>::Downsampling::max_dense_intervals()
    const {
  return downsampling_parameters_.max_dense_intervals;
}

template<typename Frame>
Length DiscreteTrajectory<Frame>::Downsampling::tolerance() const {
  return downsampling_parameters_.tolerance;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::Append(
    TimelineConstIterator const it) {
  CHECK(!full());
  if (empty()) {
    start_time_ = it->first;
  }
  dense_iterators_.push_back(it);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::ForgetAfter(Instant const& t) {
  UpdateDenseIteratorsIfNeeded();
  auto const it = std::upper_bound(
      dense_iterators_.cbegin(),
      dense_iterators_.cend(),
      t,
      [](Instant const& left, TimelineConstIterator const right) {
        return left < right->first;
      });
  dense_iterators_.erase(it, dense_iterators_.cend());
  UpdateStartTimeIfNeeded();
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::ForgetBefore(Instant const& t) {
  UpdateDenseIteratorsIfNeeded();
  auto const it = std::lower_bound(
      dense_iterators_.cbegin(),
      dense_iterators_.cend(),
      t,
      [](TimelineConstIterator const left, Instant const& right) {
        return left->first < right;
      });
  dense_iterators_.erase(dense_iterators_.cbegin(), it);
  UpdateStartTimeIfNeeded();
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::Downsampling::empty() const {
  return dense_iterators_.empty();
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::Downsampling::full() const {
  return dense_iterators_.size() >= max_dense_intervals();
}

template<typename Frame>
std::vector<typename DiscreteTrajectory<Frame>::TimelineConstIterator>
DiscreteTrajectory<Frame>::Downsampling::ExtractDenseIterators() {
  UpdateDenseIteratorsIfNeeded();
  auto returned_dense_iterators = std::move(dense_iterators_);
  UpdateStartTimeIfNeeded();
  return std::move(returned_dense_iterators);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::WriteToMessage(
    Instant const& after_time,
    not_null<serialization::DiscreteTrajectory::Downsampling*> const message)
    const {
  message->set_max_dense_intervals(max_dense_intervals());
  tolerance().WriteToMessage(message->mutable_tolerance());
  if (!empty()) {
    // We can't call UpdateDenseIteratorsIfNeeded here because this method is
    // const.
    start_time_.value().WriteToMessage(message->add_dense_timeline());
    for (int i = 1; i < dense_iterators_.size(); ++i) {
      Instant const& time = dense_iterators_[i]->first;
      if (time >= after_time) {
        time.WriteToMessage(message->add_dense_timeline());
      }
    }
  }
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::Downsampling
DiscreteTrajectory<Frame>::Downsampling::ReadFromMessage(
    serialization::DiscreteTrajectory::Downsampling const& message,
    DiscreteTrajectory const& trajectory) {
  bool const is_pre_zermelo = message.has_start_of_dense_timeline();
  LOG_IF(WARNING, is_pre_zermelo)
      << "Reading pre-Zermelo Downsampling";
  Downsampling downsampling({message.max_dense_intervals(),
                             Length::ReadFromMessage(message.tolerance())},
                            [&trajectory](Instant const& t){
                              auto const it = trajectory.Find(t);
                              CHECK(it != trajectory.end());
                              return it.current();
                            });
  if (is_pre_zermelo) {
    // No support for forks in legacy saves, so |find| will succeed and ++ is
    // safe.
    auto it = trajectory.timeline_.find(
        Instant::ReadFromMessage(message.start_of_dense_timeline()));
    CHECK(it != trajectory.timeline_.end());
    for (; it != trajectory.timeline_.end(); ++it) {
      downsampling.Append(it);
    }
  } else {
    for (auto const& dense_time : message.dense_timeline()) {
      auto const t = Instant::ReadFromMessage(dense_time);
      downsampling.Append(downsampling.iterator_for_time_(t));
    }
  }
  return downsampling;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::UpdateDenseIteratorsIfNeeded() {
  if (!dense_iterators_.empty()) {
    dense_iterators_[0] = iterator_for_time_(start_time_.value());
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::UpdateStartTimeIfNeeded() {
  if (dense_iterators_.empty()) {
    start_time_.reset();
  }
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::WriteSubTreeToMessage(
    Instant const& after_time,
    std::set<DiscreteTrajectory const*>& excluded,
    std::vector<DiscreteTrajectory const*>& tracked,
    not_null<serialization::DiscreteTrajectory*> const message) const {
  bool const included =
      Forkable<DiscreteTrajectory, Iterator, DiscreteTrajectoryTraits<Frame>>::
          WriteSubTreeToMessage(after_time, excluded, tracked, message);
  if (!included) {
    return false;
  }

  int timeline_size = timeline_.size();
  auto* const zfp = message->mutable_zfp();

  // The timeline data is made dimensionless and stored in separate arrays per
  // coordinate.  We expect strong correlations within a coordinate over time,
  // but not between coordinates.
  std::vector<double> t;
  std::vector<double> qx;
  std::vector<double> qy;
  std::vector<double> qz;
  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> pz;
  t.reserve(timeline_size);
  qx.reserve(timeline_size);
  qy.reserve(timeline_size);
  qz.reserve(timeline_size);
  px.reserve(timeline_size);
  py.reserve(timeline_size);
  pz.reserve(timeline_size);
  std::optional<Instant> previous_instant;
  Time max_Δt;
  std::string* const zfp_timeline = zfp->mutable_timeline();
  for (auto const& [time, degrees_of_freedom] : timeline_) {
    if (time < after_time) {
      --timeline_size;
    } else {
      auto const q = degrees_of_freedom.position() - Frame::origin;
      auto const p = degrees_of_freedom.velocity();
      t.push_back((time - Instant{}) / Second);
      qx.push_back(q.coordinates().x / Metre);
      qy.push_back(q.coordinates().y / Metre);
      qz.push_back(q.coordinates().z / Metre);
      px.push_back(p.coordinates().x / (Metre / Second));
      py.push_back(p.coordinates().y / (Metre / Second));
      pz.push_back(p.coordinates().z / (Metre / Second));
      if (previous_instant.has_value()) {
        max_Δt = std::max(max_Δt, time - *previous_instant);
      }
      previous_instant = time;
    }
  }
  zfp->set_timeline_size(timeline_size);

  // Times are exact.
  ZfpCompressor time_compressor(0);
  // Lengths are approximated to the downsampling tolerance if downsampling is
  // enabled, otherwise they are exact.
  Length const length_tolerance =
      downsampling_.has_value() ? downsampling_->tolerance() : Length();
  ZfpCompressor length_compressor(length_tolerance / Metre);
  // Speeds are approximated based on the length tolerance and the maximum
  // step in the timeline.
  ZfpCompressor const speed_compressor((length_tolerance / max_Δt) /
                                        (Metre / Second));

  ZfpCompressor::WriteVersion(message);
  time_compressor.WriteToMessageMultidimensional<2>(t, zfp_timeline);
  length_compressor.WriteToMessageMultidimensional<2>(qx, zfp_timeline);
  length_compressor.WriteToMessageMultidimensional<2>(qy, zfp_timeline);
  length_compressor.WriteToMessageMultidimensional<2>(qz, zfp_timeline);
  speed_compressor.WriteToMessageMultidimensional<2>(px, zfp_timeline);
  speed_compressor.WriteToMessageMultidimensional<2>(py, zfp_timeline);
  speed_compressor.WriteToMessageMultidimensional<2>(pz, zfp_timeline);

  if (downsampling_.has_value()) {
    downsampling_->WriteToMessage(after_time, message->mutable_downsampling());
  }
  return true;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::FillSubTreeFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<DiscreteTrajectory<Frame>**> const& tracked,
    Timeline const& exact) {
  bool const is_pre_frobenius = !message.has_zfp();
  LOG_IF(WARNING, is_pre_frobenius)
      << "Reading pre-Frobenius DiscreteTrajectory";
  if (is_pre_frobenius) {
    for (auto const& instantaneous_dof : message.timeline()) {
      Append(Instant::ReadFromMessage(instantaneous_dof.instant()),
             DegreesOfFreedom<Frame>::ReadFromMessage(
                 instantaneous_dof.degrees_of_freedom()));
    }
  } else {
    ZfpCompressor decompressor;
    ZfpCompressor::ReadVersion(message);

    int const timeline_size = message.zfp().timeline_size();
    std::vector<double> t(timeline_size);
    std::vector<double> qx(timeline_size);
    std::vector<double> qy(timeline_size);
    std::vector<double> qz(timeline_size);
    std::vector<double> px(timeline_size);
    std::vector<double> py(timeline_size);
    std::vector<double> pz(timeline_size);
    std::string_view zfp_timeline(message.zfp().timeline().data(),
                                  message.zfp().timeline().size());

    decompressor.ReadFromMessageMultidimensional<2>(t, zfp_timeline);
    decompressor.ReadFromMessageMultidimensional<2>(qx, zfp_timeline);
    decompressor.ReadFromMessageMultidimensional<2>(qy, zfp_timeline);
    decompressor.ReadFromMessageMultidimensional<2>(qz, zfp_timeline);
    decompressor.ReadFromMessageMultidimensional<2>(px, zfp_timeline);
    decompressor.ReadFromMessageMultidimensional<2>(py, zfp_timeline);
    decompressor.ReadFromMessageMultidimensional<2>(pz, zfp_timeline);

    for (int i = 0; i < timeline_size; ++i) {
      Position<Frame> const q =
          Frame::origin +
          Displacement<Frame>({qx[i] * Metre, qy[i] * Metre, qz[i] * Metre});
      Velocity<Frame> const p({px[i] * (Metre / Second),
                               py[i] * (Metre / Second),
                               pz[i] * (Metre / Second)});

      // See if this is a point whose degrees of freedom must be restored
      // exactly.
      Instant const time = Instant() + t[i] * Second;
      if (auto it = exact.find(time); it == exact.cend()) {
        Append(time, DegreesOfFreedom<Frame>(q, p));
      } else {
        Append(time, it->second);
      }
    }
  }
  if (message.has_downsampling()) {
    downsampling_.emplace(
        Downsampling::ReadFromMessage(message.downsampling(), *this));
  }
  Forkable<DiscreteTrajectory, Iterator, DiscreteTrajectoryTraits<Frame>>::
      FillSubTreeFromMessage(message, tracked, exact);
}

template<typename Frame>
Hermite3<Instant, Position<Frame>> DiscreteTrajectory<Frame>::GetInterpolation(
    Iterator const& upper) const {
  CHECK(upper != this->begin());
  auto const lower = --Iterator{upper};
  return Hermite3<Instant, Position<Frame>>{
      {lower->time, upper->time},
      {lower->degrees_of_freedom.position(),
       upper->degrees_of_freedom.position()},
      {lower->degrees_of_freedom.velocity(),
       upper->degrees_of_freedom.velocity()}};
}

template<typename Frame>
absl::Status DiscreteTrajectory<Frame>::UpdateDownsampling(
    TimelineConstIterator const appended) {
  this->CheckNoForksBefore(appended->first);
  downsampling_->Append(appended);
  if (downsampling_->full()) {
    auto const dense_iterators = downsampling_->ExtractDenseIterators();
    auto right_endpoints = FitHermiteSpline<Instant, Position<Frame>>(
        dense_iterators,
        [](auto&& it) -> auto&& { return it->first; },
        [](auto&& it) -> auto&& { return it->second.position(); },
        [](auto&& it) -> auto&& { return it->second.velocity(); },
        downsampling_->tolerance());
    if (!right_endpoints.ok()) {
      // Note that the actual appending took place; the propagated status only
      // reflects a lack of downsampling.
      return right_endpoints.status();
    }
    if (right_endpoints->empty()) {
      right_endpoints->push_back(dense_iterators.cend() - 1);
    }

    // Poke holes in the timeline at the places given by |right_endpoints|.
    // Note that this code carefully avoids incrementing |dense_iterators[0]| as
    // it may live in a different timeline than the others.
    CHECK_LE(1, dense_iterators.size());
    TimelineConstIterator left = dense_iterators[1];
    for (const auto& it_in_dense_iterators : right_endpoints.value()) {
      TimelineConstIterator const right = *it_in_dense_iterators;
      timeline_.erase(left, right);
      left = right;
      ++left;
    }

    // Re-append the dense iterators that have not been consumed.
    for (auto it = right_endpoints->back(); it < dense_iterators.cend(); ++it) {
      downsampling_->Append(*it);
    }
    CHECK(!downsampling_->empty());
    CHECK(!downsampling_->full());
  }
  return absl::OkStatus();
}

}  // namespace internal_discrete_trajectory
}  // namespace physics
}  // namespace principia
