#pragma once

#include <cstdint>
#include <iterator>
#include <optional>
#include <vector>

#include "absl/container/btree_map.h"
#include "absl/status/status.h"
#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "numerics/hermite3.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_iterator.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"
#include "physics/trajectory.hpp"
#include "serialization/physics.pb.h"

namespace principia {

namespace testing_utilities {
FORWARD_DECLARE(TEMPLATE(typename Frame) class,
                DiscreteTrajectoryFactoriesFriend,
                FROM(discrete_trajectory_factories));
}  // namespace testing_utilities

namespace physics {

FORWARD_DECLARE(TEMPLATE(typename Frame) class,
                DiscreteTrajectory,
                FROM(discrete_trajectory),
                INTO(discrete_trajectory_segment));

class DiscreteTrajectoryIteratorTest;
class DiscreteTrajectorySegmentIteratorTest;
class DiscreteTrajectorySegmentTest;

namespace _discrete_trajectory_segment {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::numerics::_hermite3;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory_iterator;
using namespace principia::physics::_discrete_trajectory_segment_iterator;
using namespace principia::physics::_discrete_trajectory_types;
using namespace principia::physics::_trajectory;

template<typename Frame>
class DiscreteTrajectorySegment : public Trajectory<Frame> {
  using Timeline = _discrete_trajectory_types::Timeline<Frame>;

 public:
  using key_type = typename Timeline::key_type;
  using value_type = typename Timeline::value_type;

  using iterator = DiscreteTrajectoryIterator<Frame>;
  using reference = value_type const&;
  using reverse_iterator = std::reverse_iterator<iterator>;

  using DownsamplingParameters =
      _discrete_trajectory_types::DownsamplingParameters;

  // TODO(phl): Decide which constructors should be public.
  DiscreteTrajectorySegment() = default;
  explicit DiscreteTrajectorySegment(
      DiscreteTrajectorySegmentIterator<Frame> self);

  ~DiscreteTrajectorySegment() = default;

  // Moveable.
  DiscreteTrajectorySegment(DiscreteTrajectorySegment&&) = default;
  DiscreteTrajectorySegment& operator=(DiscreteTrajectorySegment&&) = default;
  DiscreteTrajectorySegment(const DiscreteTrajectorySegment&) = delete;
  DiscreteTrajectorySegment& operator=(const DiscreteTrajectorySegment&) =
      delete;

  reference front() const;
  reference back() const;

  iterator begin() const;
  iterator end() const;

  reverse_iterator rbegin() const;
  reverse_iterator rend() const;

  // TODO(phl): We probably don't want empty segments.
  bool empty() const;
  std::int64_t size() const;

  void clear();

  iterator find(Instant const& t) const;

  iterator lower_bound(Instant const& t) const;
  iterator upper_bound(Instant const& t) const;

  Instant t_min() const override;
  Instant t_max() const override;

  Position<Frame> EvaluatePosition(Instant const& t) const override;
  Velocity<Frame> EvaluateVelocity(Instant const& t) const override;
  DegreesOfFreedom<Frame> EvaluateDegreesOfFreedom(
      Instant const& t) const override;

  // This segment must have 0 or 1 points.  Occasionally removes intermediate
  // points from the segment when |Append|ing, ensuring that positions remain
  // within the desired tolerance.
  void SetDownsampling(DownsamplingParameters const& downsampling_parameters);

  // Same as above, but can be used on a nontrivial segment.  Use exclusively
  // for compatibility deserialization.
  void SetDownsamplingUnconditionally(
      DownsamplingParameters const& downsampling_parameters);

  // Clear the downsampling parameters.  From now on, all points appended to the
  // segment are going to be retained.
  void ClearDownsampling();

  // Returns true iff this segment was downsampled at least once since its
  // creation or the last call to |clear|.  Only use for optimization purposes,
  // not to depend on the actual structure of the timeline.
  bool was_downsampled() const;

  // The points denoted by |exact| are written and re-read exactly and are not
  // affected by any errors introduced by zfp compression.  The endpoints of a
  // segment are always exact.
  void WriteToMessage(
      not_null<serialization::DiscreteTrajectorySegment*> message,
      std::vector<iterator> const& exact) const;
  // Same as above, but only the points defined by [begin, end[ are written.
  // The iterators must have been obtained by operations on this segment.
  void WriteToMessage(
      not_null<serialization::DiscreteTrajectorySegment*> message,
      iterator begin,
      iterator end,
      std::vector<iterator> const& exact) const;

  template<typename F = Frame,
           typename = std::enable_if_t<is_serializable_v<F>>>
  static DiscreteTrajectorySegment ReadFromMessage(
      serialization::DiscreteTrajectorySegment const& message,
      DiscreteTrajectorySegmentIterator<Frame> self);

 private:
  // Versions of find, lower_bound, and upper_bound that use optionals to
  // indicate absense. Useful for bypassing expensive end() call.
  std::optional<iterator> FindOrNullopt(Instant const& t) const;
  std::optional<iterator> LowerBoundOrNullopt(Instant const& t) const;
  std::optional<iterator> UpperBoundOrNullopt(Instant const& t) const;

  // Changes the |self_| iterator.  Only for use when attaching/detaching
  // segments.
  void SetSelf(DiscreteTrajectorySegmentIterator<Frame> self);

  void Prepend(Instant const& t,
               DegreesOfFreedom<Frame> const& degrees_of_freedom);

  // Removes all points with a time greater than or equal to |t| (1st overload)
  // or starting at |begin| (2nd overload).
  void ForgetAfter(Instant const& t);
  void ForgetAfter(typename Timeline::const_iterator begin);

  // Removes all points with a time strictly less than |t| (1st overload) or
  // ending at |end| (2nd overload).
  void ForgetBefore(Instant const& t);
  void ForgetBefore(typename Timeline::const_iterator end);

  absl::Status Append(Instant const& t,
                      DegreesOfFreedom<Frame> const& degrees_of_freedom);

  // Merges the points from the given |segment| into this object.  The two
  // segments must have nonoverlapping times.  The downsampling state of the
  // result is that of the latest segment (with the largest times).
  void Merge(DiscreteTrajectorySegment<Frame> segment);

  // Computes |number_of_dense_points_| based on the start of the dense
  // timeline.  Used for compatibility deserialization.
  void SetStartOfDenseTimeline(Instant const& t);

  // Inserts the given point at the beginning of the timeline.  Used for
  // compatibility deserialization.
  void SetForkPoint(value_type const& point);

  // Called by |Append| after appending a point to this segment.  If
  // appropriate, performs downsampling and deletes some of the points of the
  // segment.
  absl::Status DownsampleIfNeeded();

  // Returns the Hermite interpolation for the left-open, right-closed
  // trajectory segment bounded above by |upper|.
  Hermite3<Instant, Position<Frame>> GetInterpolation(
      typename Timeline::const_iterator upper) const;

  typename Timeline::const_iterator timeline_begin() const;
  typename Timeline::const_iterator timeline_end() const;
  bool timeline_empty() const;
  std::int64_t timeline_size() const;

  // Implementation of serialization.  The caller is expected to pass consistent
  // parameters.  |timeline_begin| and |timeline_end| define the range to write.
  // |timeline_size| is the distance from |timeline_begin| to |timeline_end|.
  // |number_of_points_to_skip_at_end| is the distance between |timeline_end|
  // and the true |end| of the timeline.
  void WriteToMessage(
      not_null<serialization::DiscreteTrajectorySegment*> message,
      typename Timeline::const_iterator timeline_begin,
      typename Timeline::const_iterator timeline_end,
      std::int64_t timeline_size,
      std::int64_t number_of_points_to_skip_at_end,
      std::vector<iterator> const& exact) const;

  std::optional<DownsamplingParameters> downsampling_parameters_;

  // The number of points at the end of the segment that are part of a "dense"
  // span, i.e., have not been subjected to downsampling yet.
  std::int64_t number_of_dense_points_ = 0;

  bool was_downsampled_ = false;

  DiscreteTrajectorySegmentIterator<Frame> self_;
  Timeline timeline_;

  template<typename F>
  friend class _discrete_trajectory::internal::DiscreteTrajectory;
  template<typename F>
  friend class _discrete_trajectory_iterator::internal::
      DiscreteTrajectoryIterator;

  // For testing.
  friend class physics::DiscreteTrajectoryIteratorTest;
  friend class physics::DiscreteTrajectorySegmentIteratorTest;
  friend class physics::DiscreteTrajectorySegmentTest;
  template<typename F>
  friend class testing_utilities::_discrete_trajectory_factories::
      DiscreteTrajectoryFactoriesFriend;
};

}  // namespace internal

using internal::DiscreteTrajectorySegment;

}  // namespace _discrete_trajectory_segment
}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_segment_body.hpp"
