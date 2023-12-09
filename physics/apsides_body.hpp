#pragma once

#include "physics/apsides.hpp"

#include <optional>

#include "base/array.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "geometry/space.hpp"
#include "numerics/hermite3.hpp"
#include "numerics/root_finders.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace _apsides {
namespace internal {

using namespace principia::base::_array;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_sign;
using namespace principia::geometry::_space;
using namespace principia::numerics::_hermite3;
using namespace principia::numerics::_root_finders;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;

template<typename Frame>
void ComputeApsides(Trajectory<Frame> const& reference,
                    Trajectory<Frame> const& trajectory,
                    typename DiscreteTrajectory<Frame>::iterator const begin,
                    typename DiscreteTrajectory<Frame>::iterator const end,
                    int const max_points,
                    DiscreteTrajectory<Frame>& apoapsides,
                    DiscreteTrajectory<Frame>& periapsides) {
  std::optional<Instant> previous_time;
  std::optional<DegreesOfFreedom<Frame>> previous_degrees_of_freedom;
  std::optional<Square<Length>> previous_squared_distance;
  std::optional<Variation<Square<Length>>>
      previous_squared_distance_derivative;

  Instant const t_min = reference.t_min();
  Instant const t_max = reference.t_max();
  for (auto it = begin; it != end; ++it) {
    auto const& [time, degrees_of_freedom] = *it;
    if (time < t_min) {
      continue;
    }
    if (time > t_max) {
      break;
    }
    DegreesOfFreedom<Frame> const body_degrees_of_freedom =
        reference.EvaluateDegreesOfFreedom(time);
    RelativeDegreesOfFreedom<Frame> const relative =
        degrees_of_freedom - body_degrees_of_freedom;
    Square<Length> const squared_distance = relative.displacement().Norm²();
    // This is the derivative of |squared_distance|.
    Variation<Square<Length>> const squared_distance_derivative =
        2.0 * InnerProduct(relative.displacement(), relative.velocity());

    if (previous_squared_distance_derivative &&
        Sign(squared_distance_derivative) !=
            Sign(*previous_squared_distance_derivative)) {
      CHECK(previous_time &&
            previous_degrees_of_freedom &&
            previous_squared_distance);

      // The derivative of |squared_distance| changed sign.  Construct a Hermite
      // approximation of |squared_distance| and find its extrema.
      Hermite3<Instant, Square<Length>> const
          squared_distance_approximation(
              {*previous_time, time},
              {*previous_squared_distance, squared_distance},
              {*previous_squared_distance_derivative,
               squared_distance_derivative});
      BoundedArray<Instant, 2> const extrema =
          squared_distance_approximation.FindExtrema();

      // Now look at the extrema and check that exactly one is in the required
      // time interval.  This is normally the case, but it can fail due to
      // ill-conditioning.
      Instant apsis_time;
      int valid_extrema = 0;
      for (auto const& extremum : extrema) {
        if (extremum >= *previous_time && extremum <= time) {
          apsis_time = extremum;
          ++valid_extrema;
        }
      }
      if (valid_extrema != 1) {
        // Something went wrong when finding the extrema of
        // |squared_distance_approximation|. Use a linear interpolation of
        // |squared_distance_derivative| instead.
        apsis_time = Barycentre<Instant, Variation<Square<Length>>>(
            {time, *previous_time},
            {*previous_squared_distance_derivative,
             -squared_distance_derivative});
      }

      // This can happen for instance if the square distance is stationary.
      // Safer to give up.
      if (!IsFinite(apsis_time - Instant{})) {
        break;
      }

      // Now that we know the time of the apsis, use a Hermite approximation to
      // derive its degrees of freedom.  Note that an extremum of
      // |squared_distance_approximation| is in general not an extremum for
      // |position_approximation|: the distance computed using the latter is a
      // 6th-degree polynomial.  However, approximating this polynomial using a
      // 3rd-degree polynomial would yield |squared_distance_approximation|, so
      // we shouldn't be far from the truth.
      DegreesOfFreedom<Frame> const apsis_degrees_of_freedom =
          trajectory.EvaluateDegreesOfFreedom(apsis_time);
      if (Sign(squared_distance_derivative).is_negative()) {
        apoapsides.Append(apsis_time, apsis_degrees_of_freedom).IgnoreError();
      } else {
        periapsides.Append(apsis_time, apsis_degrees_of_freedom).IgnoreError();
      }
      if (apoapsides.size() >= max_points && periapsides.size() >= max_points) {
        break;
      }
    }

    previous_time = time;
    previous_degrees_of_freedom = degrees_of_freedom;
    previous_squared_distance = squared_distance;
    previous_squared_distance_derivative = squared_distance_derivative;
  }
}

template<typename Frame>
std::optional<typename DiscreteTrajectory<Frame>::value_type> ComputeCollision(
    RotatingBody<Frame> const& reference_body,
    Trajectory<Frame> const& reference,
    Trajectory<Frame> const& trajectory,
    typename DiscreteTrajectory<Frame>::iterator begin,
    typename DiscreteTrajectory<Frame>::iterator end,
    std::function<Length(Angle const& latitude,
                         Angle const& longitude)> const& radius) {
  // The frame of the surface of the celestial.
  using SurfaceFrame = geometry::_frame::Frame<struct SurfaceFrameTag>;

  auto squared_distance_from_centre =
      [&reference](
          typename DiscreteTrajectory<Frame>::iterator const it) {
        Position<Frame> const body_position =
            reference.EvaluatePosition(it->time);
        Position<Frame> const trajectory_position =
            it->degrees_of_freedom.position();
        Displacement<Frame> const displacement =
            trajectory_position - body_position;
        return displacement.Norm²();
      };

  auto const max_radius² = Pow<2>(reference_body.max_radius());
  auto const min_radius² = Pow<2>(reference_body.min_radius());

  //TODO(phl): Loop over all the segments.

  // Determine if the |trajectory| crosses below |max_radius|, and if so find
  // the first segment where it is below |max_radius|.  This is local to C++.
  typename DiscreteTrajectory<Frame>::iterator below_max_radius_begin = begin;
  typename DiscreteTrajectory<Frame>::iterator below_max_radius_end = end;
  bool has_point_below_max_radius = false;
  Instant const t_min = reference.t_min();
  Instant const t_max = reference.t_max();
  for (auto it = begin; it != end; ++it) {
    auto const& [time, _] = *it;
    if (time < t_min) {
      below_max_radius_begin = it;
      continue;
    }
    if (time > t_max) {
      below_max_radius_end = it;
      break;
    }
    if (squared_distance_from_centre(it) <= max_radius²) {
      if (!has_point_below_max_radius) {
        below_max_radius_begin = it;
        has_point_below_max_radius = true;
      }
    } else if (has_point_below_max_radius) {
      below_max_radius_end = it;
      break;
    }
  }

  // If there is no point below max_radius surely there is no collision.
  if (!has_point_below_max_radius) {
    return std::nullopt;
  }

  auto height_above_terrain_at_time = [&radius,
                                       &reference,
                                       &reference_body,
                                       &trajectory](Instant const& t) {
    auto const reference_position = reference.EvaluatePosition(t);
    auto const trajectory_position = trajectory.EvaluatePosition(t);
    Displacement<Frame> const displacement_in_frame =
        trajectory_position - reference_position;

    auto const to_surface_frame =
        reference_body.template ToSurfaceFrame<SurfaceFrame>(t);
    Displacement<SurfaceFrame> const displacement_in_surface =
        to_surface_frame(displacement_in_frame);

    SphericalCoordinates<Length> const spherical_coordinates =
        displacement_in_surface.coordinates().ToSpherical();

    return spherical_coordinates.radius -
           radius(spherical_coordinates.latitude,
                  spherical_coordinates.longitude);
  };

  // Determine if there is a point where the |trajectory| is below the terrain.
  // This may be expensive as it calls into C#.
  for (auto it = below_max_radius_begin; it != below_max_radius_end; ++it) {
    //TODO(phl): Optimize using mix_radius?
    if (height_above_terrain_at_time(it->time) < Length{}) {
      // The trajectory is below the terrain.  Find the exact crossing point.
      CHECK(it != begin) << "Below terrain at " << begin->time << ": "
                         << begin->degrees_of_freedom;
      auto const previous_it = std::prev(it);
      Instant const collision_time =
          Brent(height_above_terrain_at_time, previous_it->time, it->time);
      return typename DiscreteTrajectory<Frame>::value_type(
          collision_time,
          trajectory.EvaluateDegreesOfFreedom(collision_time));
    }
  }

  return std::nullopt;
}

template<typename Frame, typename Predicate>
absl::Status ComputeNodes(
    Trajectory<Frame> const& trajectory,
    typename DiscreteTrajectory<Frame>::iterator const begin,
    typename DiscreteTrajectory<Frame>::iterator const end,
    Vector<double, Frame> const& north,
    int const max_points,
    DiscreteTrajectory<Frame>& ascending,
    DiscreteTrajectory<Frame>& descending,
    Predicate predicate) {
  static_assert(
      std::is_convertible<decltype(predicate(
                              std::declval<DegreesOfFreedom<Frame>>())),
                          bool>::value,
      "|predicate| must be a predicate on |DegreesOfFreedom<Frame>|");

  std::optional<Instant> previous_time;
  std::optional<Length> previous_z;
  std::optional<Speed> previous_z_speed;

  for (auto it = begin; it != end; ++it) {
    RETURN_IF_STOPPED;
    auto const& [time, degrees_of_freedom] = *it;
    Length const z =
        (degrees_of_freedom.position() - Frame::origin).coordinates().z;
    Speed const z_speed = degrees_of_freedom.velocity().coordinates().z;

    if (previous_z && Sign(z) != Sign(*previous_z)) {
      CHECK(previous_time && previous_z_speed);

      // |z| changed sign.  Construct a Hermite approximation of |z| and find
      // its zeros.
      Hermite3<Instant, Length> const z_approximation(
          {*previous_time, time},
          {*previous_z, z},
          {*previous_z_speed, z_speed});

      Instant node_time;
      if (Sign(z_approximation.Evaluate(*previous_time)) ==
          Sign(z_approximation.Evaluate(time))) {
        // The Hermite approximation is poorly conditioned, let's use a linear
        // approximation
        node_time = Barycentre<Instant, Length>({*previous_time, time},
                                                {z, -*previous_z});
      } else {
        // The normal case, find the intersection with z = 0 using Brent's
        // method.
        node_time = Brent(
            [&z_approximation](Instant const& t) {
              return z_approximation.Evaluate(t);
            },
            *previous_time,
            time);
      }

      DegreesOfFreedom<Frame> const node_degrees_of_freedom =
          trajectory.EvaluateDegreesOfFreedom(node_time);
      if (predicate(node_degrees_of_freedom)) {
        if (Sign(InnerProduct(north, Vector<double, Frame>({0, 0, 1}))) ==
            Sign(z_speed)) {
          // |north| is up and we are going up, or |north| is down and we are
          // going down.
          ascending.Append(node_time, node_degrees_of_freedom).IgnoreError();
        } else {
          descending.Append(node_time, node_degrees_of_freedom).IgnoreError();
        }
        if (ascending.size() >= max_points && descending.size() >= max_points) {
          break;
        }
      }
    }

    previous_time = time;
    previous_z = z;
    previous_z_speed = z_speed;
  }
  return absl::OkStatus();
}

}  // namespace internal
}  // namespace _apsides
}  // namespace physics
}  // namespace principia
