#pragma once

#include "physics/apsides.hpp"

#include <vector>

#include "base/array.hpp"
#include "numerics/root_finders.hpp"

namespace principia {
namespace physics {
namespace internal_apsides {

using base::BoundedArray;
using geometry::Barycentre;
using geometry::Instant;
using geometry::Position;
using geometry::Sign;
using numerics::Bisect;
using numerics::Hermite3;
using quantities::Length;
using quantities::Speed;
using quantities::Square;
using quantities::Variation;

template<typename Frame>
void ComputeApsides(Trajectory<Frame> const& reference,
                    typename DiscreteTrajectory<Frame>::Iterator const begin,
                    typename DiscreteTrajectory<Frame>::Iterator const end,
                    DiscreteTrajectory<Frame>& apoapsides,
                    DiscreteTrajectory<Frame>& periapsides) {
  std::experimental::optional<Instant> previous_time;
  std::experimental::optional<DegreesOfFreedom<Frame>>
      previous_degrees_of_freedom;
  std::experimental::optional<Square<Length>> previous_squared_distance;
  std::experimental::optional<Variation<Square<Length>>>
      previous_squared_distance_derivative;

  Instant const t_min = reference.t_min();
  Instant const t_max = reference.t_max();
  for (auto it = begin; it != end; ++it) {
    Instant const& time = it.time();
    if (time < t_min) {
      continue;
    }
    if (time > t_max) {
      break;
    }
    DegreesOfFreedom<Frame> const degrees_of_freedom = it.degrees_of_freedom();
    DegreesOfFreedom<Frame> const body_degrees_of_freedom =
        reference.EvaluateDegreesOfFreedom(time);
    RelativeDegreesOfFreedom<Frame> const relative =
        degrees_of_freedom - body_degrees_of_freedom;
    Square<Length> const squared_distance =
        InnerProduct(relative.displacement(), relative.displacement());
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

      // Now that we know the time of the apsis, use a Hermite approximation to
      // derive its degrees of freedom.  Note that an extremum of
      // |squared_distance_approximation| is in general not an extremum for
      // |position_approximation|: the distance computed using the latter is a
      // 6th-degree polynomial.  However, approximating this polynomial using a
      // 3rd-degree polynomial would yield |squared_distance_approximation|, so
      // we shouldn't be far from the truth.
      DegreesOfFreedom<Frame> const apsis_degrees_of_freedom =
          begin.trajectory()->EvaluateDegreesOfFreedom(apsis_time);
      if (Sign(squared_distance_derivative).Negative()) {
        apoapsides.Append(apsis_time, apsis_degrees_of_freedom);
      } else {
        periapsides.Append(apsis_time, apsis_degrees_of_freedom);
      }
    }

    previous_time = time;
    previous_degrees_of_freedom = degrees_of_freedom;
    previous_squared_distance = squared_distance;
    previous_squared_distance_derivative = squared_distance_derivative;
  }
}

template<typename Frame>
void ComputeNodes(typename DiscreteTrajectory<Frame>::Iterator begin,
                  typename DiscreteTrajectory<Frame>::Iterator end,
                  Vector<double, Frame> const& north,
                  DiscreteTrajectory<Frame>& ascending,
                  DiscreteTrajectory<Frame>& descending) {
  std::experimental::optional<Instant> previous_time;
  std::experimental::optional<Length> previous_z;
  std::experimental::optional<Speed> previous_z_speed;

  for (auto it = begin; it != end; ++it) {
    Instant const time = it.time();
    DegreesOfFreedom<Frame> const& degrees_of_freedom = it.degrees_of_freedom();
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
        // The normal case, find the intersection with z = 0 using bisection.
        // TODO(egg): Bisection on a polynomial seems daft; we should have
        // Newton's method.
        node_time = Bisect(
            [&z_approximation](Instant const& t) {
              return z_approximation.Evaluate(t);
            },
            *previous_time,
            time);
      }

      DegreesOfFreedom<Frame> const node_degrees_of_freedom =
          begin.trajectory()->EvaluateDegreesOfFreedom(node_time);
      if (Sign(InnerProduct(north, Vector<double, Frame>({0, 0, 1}))) ==
          Sign(z_speed)) {
        // |north| is up and we are going up, or |north| is down and we are
        // going down.
        ascending.Append(node_time, node_degrees_of_freedom);
      } else {
        descending.Append(node_time, node_degrees_of_freedom);
      }
    }

    previous_time = time;
    previous_z = z;
    previous_z_speed = z_speed;
  }
}

}  // namespace internal_apsides
}  // namespace physics
}  // namespace principia
