﻿#pragma once

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "physics/oblate_body.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace physics {
namespace internal_geopotential {

using base::not_null;
using geometry::Displacement;
using geometry::Instant;
using geometry::Rotation;
using geometry::Vector;
using quantities::Acceleration;
using quantities::Angle;
using quantities::Exponentiation;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Quotient;
using quantities::Square;
using serialization::Numerics;

// Representation of the geopotential model of an oblate body.
template<typename Frame>
class Geopotential {
 public:
  Geopotential(not_null<OblateBody<Frame> const*> body,
               serialization::Numerics::Mode mode);

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  SphericalHarmonicsAcceleration(
      Instant const& t,
      Displacement<Frame> const& r,
      Square<Length> const& r²,
      Exponentiation<Length, -3> const& one_over_r³) const;

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  GeneralSphericalHarmonicsAcceleration(
      Instant const& t,
      Displacement<Frame> const& r,
      Length const& r_norm,
      Square<Length> const& r²,
      Exponentiation<Length, -3> const& one_over_r³) const;

 private:
  // The frame of the surface of the celestial.
  struct SurfaceFrame;
  static const Vector<double, SurfaceFrame> x_;
  static const Vector<double, SurfaceFrame> y_;

  // This is the type that we return, so better have a name for it.
  using ReducedAcceleration = Quotient<Acceleration, GravitationalParameter>;

  // List of reduced accelerations computed for all degrees or orders.
  template<int size>
  using ReducedAccelerations =
      std::array<Vector<ReducedAcceleration, Frame>, size>;

  using UnitVector = Vector<double, Frame>;

  // Holds precomputed data for one evaluation of the acceleration.
  template<int size>
  struct Precomputations;

  // Helper templates for iterating over the degrees/orders of the geopotential.
  template<int size, int degree, int order>
  struct DegreeNOrderM;
  template<int size, int degree, typename>
  struct DegreeNAllOrders;
  template<typename>
  struct AllDegrees;

  // If z is a unit vector along the axis of rotation, and r a vector from the
  // center of |body_| to some point in space, the acceleration computed here
  // is:
  //
  //   -(J₂ / (μ ‖r‖⁵)) (3 z (r.z) + r (3 - 15 (r.z)² / ‖r‖²) / 2)
  //
  // Where ‖r‖ is the norm of r and r.z is the inner product.  It is the
  // additional acceleration exerted by the oblateness of |body| on a point at
  // position r.  J₂, J̃₂ and J̄₂ are normally positive and C̃₂₀ and C̄₂₀ negative
  // because the planets are oblate, not prolate.  Note that this follows IERS
  // Technical Note 36 and it differs from
  // https://en.wikipedia.org/wiki/Geopotential_model which seems to want J̃₂ to
  // be negative.
  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  Degree2ZonalAcceleration(UnitVector const& axis,
                           Displacement<Frame> const& r,
                           Exponentiation<Length, -2> const& one_over_r²,
                           Exponentiation<Length, -3> const& one_over_r³) const;

  static Rotation<SurfaceFrame, Frame> FromSurfaceFrame(
      OblateBody<Frame> const& body,
      Numerics::Mode const mode,
      Instant const& t);

  not_null<OblateBody<Frame> const*> body_;
  serialization::Numerics::Mode mode_;
};

}  // namespace internal_geopotential

using internal_geopotential::Geopotential;

}  // namespace physics
}  // namespace principia

#include "physics/geopotential_body.hpp"
