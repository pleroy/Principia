﻿
#pragma once

#include "geometry/rotation.hpp"

#include <algorithm>

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace internal_rotation {

using base::not_null;
using quantities::Cos;
using quantities::Sin;

// Well-conditioned conversion of a rotation matrix to a quaternion.  See
// http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion and
// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/.
FORCE_INLINE Quaternion ToQuaternion(R3x3Matrix const& matrix) {
  // TODO(egg): this should probably contain some checks that |matrix| has
  // positive determinant...
  double const t = matrix.Trace();
  double real_part;
  R3Element<double> imaginary_part;
  if (t > 0) {
    double const r = sqrt(1.0 + t);
    double const s = 0.5 / r;
    real_part = 0.5 * r;
    imaginary_part.x = (matrix(2, 1) - matrix(1, 2)) * s;
    imaginary_part.y = (matrix(0, 2) - matrix(2, 0)) * s;
    imaginary_part.z = (matrix(1, 0) - matrix(0, 1)) * s;
  } else if (matrix(0, 0) > std::max(matrix(1, 1), matrix(2, 2))) {
    double const r =
        sqrt(1.0 + matrix(0, 0) - matrix(1, 1) - matrix(2, 2));
    double const s = 0.5 / r;
    real_part = (matrix(2, 1) - matrix(1, 2)) * s;
    imaginary_part.x = 0.5 * r;
    imaginary_part.y = (matrix(0, 1) + matrix(1, 0)) * s;
    imaginary_part.z = (matrix(2, 0) + matrix(0, 2)) * s;
  } else if (matrix(1, 1) > matrix(2, 2)) {
    double const r =
        sqrt(1.0 - matrix(0, 0) + matrix(1, 1) - matrix(2, 2));
    double const s = 0.5 / r;
    real_part = (matrix(0, 2) - matrix(2, 0)) * s;
    imaginary_part.x = (matrix(0, 1) + matrix(1, 0)) * s;
    imaginary_part.y = 0.5 * r;
    imaginary_part.z = (matrix(1, 2) + matrix(2, 1)) * s;
  } else {
    double const r =
        sqrt(1.0 - matrix(0, 0) - matrix(1, 1) + matrix(2, 2));
    double const s = 0.5 / r;
    real_part = (matrix(1, 0) - matrix(0, 1)) * s;
    imaginary_part.x = (matrix(0, 2) + matrix(2, 0)) * s;
    imaginary_part.y = (matrix(1, 2) + matrix(2, 1)) * s;
    imaginary_part.z = 0.5 * r;
  }
  return Quaternion(real_part, imaginary_part);
}

// Returns a rotation of |angle| around |axis|.  |axis| must be normalized.
inline Quaternion AngleAxis(Angle const& angle, R3Element<double> const& axis) {
  quantities::Angle const half_angle = 0.5 * angle;
  return Quaternion(Cos(half_angle), Sin(half_angle) * axis);
}

// Returns the digits of the 3ⁿs from the given |BinaryCodedTernary number|.
// Note that this does not check that |number| is valid binary-coded ternary,
// nor that the result is between 1 and 2.
template<typename BinaryCodedTernary>
int BinaryCodedTernaryDigit(int const n, BinaryCodedTernary const number) {
  return (static_cast<int>(number) >> (2 * n)) & 0b11;
}

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame>::Rotation(Quaternion const& quaternion)
    : quaternion_(quaternion) {}

template<typename FromFrame, typename ToFrame>
template<typename Scalar, typename F, typename T, typename>
Rotation<FromFrame, ToFrame>::Rotation(quantities::Angle const& angle,
                                       Bivector<Scalar, FromFrame> const& axis)
    : Rotation(AngleAxis(angle, Normalize(axis).coordinates())) {}

template<typename FromFrame, typename ToFrame>
template<int rank_x, int rank_y, int rank_z, typename F, typename T, typename>
Rotation<FromFrame, ToFrame>::Rotation(
    Multivector<double, FromFrame, rank_x> x_to_frame,
    Multivector<double, FromFrame, rank_y> y_to_frame,
    Multivector<double, FromFrame, rank_z> z_to_frame)
    : Rotation<FromFrame, ToFrame>(
          ToQuaternion(R3x3Matrix(x_to_frame.coordinates(),
                                  y_to_frame.coordinates(),
                                  z_to_frame.coordinates()))) {
  static_assert((rank_x + rank_y + rank_z) % 2 == 0, "chiral basis");
  static_assert(rank_x < 3 && rank_y < 3 && rank_z < 3, "bad dimension");
}

template<typename FromFrame, typename ToFrame>
template<int rank_x, int rank_y, int rank_z,
         typename F, typename T, typename, typename>  // typename and spam.
Rotation<FromFrame, ToFrame>::Rotation(
    Multivector<double, ToFrame, rank_x> x_from_frame,
    Multivector<double, ToFrame, rank_y> y_from_frame,
    Multivector<double, ToFrame, rank_z> z_from_frame)
    : Rotation<FromFrame, ToFrame>(
          ToQuaternion(R3x3Matrix(x_from_frame.coordinates(),
                                  y_from_frame.coordinates(),
                                  z_from_frame.coordinates()).Transpose())) {
  static_assert((rank_x + rank_y + rank_z) % 2 == 0, "chiral basis");
  static_assert(rank_x < 3 && rank_y < 3 && rank_z < 3, "bad dimension");
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar, typename F, typename T, typename>
Rotation<FromFrame, ToFrame>::Rotation(Angle const& angle,
                                       Bivector<Scalar, FromFrame> const& axis,
                                       DefinesFrame<ToFrame> tag)
    : Rotation(AngleAxis(-angle, Normalize(axis).coordinates())) {}

template<typename FromFrame, typename ToFrame>
template<typename Scalar, typename F, typename T, typename, typename>
Rotation<FromFrame, ToFrame>::Rotation(Angle const& angle,
                                       Bivector<Scalar, ToFrame> const& axis,
                                       DefinesFrame<FromFrame> tag)
    : Rotation(AngleAxis(angle, Normalize(axis).coordinates())) {}

template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename>
Rotation<FromFrame, ToFrame>::Rotation(
    Angle const& α,
    Angle const& β,
    Angle const& γ,
    EulerAngles const axes,
    DefinesFrame<ToFrame> tag)
    : Rotation(Rotation<ToFrame, FromFrame>(α, β, γ, axes, tag).Inverse()) {}

template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename, typename>
Rotation<FromFrame, ToFrame>::Rotation(
    Angle const& α,
    Angle const& β,
    Angle const& γ,
    EulerAngles const axes,
    DefinesFrame<FromFrame> tag)
    : Rotation(AngleAxis(α, BasisVector(BinaryCodedTernaryDigit(2, axes))) *
               AngleAxis(β, BasisVector(BinaryCodedTernaryDigit(1, axes))) *
               AngleAxis(γ, BasisVector(BinaryCodedTernaryDigit(0, axes)))) {}

template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename>
Rotation<FromFrame, ToFrame>::Rotation(
    Angle const& α,
    Angle const& β,
    Angle const& γ,
    CardanoAngles const axes,
    DefinesFrame<ToFrame> tag)
    : Rotation(Rotation<ToFrame, FromFrame>(α, β, γ, axes, tag).Inverse()) {}

template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename, typename>
Rotation<FromFrame, ToFrame>::Rotation(
    Angle const& α,
    Angle const& β,
    Angle const& γ,
    CardanoAngles const axes,
    DefinesFrame<FromFrame> tag)
    : Rotation(AngleAxis(α, BasisVector(BinaryCodedTernaryDigit(2, axes))) *
               AngleAxis(β, BasisVector(BinaryCodedTernaryDigit(1, axes))) *
               AngleAxis(γ, BasisVector(BinaryCodedTernaryDigit(0, axes)))) {}

template<typename FromFrame, typename ToFrame>
Sign Rotation<FromFrame, ToFrame>::Determinant() const {
  return Sign(1);
}

template<typename FromFrame, typename ToFrame>
Rotation<ToFrame, FromFrame> Rotation<FromFrame, ToFrame>::Inverse() const {
  // Because |quaternion_| has norm 1, its inverse is just its conjugate.
  return Rotation<ToFrame, FromFrame>(quaternion_.Conjugate());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> Rotation<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return Vector<Scalar, ToFrame>((*this)(vector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> Rotation<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>((*this)(bivector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> Rotation<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return trivector;
}

template<typename FromFrame, typename ToFrame>
template<typename T>
typename base::Mappable<Rotation<FromFrame, ToFrame>, T>::type
Rotation<FromFrame, ToFrame>::operator()(T const& t) const {
  return base::Mappable<Rotation, T>::Do(*this, t);
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> Rotation<FromFrame, ToFrame>::Forget() const {
  return OrthogonalMap<FromFrame, ToFrame>(Sign(1), *this);
}

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> Rotation<FromFrame, ToFrame>::Identity() {
  return Rotation(Quaternion(1));
}

template<typename FromFrame, typename ToFrame>
Quaternion const& Rotation<FromFrame, ToFrame>::quaternion() const {
  return quaternion_;
}

template<typename FromFrame, typename ToFrame>
void Rotation<FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::LinearMap*> const message) const {
  LinearMap<FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(message->MutableExtension(serialization::Rotation::extension));
}

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> Rotation<FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::Rotation::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::Rotation::extension));
}

template<typename FromFrame, typename ToFrame>
void Rotation<FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::Rotation*> const message) const {
  quaternion_.WriteToMessage(message->mutable_quaternion());
}

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> Rotation<FromFrame, ToFrame>::ReadFromMessage(
    serialization::Rotation const& message) {
  return Rotation(Quaternion::ReadFromMessage(message.quaternion()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
R3Element<Scalar> Rotation<FromFrame, ToFrame>::operator()(
    R3Element<Scalar> const& r3_element) const {
  double const real_part = quaternion_.real_part();
  R3Element<double> const& imaginary_part = quaternion_.imaginary_part();
  return r3_element + 2 * Cross(imaginary_part,
                                Cross(imaginary_part, r3_element) +
                                      real_part * r3_element);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> operator*(
    Rotation<ThroughFrame, ToFrame> const& left,
    Rotation<FromFrame, ThroughFrame> const& right) {
  return Rotation<FromFrame, ToFrame>(left.quaternion_ * right.quaternion_);
}

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Rotation<FromFrame, ToFrame> const& rotation) {
  return out << rotation.quaternion_;
}

}  // namespace internal_rotation
}  // namespace geometry
}  // namespace principia
