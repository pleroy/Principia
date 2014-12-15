#pragma once

#include <map>
#include <utility>

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {

namespace {

// The pair<> is in diagrammatic order: right is applied first and is the first
// element of the pair, left is applied second and is the second  element.
template<typename Right, typename Left, typename Result>
struct MultiplicationMap {
  static std::map<
      std::pair<typename Right::CoordinatePermutation,
                typename Left::CoordinatePermutation>,
      typename Result::CoordinatePermutation> const map;
};

template<typename Right, typename Left, typename Result>
std::map<std::pair<typename Right::CoordinatePermutation,
                   typename Left::CoordinatePermutation>,
         typename Result::CoordinatePermutation> const
    MultiplicationMap<Right, Left, Result>::map = {
    {{Right::XYZ, Left::XYZ}, Result::XYZ},
    {{Right::XYZ, Left::YZX}, Result::YZX},
    {{Right::XYZ, Left::ZXY}, Result::ZXY},
    {{Right::XYZ, Left::XZY}, Result::XZY},
    {{Right::XYZ, Left::ZYX}, Result::ZYX},
    {{Right::XYZ, Left::YXZ}, Result::YXZ},

    {{Right::YZX, Left::XYZ}, Result::YZX},
    {{Right::YZX, Left::YZX}, Result::ZXY},
    {{Right::YZX, Left::ZXY}, Result::XYZ},
    {{Right::YZX, Left::XZY}, Result::YXZ},
    {{Right::YZX, Left::ZYX}, Result::XZY},
    {{Right::YZX, Left::YXZ}, Result::ZYX},

    {{Right::ZXY, Left::XYZ}, Result::ZXY},
    {{Right::ZXY, Left::YZX}, Result::XYZ},
    {{Right::ZXY, Left::ZXY}, Result::YZX},
    {{Right::ZXY, Left::XZY}, Result::ZYX},
    {{Right::ZXY, Left::ZYX}, Result::YXZ},
    {{Right::ZXY, Left::YXZ}, Result::XZY},

    {{Right::XZY, Left::XYZ}, Result::XZY},
    {{Right::XZY, Left::YZX}, Result::ZYX},
    {{Right::XZY, Left::ZXY}, Result::YXZ},
    {{Right::XZY, Left::XZY}, Result::XYZ},
    {{Right::XZY, Left::ZYX}, Result::YZX},
    {{Right::XZY, Left::YXZ}, Result::ZXY},

    {{Right::ZYX, Left::XYZ}, Result::ZYX},
    {{Right::ZYX, Left::YZX}, Result::YXZ},
    {{Right::ZYX, Left::ZXY}, Result::XZY},
    {{Right::ZYX, Left::XZY}, Result::ZXY},
    {{Right::ZYX, Left::ZYX}, Result::XYZ},
    {{Right::ZYX, Left::YXZ}, Result::YZX},

    {{Right::YXZ, Left::XYZ}, Result::YXZ},
    {{Right::YXZ, Left::YZX}, Result::XZY},
    {{Right::YXZ, Left::ZXY}, Result::ZYX},
    {{Right::YXZ, Left::XZY}, Result::YZX},
    {{Right::YXZ, Left::ZYX}, Result::ZXY},
    {{Right::YXZ, Left::YXZ}, Result::XYZ}};

}  // namespace

template<typename FromFrame, typename ToFrame>
Permutation<FromFrame, ToFrame>::Permutation(
    CoordinatePermutation const coordinate_permutation)
    : coordinate_permutation_(coordinate_permutation) {}

template<typename FromFrame, typename ToFrame>
inline Sign Permutation<FromFrame, ToFrame>::Determinant() const {
  return Sign(coordinate_permutation_);
}

template<typename FromFrame, typename ToFrame>
Permutation<ToFrame, FromFrame>
Permutation<FromFrame, ToFrame>::Inverse() const {
  return PTF(inverse_.at(coordinate_permutation_));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> Permutation<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return Vector<Scalar, ToFrame>((*this)(vector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> Permutation<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>(
      Determinant() * (*this)(bivector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> Permutation<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return Trivector<Scalar, ToFrame>(Determinant() * trivector.coordinates());
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame>
Permutation<FromFrame, ToFrame>::Forget() const {
  return OrthogonalMap<FromFrame, ToFrame>(
      Determinant(),
      Rotation<FromFrame, ToFrame>(quaternion_.at(coordinate_permutation_)));
}

template<typename FromFrame, typename ToFrame>
Permutation<FromFrame, ToFrame> Permutation<FromFrame, ToFrame>::Identity() {
  return Permutation(XYZ);
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
R3Element<Scalar> Permutation<FromFrame, ToFrame>::operator()(
    R3Element<Scalar> const& r3_element) const {
  R3Element<Scalar> result;
  for (int coordinate = x; coordinate <= z; ++coordinate) {
    result[coordinate] = r3_element[
        0x3 &
        (static_cast<int>(coordinate_permutation_) >> (coordinate * 2))];
  }
  return result;
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Permutation<FromFrame, ToFrame> operator*(
    Permutation<ThroughFrame, ToFrame> const& left,
    Permutation<FromFrame, ThroughFrame> const& right) {
  using Left = Permutation<ThroughFrame, ToFrame>;
  using Right = Permutation<FromFrame, ThroughFrame>;
  using Result = Permutation<FromFrame, ToFrame>;
  return Result(MultiplicationMap<Right, Left, Result>::map.at(
                    {right.coordinate_permutation_,
                     left.coordinate_permutation_}));
}

template<typename FromFrame, typename ToFrame>
std::map<typename Permutation<FromFrame,
                              ToFrame>::PFT::CoordinatePermutation,
         typename Permutation<FromFrame,
                              ToFrame>::PTF::CoordinatePermutation> const
    Permutation<FromFrame, ToFrame>::inverse_ = {
  {PFT::XYZ, PTF::XYZ},
  {PFT::YZX, PTF::ZXY},
  {PFT::ZXY, PTF::YZX},
  {PFT::XZY, PTF::XZY},
  {PFT::ZYX, PTF::ZYX},
  {PFT::YXZ, PTF::YXZ}};

template<typename FromFrame, typename ToFrame>
std::map<typename Permutation<FromFrame, ToFrame>::CoordinatePermutation,
         Quaternion> const Permutation<FromFrame, ToFrame>::quaternion_ = {
    {XYZ, Quaternion(1, {0, 0, 0})},
    {YZX, Quaternion(0.5, {-0.5, -0.5, -0.5})},
    {ZXY, Quaternion(0.5, {0.5, 0.5, 0.5})},
    {XZY, Quaternion(0, {0, -quantities::Sqrt(0.5), quantities::Sqrt(0.5)})},
    {ZYX, Quaternion(0, {-quantities::Sqrt(0.5), 0, quantities::Sqrt(0.5)})},
    {YXZ, Quaternion(0, {-quantities::Sqrt(0.5), quantities::Sqrt(0.5), 0})}};

}  // namespace geometry
}  // namespace principia
