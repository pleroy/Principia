﻿
#pragma once

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace internal_elementary_functions {

// Equivalent to |std::fma(x, y, z)|.
template<typename Q1, typename Q2,
         typename = std::enable_if<is_quantity<Q1>::value>,
         typename = std::enable_if<is_quantity<Q2>::value>>
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z);

// Equivalent to |std::abs(x)|.
template<typename Q, typename = std::enable_if<is_quantity<Q>::value>>
Q Abs(Q const& x);

// Equivalent to |std::sqrt(x)|.
template<typename Q, typename = std::enable_if<is_quantity<Q>::value>>
SquareRoot<Q> Sqrt(Q const& x);

// Equivalent to |std::cbrt(x)|.
template<typename Q, typename = std::enable_if<is_quantity<Q>::value>>
CubeRoot<Q> Cbrt(Q const& x);

// Equivalent to |std::pow(x, exponent)| unless -3 ≤ x ≤ 3, in which case
// explicit specialization yields multiplications statically.
template<int exponent, typename Q,
         typename = std::enable_if<is_quantity<Q>::value>>
constexpr Exponentiation<Q, exponent> Pow(Q const& x);

double Sin(Angle const& α);
double Cos(Angle const& α);
double Tan(Angle const& α);

Angle ArcSin(double x);
Angle ArcCos(double x);
Angle ArcTan(double y, double x = 1);
template<typename D>
Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x);

// We consider hyperbolic functions as dealing with quotients of arc length to
// curvature radius in the hyperbolic plane, which are angles. This explains the
// use of "arc" for inverse functions.

double Sinh(Angle const& α);
double Cosh(Angle const& α);
double Tanh(Angle const& α);

Angle ArcSinh(double x);
Angle ArcCosh(double x);
Angle ArcTanh(double x);

}  // namespace internal_elementary_functions

using internal_elementary_functions::Abs;
using internal_elementary_functions::ArcCos;
using internal_elementary_functions::ArcCosh;
using internal_elementary_functions::ArcSin;
using internal_elementary_functions::ArcSinh;
using internal_elementary_functions::ArcTan;
using internal_elementary_functions::ArcTanh;
using internal_elementary_functions::Cbrt;
using internal_elementary_functions::Cos;
using internal_elementary_functions::Cosh;
using internal_elementary_functions::FusedMultiplyAdd;
using internal_elementary_functions::Pow;
using internal_elementary_functions::Sin;
using internal_elementary_functions::Sinh;
using internal_elementary_functions::Sqrt;
using internal_elementary_functions::Tan;
using internal_elementary_functions::Tanh;

}  // namespace quantities
}  // namespace principia

#include "quantities/elementary_functions_body.hpp"
