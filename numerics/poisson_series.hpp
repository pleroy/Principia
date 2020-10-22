﻿
#pragma once

#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/complexification.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/polynomial.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
FORWARD_DECLARE_FROM(poisson_series,
                     TEMPLATE(typename Value, int degree_,
                              template<typename, typename, int> class Evaluator)
                              class,
                     PoissonSeries);
FORWARD_DECLARE_FROM(poisson_series,
                     TEMPLATE(typename Value, int degree_,
                              template<typename, typename, int> class Evaluator)
                              class,
                     PiecewisePoissonSeries);
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION_FROM(
    mathematica,
    TEMPLATE(typename Value, int degree_,
             template<typename, typename, int> class Evaluator,
             typename OptionalExpressIn) std::string,
    ToMathematicaExpression,
    (numerics::PoissonSeries<Value, degree_, Evaluator> const& series,
     OptionalExpressIn express_in));
FORWARD_DECLARE_FUNCTION_FROM(
    mathematica,
    TEMPLATE(typename Value, int degree_,
             template<typename, typename, int> class Evaluator,
             typename OptionalExpressIn) std::string,
    ToMathematicaExpression,
    (numerics::PiecewisePoissonSeries<Value, degree_, Evaluator> const& series,
     OptionalExpressIn express_in));
}  // namespace mathematica

namespace numerics {
namespace internal_poisson_series {

using base::not_null;
using geometry::Complexification;
using geometry::Hilbert;
using geometry::Instant;
using geometry::Interval;
using quantities::AngularFrequency;
using quantities::Primitive;
using quantities::Product;
using quantities::Quotient;
using quantities::Time;

// The trigonometric functions are by default assumed to look like a polynomial
// of this degree over an interval of a piecewise series.
constexpr int estimated_trigonometric_degree = 1;

// A Poisson series is the sum of terms of the form:
//   aₙtⁿ      aₙₖ tⁿ sin ωₖ t      aₙₖ tⁿ cos ωₖ t
// Terms of the first kind are called aperiodic, terms of the second and third
// kind are called periodic.  Poisson series form an algebra that is stable by
// derivation and integration.
template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
class PoissonSeries {
 public:
  static constexpr int degree = degree_;
  using Polynomial =
      numerics::PolynomialInMonomialBasis<Value, Instant, degree_, Evaluator>;

  // TODO(phl): Use designated initializers for this struct once this project
  // can be compiled using c++latest.
  // TODO(phl): If we wanted to have Poisson series returning affine values,
  // these polynomials should be changed to return Difference<Value>.
  struct Polynomials {
    Polynomial sin;
    Polynomial cos;
  };

  using PolynomialsByAngularFrequency =
      std::vector<std::pair<AngularFrequency, Polynomials>>;

  // The |periodic| vector may contain frequencies that are negative or zero, as
  // well as repeated frequencies.  It is not expected to be ordered.
  PoissonSeries(Polynomial const& aperiodic,
                PolynomialsByAngularFrequency const& periodic);

  // A Poisson series may be explicitly converted to a higher degree (possibly
  // with a different evaluator).
  template<int higher_degree_,
           template<typename, typename, int> class HigherEvaluator>
  explicit operator
      PoissonSeries<Value, higher_degree_, HigherEvaluator>() const;

  Instant const& origin() const;

  Value operator()(Instant const& t) const;

  // Returns a copy of this series adjusted to the given origin.
  PoissonSeries AtOrigin(Instant const& origin) const;

  // The constant term of the result is zero.
  PoissonSeries<quantities::Primitive<Value, Time>, degree_ + 1, Evaluator>
  Primitive() const;

  quantities::Primitive<Value, Time> Integrate(Instant const& t1,
                                               Instant const& t2) const;

  template<int wdegree_>
  typename Hilbert<Value>::NormType Norm(
      PoissonSeries<double, wdegree_, Evaluator> const& weight,
      Instant const& t_min,
      Instant const& t_max) const;

  template<int d>
  PoissonSeries& operator+=(PoissonSeries<Value, d, Evaluator> const& right);
  template<int d>
  PoissonSeries& operator-=(PoissonSeries<Value, d, Evaluator> const& right);

  void WriteToMessage(not_null<serialization::PoissonSeries*> message) const;
  static PoissonSeries ReadFromMessage(
      serialization::PoissonSeries const& message);

 private:
  // Similar to the public constructor, but passing by copy allows moves, which
  // is useful for internal algorithms.
  struct PrivateConstructor {};
  PoissonSeries(PrivateConstructor,
                Polynomial aperiodic,
                PolynomialsByAngularFrequency periodic);

  // Similar to the previous constructor, except that the |periodic| vector is
  // used verbatim, without sorting or normalization, which is useful for
  // internal algorithms which produce positive, ordered frequencies.
  struct TrustedPrivateConstructor {};
  PoissonSeries(TrustedPrivateConstructor,
                Polynomial aperiodic,
                PolynomialsByAngularFrequency periodic);

  //TODO(phl): Comment.
  struct SplitPoissonSeries {
    PoissonSeries slow;
    PoissonSeries fast;
  };
  SplitPoissonSeries Split(AngularFrequency const& ω_cutoff) const;

  Instant origin_;  // Common to all polynomials.
  Polynomial aperiodic_;
  // The frequencies in this vector are positive, distinct and in increasing
  // order.
  PolynomialsByAngularFrequency periodic_;

  template<typename V, int r, template<typename, typename, int> class E>
  PoissonSeries<V, r, E>
  friend operator-(PoissonSeries<V, r, E> const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PoissonSeries<V, std::max(l, r), E>
  friend operator+(PoissonSeries<V, l, E> const& left,
                   PoissonSeries<V, r, E> const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PoissonSeries<V, std::max(l, r), E>
  friend operator-(PoissonSeries<V, l, E> const& left,
                   PoissonSeries<V, r, E> const& right);
  template<typename Scalar,
           typename V, int d,
           template<typename, typename, int> class E>
  PoissonSeries<Product<Scalar, V>, d, E>
  friend operator*(Scalar const& left,
                   PoissonSeries<V, d, E> const& right);
  template<typename Scalar,
           typename V, int d,
           template<typename, typename, int> class E>
  PoissonSeries<Product<V, Scalar>, d, E>
  friend operator*(PoissonSeries<V, d, E> const& left,
                   Scalar const& right);
  template<typename Scalar,
           typename V, int d,
           template<typename, typename, int> class E>
  PoissonSeries<Quotient<V, Scalar>, d, E>
  friend operator/(PoissonSeries<V, d, E> const& left,
                   Scalar const& right);
  template<typename L, typename R,
           int l, int r,
           template<typename, typename, int> class E,
           typename P>
  auto friend Multiply(PoissonSeries<L, l, E> const& left,
                       PoissonSeries<R, r, E> const& right,
                       P const& product);
  template<typename V, int d, template<typename, typename, int> class E>
  friend std::ostream& operator<<(std::ostream& out,
                                  PoissonSeries<V, d, E> const& series);
  template<typename L, typename R,
         int l, int r, int w,
         template<typename, typename, int> class E>
  friend typename Hilbert<L, R>::InnerProductType InnerProduct(
      PoissonSeries<L, l, E> const& left,
      PoissonSeries<R, r, E> const& right,
      PoissonSeries<double, w, E> const& weight,
      Instant const& t_min,
      Instant const& t_max);
  template<typename V, int d,
           template<typename, typename, int> class E,
           typename O>
  friend std::string mathematica::internal_mathematica::ToMathematicaExpression(
      PoissonSeries<V, d, E> const& polynomial,
      O express_in);
  template<typename V, int d,
           template<typename, typename, int> class E>
  friend class PoissonSeries;
};

// Vector space of Poisson series.

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, rdegree_, Evaluator>
operator+(PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, rdegree_, Evaluator>
operator-(PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Scalar, Value>, degree_, Evaluator>
operator*(Scalar const& left,
          PoissonSeries<Value, degree_, Evaluator> const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<Value, Scalar>, degree_, Evaluator>
operator*(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>
operator/(PoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right);

// Algebra of Poisson series.

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>
operator*(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
          PoissonSeries<RValue, rdegree_, Evaluator> const& right);

// Returns a scalar-valued Poisson series obtained by pointwise inner product of
// two vector-valued series.
template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<typename Hilbert<LValue, RValue>::InnerProductType,
              ldegree_ + rdegree_,
              Evaluator>
PointwiseInnerProduct(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
                      PoissonSeries<RValue, rdegree_, Evaluator> const& right);

// Output.

template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
std::ostream& operator<<(
    std::ostream& out,
    PoissonSeries<Value, degree_, Evaluator> const& series);

// Inner product space of Poisson series.

// Technically the weight function must be nonnegative for this to be an inner
// product.  Not sure how this works with the flat-top windows, which can be
// negative.  Note that the result is normalized by dividing by (t_max - t_min).
template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
             PoissonSeries<RValue, rdegree_, Evaluator> const& right,
             PoissonSeries<double, wdegree_, Evaluator> const& weight,
             Instant const& t_min,
             Instant const& t_max);

// A function defined by Poisson series piecewise.  Each of the Poisson series
// making up the function applies over the semi-open interval
// [internal.min, interval.max[.  It's not required that the function be
// continuous.
template<typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
class PiecewisePoissonSeries {
 public:
  static constexpr int degree = degree_;
  using Series = PoissonSeries<Value, degree_, Evaluator>;
  using Spectrum = std::function<Complexification<Primitive<Value, Time>>(
      AngularFrequency const&)>;

  PiecewisePoissonSeries(Interval<Instant> const& interval,
                         Series const& series);

  // The intervals for successive calls to Append must be consecutive.  For the
  // first call, the interval must be consecutive with the one passed at
  // construction.
  void Append(Interval<Instant> const& interval,
              Series const& series);

  Instant t_min() const;
  Instant t_max() const;

  // t must be in the interval [t_min, t_max].
  Value operator()(Instant const& t) const;

  // Returns the Fourier transform of this piecewise Poisson series.
  // The function is taken to be 0 outside [t_min, t_max].
  // The convention used is ∫ f(t) exp(-iωt) dt, corresponding to Mathematica’s
  // FourierParameters -> {1, -1} for FourierTransform (the “pure mathematics;
  // systems engineering”).
  // When evaluated at a given frequency, the Fourier transform is computed by
  // Gauss-Legendre quadrature on each subinterval, where the number of points
  // is chosen assuming that the periods of periodic terms are all large
  // compared to the subintervals.
  // If apodization is desired, |*this| should be multiplied by an apodization
  // function, and |FourierTransform| should be called on the product.
  // |*this| must outlive the resulting function.
  Spectrum FourierTransform() const;

  template<int d>
  PiecewisePoissonSeries& operator+=(
      PoissonSeries<Value, d, Evaluator> const& right);
  template<int d>
  PiecewisePoissonSeries& operator-=(
      PoissonSeries<Value, d, Evaluator> const& right);

  void WriteToMessage(
      not_null<serialization::PiecewisePoissonSeries*> message) const;
  static PiecewisePoissonSeries ReadFromMessage(
      serialization::PiecewisePoissonSeries const& message);

 private:
  PiecewisePoissonSeries(std::vector<Instant> const& bounds,
                         std::vector<Series> const& series,
                         std::optional<Series> const& addend);

  Value EvaluateAddend(Instant const& t) const;

  std::vector<Instant> bounds_;
  std::vector<Series> series_;
  std::optional<Series> addend_;

  template<typename V, int r, template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, r, E>
  friend operator-(PiecewisePoissonSeries<V, r, E> const& right);
  template<typename S, typename V, int d,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<S, V>, d, E>
  friend operator*(S const& left,
                   PiecewisePoissonSeries<V, d, E> const& right);
  template<typename S, typename V, int d,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<V, S>, d, E>
  friend operator*(PiecewisePoissonSeries<V, d, E> const& left,
                   S const& right);
  template<typename S, typename V, int d,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Quotient<V, S>, d, E>
  friend operator/(PiecewisePoissonSeries<V, d, E> const& left,
                   S const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, std::max(l, r), E>
  friend operator+(PoissonSeries<V, l, E> const& left,
                   PiecewisePoissonSeries<V, r, E> const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, std::max(l, r), E>
  friend operator+(PiecewisePoissonSeries<V, l, E> const& left,
                   PoissonSeries<V, r, E> const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, std::max(l, r), E>
  friend operator-(PoissonSeries<V, l, E> const& left,
                   PiecewisePoissonSeries<V, r, E> const& right);
  template<typename V, int l, int r, template<typename, typename, int> class E>
  PiecewisePoissonSeries<V, std::max(l, r), E>
  friend operator-(PiecewisePoissonSeries<V, l, E> const& left,
                   PoissonSeries<V, r, E> const& right);
  template<typename L, typename R, int l, int r,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<L, R>, l + r, E>
  friend operator*(PoissonSeries<L, l, E> const& left,
                   PiecewisePoissonSeries<R, r, E> const& right);
  template<typename L, typename R, int l, int r,
           template<typename, typename, int> class E>
  PiecewisePoissonSeries<Product<L, R>, l + r, E>
  friend operator*(PiecewisePoissonSeries<L, l, E> const& left,
                   PoissonSeries<R, r, E> const& right);
  template<typename L, typename R, int l, int r, int w,
           template<typename, typename, int> class E, int p>
  typename Hilbert<L, R>::InnerProductType
  friend InnerProduct(PoissonSeries<L, l, E> const& left,
                      PiecewisePoissonSeries<R, r, E> const& right,
                      PoissonSeries<double, w, E> const& weight,
                      Instant const& t_min,
                      Instant const& t_max);
  template<typename L, typename R, int l, int r, int w,
           template<typename, typename, int> class E, int p>
  typename Hilbert<L, R>::InnerProductType
  friend InnerProduct(PiecewisePoissonSeries<L, l, E> const& left,
                      PoissonSeries<R, r, E> const& right,
                      PoissonSeries<double, w, E> const& weight,
                      Instant const& t_min,
                      Instant const& t_max);
  template<typename V, int d,
           template<typename, typename, int> class E,
           typename O>
  friend std::string mathematica::internal_mathematica::ToMathematicaExpression(
      PiecewisePoissonSeries<V, d, E> const& polynomial,
      O express_in);
};

// Some of the vector space operations for piecewise Poisson series.  Note that
// while piecewise Poisson series have an algebra structure, we do not need it.

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, rdegree_, Evaluator>
operator+(PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, rdegree_, Evaluator>
operator-(PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<Scalar, Value>, degree_, Evaluator>
operator*(Scalar const& left,
          PiecewisePoissonSeries<Value, degree_, Evaluator> const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<Value, Scalar>, degree_, Evaluator>
operator*(PiecewisePoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right);

template<typename Scalar,
         typename Value, int degree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Quotient<Value, Scalar>, degree_, Evaluator>
operator/(PiecewisePoissonSeries<Value, degree_, Evaluator> const& left,
          Scalar const& right);

// Action of Poisson series on piecewise Poisson series.

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator+(PiecewisePoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PoissonSeries<Value, ldegree_, Evaluator> const& left,
          PiecewisePoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename Value, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Value, std::max(ldegree_, rdegree_), Evaluator>
operator-(PiecewisePoissonSeries<Value, ldegree_, Evaluator> const& left,
          PoissonSeries<Value, rdegree_, Evaluator> const& right);

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>
operator*(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
          PiecewisePoissonSeries<RValue, rdegree_, Evaluator> const& right);

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
PiecewisePoissonSeries<Product<LValue, RValue>, ldegree_ + rdegree_, Evaluator>
operator*(PiecewisePoissonSeries<LValue, ldegree_, Evaluator> const& left,
          PoissonSeries<RValue, rdegree_, Evaluator> const& right);

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator,
         int points = (ldegree_ + estimated_trigonometric_degree +
                       rdegree_ + estimated_trigonometric_degree +
                       wdegree_ + estimated_trigonometric_degree) / 2>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
             PiecewisePoissonSeries<RValue, rdegree_, Evaluator> const& right,
             PoissonSeries<double, wdegree_, Evaluator> const& weight);

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator,
         int points = (ldegree_ + estimated_trigonometric_degree +
                       rdegree_ + estimated_trigonometric_degree +
                       wdegree_ + estimated_trigonometric_degree) / 2>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PoissonSeries<LValue, ldegree_, Evaluator> const& left,
             PiecewisePoissonSeries<RValue, rdegree_, Evaluator> const& right,
             PoissonSeries<double, wdegree_, Evaluator> const& weight,
             Instant const& t_min,
             Instant const& t_max);

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator,
         int points = (ldegree_ + estimated_trigonometric_degree +
                       rdegree_ + estimated_trigonometric_degree +
                       wdegree_ + estimated_trigonometric_degree) / 2>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PiecewisePoissonSeries<LValue, ldegree_, Evaluator> const& left,
             PoissonSeries<RValue, rdegree_, Evaluator> const& right,
             PoissonSeries<double, wdegree_, Evaluator> const& weight);

template<typename LValue, typename RValue,
         int ldegree_, int rdegree_, int wdegree_,
         template<typename, typename, int> class Evaluator,
         int points = (ldegree_ + estimated_trigonometric_degree +
                       rdegree_ + estimated_trigonometric_degree +
                       wdegree_ + estimated_trigonometric_degree) / 2>
typename Hilbert<LValue, RValue>::InnerProductType
InnerProduct(PiecewisePoissonSeries<LValue, ldegree_, Evaluator> const& left,
             PoissonSeries<RValue, rdegree_, Evaluator> const& right,
             PoissonSeries<double, wdegree_, Evaluator> const& weight,
             Instant const& t_min,
             Instant const& t_max);

}  // namespace internal_poisson_series

using internal_poisson_series::InnerProduct;
using internal_poisson_series::PiecewisePoissonSeries;
using internal_poisson_series::PoissonSeries;

}  // namespace numerics
}  // namespace principia

#include "numerics/poisson_series_body.hpp"
