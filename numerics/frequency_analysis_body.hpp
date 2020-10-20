
#pragma once

#include "numerics/frequency_analysis.hpp"

#include <algorithm>
#include <functional>
#include <vector>

#include "absl/strings/str_cat.h"
#include "base/tags.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/poisson_series_basis.hpp"
#include "numerics/root_finders.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using base::uninitialized;
using geometry::Hilbert;
using geometry::Vector;
using quantities::Inverse;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

template<typename Function,
         int wdegree_,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight) {
  auto const weighted_function = weight * function;
  auto const weighted_function_spectrum = weighted_function.FourierTransform();

  auto power =
      [&weighted_function_spectrum](AngularFrequency const& ω) {
        return weighted_function_spectrum(ω).Norm²();
      };

  return Brent(power,
               fft_mode.min,
               fft_mode.max,
               std::greater<>());
}

template<int degree_,
         typename Function,
         int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
Projection(Function const& function,
           AngularFrequency const& ω,
           PoissonSeries<double, wdegree_, Evaluator> const& weight,
           Instant const& t_min,
           Instant const& t_max) {
  std::optional<AngularFrequency> optional_ω = ω;

  // A calculator that returns optional_ω once and then stops.
  auto angular_frequency_calculator = [&optional_ω](auto const& residual) {
    auto const result = optional_ω;
    optional_ω = std::nullopt;
    return result;
  };

  return IncrementalProjection<degree_>(function,
                                        angular_frequency_calculator,
                                        weight,
                                        t_min, t_max);
}

int iter = 0;
mathematica::Logger logger(TEMP_DIR / "frequency_analysis.wl",
                           /*make_unique=*/false);

template<int degree_,
         typename Function,
         typename AngularFrequencyCalculator, int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
IncrementalProjection(Function const& function,
                      AngularFrequencyCalculator const& calculator,
                      PoissonSeries<double, wdegree_, Evaluator> const& weight,
                      Instant const& t_min,
                      Instant const& t_max) {
  logger.Set("tMin", t_min, mathematica::ExpressIn(Metre, Second, Radian));
  logger.Set("tMax", t_max, mathematica::ExpressIn(Metre, Second, Radian));

  using Value = std::invoke_result_t<Function, Instant>;
  using Norm = typename Hilbert<Value>::NormType;
  using Normalized = typename Hilbert<Value>::NormalizedType;
  using Series = PoissonSeries<Value, degree_, Evaluator>;

  // This code follows [Kud07], section 2.  Our indices start at 0, unlike those
  // of Кудрявцев which start at 1.

  Instant const& t0 = weight.origin();

  std::optional<AngularFrequency> ω = calculator(function);
  CHECK(ω.has_value());

  std::vector<Series> basis;

  int basis_size;
  if (ω.value() == AngularFrequency{}) {
    auto const ω_basis =
        PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
            t0);
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
  } else {
    auto const ω_basis =
        PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
            ω.value(), t0);
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
  }
  logger.Append(absl::StrCat("frequency[", iter, "]"),
                ω.value(),
                mathematica::ExpressIn(Metre, Second, Radian));
  logger.Append(absl::StrCat("trajectory[", iter, "]"),
                function,
                mathematica::ExpressIn(Metre, Second, Radian));

  // This is logically Q in the QR decomposition of basis.
  std::vector<PoissonSeries<Normalized, degree_, Evaluator>> q;

#define USE_INTEGRATE 0
  auto const a₀ = basis[0];
#if USE_INTEGRATE
  auto const r₀₀ = Sqrt((PointwiseInnerProduct(a₀, a₀) * weight)
                         .Integrate(t_min, t_max) /
                     (t_max - t_min));
#else
  auto const r₀₀ = a₀.Norm(weight, t_min, t_max);
#endif
  q.push_back(a₀ / r₀₀);

  logger.Append(absl::StrCat("basis[", iter, "]"),
                basis[0],
                mathematica::ExpressIn(Metre, Second, Radian));
  logger.Append(absl::StrCat("normFn[", iter, "]"),
                a₀,
                mathematica::ExpressIn(Metre, Second, Radian));
  logger.Append(absl::StrCat("normIntegrand[", iter, "]"),
                PointwiseInnerProduct(a₀, a₀) * weight,
                mathematica::ExpressIn(Metre, Second, Radian));
  logger.Append(absl::StrCat("norm[", iter, "]"),
                a₀.Norm(weight, t_min, t_max),
                mathematica::ExpressIn(Metre, Second, Radian));
  logger.Append(absl::StrCat("normIntegrate[", iter, "]"),
                Sqrt((PointwiseInnerProduct(a₀, a₀) * weight)
                         .Integrate(t_min, t_max) /
                     (t_max - t_min)),
                mathematica::ExpressIn(Metre, Second, Radian));
  logger.Append(absl::StrCat("q[", iter, "]"),
                q[0],
                mathematica::ExpressIn(Metre, Second, Radian));

  auto const A₀ = InnerProduct(function, q[0], weight, t_min, t_max);

  PoissonSeries<Value, degree_, Evaluator> F = A₀ * q[0];
  auto f = function - F;

  int m_begin = 1;
  for (;;) {
    for (int m = m_begin; m < basis_size; ++m) {
      auto aₘ⁽ᵏ⁾ = basis[m];
      for (int k = 0; k < m; ++k) {
#if USE_INTEGRATE
        auto const rₖₘ = (PointwiseInnerProduct(q[k], aₘ⁽ᵏ⁾) * weight)
                          .Integrate(t_min, t_max) / (t_max - t_min);
#else
        auto const rₖₘ = InnerProduct(q[k], aₘ⁽ᵏ⁾, weight, t_min, t_max);
#endif
        logger.Append(absl::StrCat("innerProduct[", iter, "]"),
                      InnerProduct(q[k], aₘ⁽ᵏ⁾, weight, t_min, t_max),
                      mathematica::ExpressIn(Metre, Second, Radian));
        logger.Append(absl::StrCat("innerProductIntegrate[", iter, "]"),
                      (PointwiseInnerProduct(q[k], aₘ⁽ᵏ⁾) * weight)
                          .Integrate(t_min, t_max) / (t_max - t_min),
                      mathematica::ExpressIn(Metre, Second, Radian));
        aₘ⁽ᵏ⁾ -= rₖₘ * q[k];
      }

#if USE_INTEGRATE
      auto const rₘₘ = Sqrt((PointwiseInnerProduct(aₘ⁽ᵏ⁾, aₘ⁽ᵏ⁾) * weight)
                             .Integrate(t_min, t_max) /
                         (t_max - t_min));
#else
      auto const rₘₘ = aₘ⁽ᵏ⁾.Norm(weight, t_min, t_max);
#endif
      logger.Append(absl::StrCat("normFn[", iter, "]"),
                    aₘ⁽ᵏ⁾,
                    mathematica::ExpressIn(Metre, Second, Radian));
      logger.Append(absl::StrCat("normIntegrand[", iter, "]"),
                    PointwiseInnerProduct(aₘ⁽ᵏ⁾, aₘ⁽ᵏ⁾) * weight,
                    mathematica::ExpressIn(Metre, Second, Radian));
      logger.Append(absl::StrCat("norm[", iter, "]"),
                    aₘ⁽ᵏ⁾.Norm(weight, t_min, t_max),
                    mathematica::ExpressIn(Metre, Second, Radian));
      logger.Append(absl::StrCat("normIntegrate[", iter, "]"),
                    Sqrt((PointwiseInnerProduct(aₘ⁽ᵏ⁾, aₘ⁽ᵏ⁾) * weight)
                             .Integrate(t_min, t_max) /
                         (t_max - t_min)),
                    mathematica::ExpressIn(Metre, Second, Radian));
      q.push_back(aₘ⁽ᵏ⁾ / rₘₘ);
      logger.Append(absl::StrCat("q[", iter, "]"),
                    q.back(),
                    mathematica::ExpressIn(Metre, Second, Radian));
      DCHECK_EQ(m + 1, q.size());

      Norm const Aₘ = InnerProduct(f, q[m], weight, t_min, t_max);

      f -= Aₘ * q[m];
      F += Aₘ * q[m];
    }

    logger.Append(absl::StrCat("solution[", iter, "]"),
                  F,
                  mathematica::ExpressIn(Metre, Second, Radian));
    ω = calculator(f);
    if (!ω.has_value()) {
      ++iter;
      return F;
    }
    logger.Append(absl::StrCat("frequency[", iter, "]"),
                  ω.value(),
                  mathematica::ExpressIn(Metre, Second, Radian));

    int ω_basis_size;
    if (ω.value() == AngularFrequency{}) {
      auto const ω_basis =
          PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
              t0);
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    } else {
      auto const ω_basis =
          PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
              ω.value(), t0);
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    }
    m_begin = basis_size;
    basis_size += ω_basis_size;
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
