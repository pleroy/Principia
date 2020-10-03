﻿
#pragma once

#include "numerics/frequency_analysis.hpp"

#include <algorithm>
#include <functional>
#include <vector>

#include "base/tags.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
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
using quantities::Abs;
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

template<int degree_,
         typename Function,
         typename AngularFrequencyCalculator, int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
IncrementalProjection(Function const& function,
                      AngularFrequencyCalculator const& calculator,
                      PoissonSeries<double, wdegree_, Evaluator> const& weight,
                      Instant const& t_min,
                      Instant const& t_max,
                      std::string_view const celestial,
                      mathematica::Logger& logger) {
  using Value = std::invoke_result_t<Function, Instant>;
  using Norm = typename Hilbert<Value>::NormType;
  using Norm² = typename Hilbert<Value>::InnerProductType;
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

  UnboundedLowerTriangularMatrix<Inverse<Norm>> α(basis_size, uninitialized);

  // Only indices 0 to m - 1 are used in this array.  At the beginning of
  // iteration m it contains Aⱼ⁽ᵐ⁻¹⁾.
  UnboundedVector<double> A(basis_size, uninitialized);

  Norm² const F₀ = InnerProduct(function, basis[0], weight, t_min, t_max);
  Norm² const Q₀₀ = InnerProduct(basis[0], basis[0], weight, t_min, t_max);
  α[0][0] = 1 / Sqrt(Q₀₀);
  A[0] = F₀ / Q₀₀;

  // At the beginning of iteration m this contains fₘ₋₁.
  auto f = function - A[0] * basis[0];

  int m_begin = 1;
  for (;;) {
    for (int m = m_begin; m < basis_size; ++m) {
      // Contains Fₘ.
      Norm² const F = InnerProduct(f, basis[m], weight, t_min, t_max);

      // This vector contains Qₘⱼ.
      UnboundedVector<Norm²> Q(m + 1, uninitialized);
      for (int j = 0; j <= m; ++j) {
        Q[j] = InnerProduct(basis[m], basis[j], weight, t_min, t_max);
      }

      // This vector contains Bⱼ⁽ᵐ⁾.
      UnboundedVector<Norm> B(m, uninitialized);
      for (int j = 0; j < m; ++j) {
        Norm Σ_αⱼₛ_Qₘₛ{};
        for (int s = 0; s <= j; ++s) {
          Σ_αⱼₛ_Qₘₛ += α[j][s] * Q[s];
        }
        B[j] = -Σ_αⱼₛ_Qₘₛ;
      }

      {
        Norm² Σ_Bₛ⁽ᵐ⁾²{};
        for (int s = 0; s < m; ++s) {
          Σ_Bₛ⁽ᵐ⁾² += B[s] * B[s];
        }
        const int max_degree = 2;
        double const sin2 = (Q[m] - Σ_Bₛ⁽ᵐ⁾²) / Q[m];
        if (sin2 <= 0/* ||
            m - m_begin >= 2 * 3 * (max_degree + 1)*/) {
          //TODO(phl):Comment.
          // We arrive here when the norm of Σₛ Bₛ⁽ᵐ⁾bₛ + eₘ is small (see
          // [SN97] for the notation) and, due to rounding errors, the computed
          // value of the square of that norm ends up negative, zero, or very
          // small.  It makes no sense to have complex numbers (or infinities)
          // here because our function is real and bounded.  But even if the
          // norm could be computed but was very small, we would end up with an
          // ill-conditioned solution.  Geometrically, we are in a situation
          // where eₘ is very close to the space spanned by the (bₛ), that is,
          // by the (eₛ) for i < m.  The fact that the basis elements are no
          // longer independent when the degree increases is duely noted by
          // [CV84].  Given that eₘ effectively doesn't have benefit for the
          // projection, we just drop it and continue with the algorithm.
          //LOG(ERROR) << "Q[m]: " << Q[m] << " Σ_Bₛ⁽ᵐ⁾² " << Σ_Bₛ⁽ᵐ⁾²
          //           << " difference: " << Q[m] - Σ_Bₛ⁽ᵐ⁾²;
          LOG(ERROR) << "Dropping basis element " << m << " due to sin=" << sin2
                     << " element=" << basis[m];
          int const basis_remaining = basis_size - m - 1;
          basis.erase(basis.begin() + m);
          α.EraseToEnd(m);
          α.Extend(basis_remaining, uninitialized);
          A.EraseToEnd(m);
          A.Extend(basis_remaining, uninitialized);
          --basis_size;
          --m;
          continue;
        } else {
          α[m][m] = 1 / Sqrt(Q[m] - Σ_Bₛ⁽ᵐ⁾²);
        }
      }

      for (int j = 0; j < m; ++j) {
        double Σ_Bₛ⁽ᵐ⁾_αₛⱼ = 0;
        for (int s = j; s < m; ++s) {
          Σ_Bₛ⁽ᵐ⁾_αₛⱼ += B[s] * α[s][j];
        }
        α[m][j] = α[m][m] * Σ_Bₛ⁽ᵐ⁾_αₛⱼ;
      }

      PoissonSeries<Normalized, degree_, Evaluator> Σ_αₘᵢ_eᵢ =
          α[m][0] * basis[0];
      for (int i = 1; i <= m; ++i) {
        Σ_αₘᵢ_eᵢ += α[m][i] * basis[i];
      }

      auto const norm2 = InnerProduct(Σ_αₘᵢ_eᵢ, Σ_αₘᵢ_eᵢ, weight, t_min, t_max);
      if (Abs(norm2 - 1) > 0x1.0p-8) {
        LOG(ERROR) << "Dropping basis element " << m << " due to norm=" << norm2
                   << " element=" << basis[m];
        int const basis_remaining = basis_size - m - 1;
        basis.erase(basis.begin() + m);
        α.EraseToEnd(m);
        α.Extend(basis_remaining, uninitialized);
        A.EraseToEnd(m);
        A.Extend(basis_remaining, uninitialized);
        --basis_size;
        --m;
        continue;
      }

      A[m] = α[m][m] * α[m][m] * F;

      for (int j = 0; j < m; ++j) {
        A[j] += α[m][m] * α[m][j] * F;
      }

      {
        PoissonSeries<Value, degree_, Evaluator> result = A[0] * basis[0];
        for (int i = 1; i < basis_size; ++i) {
          result += A[i] * basis[i];
        }
        logger.Append(absl::StrCat("intermediateSolution", celestial),
                      result,
                      mathematica::ExpressIn(Metre, Second, Radian));
      }

      f -= α[m][m] * F * Σ_αₘᵢ_eᵢ;
    }

    PoissonSeries<Value, degree_, Evaluator> result = A[0] * basis[0];
    for (int i = 1; i < basis_size; ++i) {
      result += A[i] * basis[i];
    }
    //logger.Append(absl::StrCat("residual", celestial),
    //              f,
    //              mathematica::ExpressIn(Metre, Second, Radian));
    logger.Append(absl::StrCat("solution", celestial),
                  result,
                  mathematica::ExpressIn(Metre, Second, Radian));

    ω = calculator(f);
    if (!ω.has_value()) {
      return result;
    }

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
    α.Extend(ω_basis_size, uninitialized);
    A.Extend(ω_basis_size, uninitialized);
    m_begin = basis_size;
    basis_size += ω_basis_size;
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
