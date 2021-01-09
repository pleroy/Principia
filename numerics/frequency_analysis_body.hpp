
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

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using base::uninitialized;
using geometry::Hilbert;
using geometry::Vector;
using quantities::Abs;
using quantities::Inverse;
using quantities::IsFinite;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;

// Appends basis elements for |ω| to |basis| and |basis_subspaces|.  Returns the
// number of elements that were appended.
template<int aperiodic_degree, int periodic_degree,
         typename BasisSeries>
int MakeBasis(std::optional<AngularFrequency> const& ω,
              Instant const& t_min,
              Instant const& t_max,
              std::vector<BasisSeries>& basis,
              std::vector<PoissonSeriesSubspace>& basis_subspaces) {
  if (ω.value() == AngularFrequency{}) {
    using Generator =
        PoissonSeriesBasisGenerator<BasisSeries, aperiodic_degree>;
    auto const ω_basis = Generator::Basis(t_min, t_max);
    auto const ω_basis_subspaces = Generator::Subspaces();
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    std::move(ω_basis_subspaces.begin(),
              ω_basis_subspaces.end(),
              std::back_inserter(basis_subspaces));
    return std::tuple_size_v<decltype(ω_basis)>;
  } else {
    using Generator = PoissonSeriesBasisGenerator<BasisSeries, periodic_degree>;
    auto const ω_basis = Generator::Basis(ω.value(), t_min, t_max);
    auto const ω_basis_subspaces = Generator::Subspaces(ω.value());
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    std::move(ω_basis_subspaces.begin(),
              ω_basis_subspaces.end(),
              std::back_inserter(basis_subspaces));
    return std::tuple_size_v<decltype(ω_basis)>;
  }
}

template<typename Function,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree,
                  Evaluator> const& weight) {
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

template<int aperiodic_degree, int periodic_degree,
         typename Function,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>,
              aperiodic_degree, periodic_degree,
              Evaluator>
Projection(Function const& function,
           AngularFrequency const& ω,
           PoissonSeries<double,
                         aperiodic_wdegree, periodic_wdegree,
                         Evaluator> const& weight,
           Instant const& t_min,
           Instant const& t_max) {
  std::optional<AngularFrequency> optional_ω = ω;

  // A calculator that returns optional_ω once and then stops.
  auto angular_frequency_calculator = [&optional_ω](auto const& residual) {
    auto const result = optional_ω;
    optional_ω = std::nullopt;
    return result;
  };

  return IncrementalProjection<aperiodic_degree, periodic_degree>(
      function,
      angular_frequency_calculator,
      weight,
      t_min, t_max);
}

#define DO_THE_LOGGING 1
#define USE_CGS 0
#define USE_INTEGRATE2 1
#define USE_INTEGRATE3 0

template<int aperiodic_degree, int periodic_degree,
         typename Function,
         typename AngularFrequencyCalculator,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>,
              aperiodic_degree, periodic_degree,
              Evaluator>
IncrementalProjection(Function const& function,
                      AngularFrequencyCalculator const& calculator,
                      PoissonSeries<double,
                                    aperiodic_wdegree, periodic_wdegree,
                                    Evaluator> const& weight,
                      Instant const& t_min,
                      Instant const& t_max) {
  using Value = std::invoke_result_t<Function, Instant>;
  using Norm = typename Hilbert<Value>::NormType;
  using Normalized = typename Hilbert<Value>::NormalizedType;
  using BasisSeries = PoissonSeries<Normalized,
                                    aperiodic_degree, periodic_degree,
                                    Evaluator>;
  using ResultSeries = PoissonSeries<Value,
                                     aperiodic_degree, periodic_degree,
                                     Evaluator>;

  // This code follows [Kud07], section 2.  Our indices start at 0, unlike those
  // of Кудрявцев which start at 1.

  Instant const& t0 = weight.origin();

  std::optional<AngularFrequency> ω = calculator(function);
  CHECK(ω.has_value());

  std::vector<BasisSeries> basis;
  // The Poisson series basis[k] belongs to the subspace basis_subspaces[k];
  // this remains true after orthonormalization, i.e., q[k] belongs to the
  // subspace basis_subspaces[k] below.
  std::vector<PoissonSeriesSubspace> basis_subspaces;

  int basis_size = MakeBasis<aperiodic_degree, periodic_degree>(
      ω, t_min, t_max, basis, basis_subspaces);

  // This is logically Q in the QR decomposition of basis.
  std::vector<BasisSeries> q;

  auto const& a₀ = basis[0];
  auto const r₀₀ = a₀.Norm(weight, t_min, t_max);
  CHECK(IsFinite(r₀₀)) << a₀;
  q.push_back(a₀ / r₀₀);

  auto const A₀ = InnerProduct(function, q[0], weight, t_min, t_max);

  ResultSeries F = A₀ * q[0];
  auto f = function - F;

  auto do_the_logging = [&](int m, BasisSeries const& qm) {
    double max = 0;
    double max_orth = 0;
    int i_max = 0;
    int i_max_orth = 0;
    for (int i = 0; i < q.size(); ++i) {
#if USE_INTEGRATE2
      double ip = (PointwiseInnerProduct(q[i], qm) * weight)
                      .Integrate(t_min, t_max) /
                  (t_max - t_min);
#else
      double ip = InnerProduct(q[i], qm, weight, t_min, t_max);
#endif
      if (i == m) {
        ip -= 1;
      }
      if (PoissonSeriesSubspace::orthogonal(basis_subspaces[i],
                                            basis_subspaces[m])) {
        if (Abs(ip) > max_orth) {
          max_orth = Abs(ip);
          i_max_orth = i;
        }
      } else {
        if (Abs(ip) > max) {
          max = Abs(ip);
          i_max = i;
        }
      }
    }
    LOG(ERROR) << "Max err: " << max << " at: " << i_max
               << " max orth err: " << max_orth << " at: " << i_max_orth;
  };
  int m_begin = 1;
  for (;;) {
    for (int m = m_begin; m < basis_size; ++m) {
      auto const& aₘ = basis[m];

      // This code follows Björk, Numerics of Gram-Schmidt Orthogonalization,
      // Algorithm 6.1.
      static constexpr double α = 0.5;
      BasisSeries q̂ₘ = aₘ;

      // Formal integration works for a single basis element.
      double q̂ₘ_norm = Sqrt(
          (PointwiseInnerProduct(q̂ₘ, q̂ₘ) * weight).Integrate(t_min, t_max) /
          (t_max - t_min));

      // Loop on p.
      BasisSeries previous_q̂ₘ = q̂ₘ;
      double previous_q̂ₘ_norm;
      do {
        previous_q̂ₘ = q̂ₘ;
        previous_q̂ₘ_norm = q̂ₘ_norm;
        for (int i = 0; i < m; ++i) {
          if (!PoissonSeriesSubspace::orthogonal(basis_subspaces[i],
                                                 basis_subspaces[m])) {
#if USE_INTEGRATE2
           auto const sᵖₘ = (PointwiseInnerProduct(q[i], previous_q̂ₘ) * weight)
                                .Integrate(t_min, t_max) /
                            (t_max - t_min);
#else
           auto const sᵖₘ =
                InnerProduct(q[i], previous_q̂ₘ, weight, t_min, t_max);
#endif
            q̂ₘ -= sᵖₘ * q[i];
          }
        }
#if USE_INTEGRATE3
        q̂ₘ_norm = Sqrt(
            (PointwiseInnerProduct(q̂ₘ, q̂ₘ) * weight).Integrate(t_min, t_max) /
            (t_max - t_min));
#else
        q̂ₘ_norm = q̂ₘ.Norm(weight, t_min, t_max);
#endif
#if DO_THE_LOGGING
        do_the_logging(m, q̂ₘ / q̂ₘ_norm);
        LOG(ERROR) << "m: " << m << " previous q̂ₘ: " << previous_q̂ₘ_norm
                   << " q̂ₘ: " << q̂ₘ_norm;
#endif

        if (!IsFinite(q̂ₘ_norm)) {
          return F;
        }
      } while (q̂ₘ_norm < α * previous_q̂ₘ_norm);

      q.push_back(q̂ₘ / q̂ₘ_norm);
#if DO_THE_LOGGING
      do_the_logging(m, q[m]);
#endif
#if USE_CGS
      auto const r₀ₘ = InnerProduct(q[0], aₘ, weight, t_min, t_max);
      Series Σrᵢₘqᵢ = r₀ₘ * q[0];
      for (int i = 1; i < m; ++i) {
        auto const rᵢₘ = InnerProduct(q[i], aₘ, weight, t_min, t_max);
        Σrᵢₘqᵢ += rᵢₘ * q[i];
      }

      auto const qʹₘ = aₘ - Σrᵢₘqᵢ;
      auto rₘₘ = qʹₘ.Norm(weight, t_min, t_max);
      q.push_back(qʹₘ / rₘₘ);
#endif
      DCHECK_EQ(m + 1, q.size());

      Norm const Aₘ = InnerProduct(f, q[m], weight, t_min, t_max);

      auto const Aₘqₘ = Aₘ * q[m];
      f -= Aₘqₘ;
      F += Aₘqₘ;
    }

    ω = calculator(f);
    if (!ω.has_value()) {
      return F;
    }

    int const ω_basis_size = MakeBasis<aperiodic_degree, periodic_degree>(
        ω, t_min, t_max, basis, basis_subspaces);
    m_begin = basis_size;
    basis_size += ω_basis_size;
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
