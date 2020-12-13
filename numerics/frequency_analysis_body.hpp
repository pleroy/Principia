
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
using quantities::Inverse;
using quantities::IsFinite;
using quantities::Pow;
using quantities::Quotient;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;

//TODO(phl): Make all this incremental.
template<typename Scalar>
UnboundedUpperTriangularMatrix<SquareRoot<Scalar>> CholeskyFactorization(
    UnboundedUpperTriangularMatrix<Scalar> const& a) {
  UnboundedUpperTriangularMatrix<SquareRoot<Scalar>> r(a.columns(),
                                                       uninitialized);
  for (int j = 0; j < a.columns(); ++j) {
    for (int i = 0; i < j; ++i) {
      Scalar Σrkirkj{};  //TODO(phl):Unicode.
      for (int k = 0; k < i; ++k) {
        Σrkirkj += r[k][i] * r[k][j];
      }
      r[i][j] = (a[i][j] - Σrkirkj) / r[i][i];
    }
    Scalar Σrkj2{};  //TODO(phl):Unicode.
    for (int k = 0; k < j; ++k) {
      Σrkj2 += Pow<2>(r[k][j]);
    }
    r[j][j] = Sqrt(a[j][j] - Σrkj2);  //TODO(phl):NaN?
  }
  return r;
}

//TODO(phl): Krishnamoorthy
template<typename Scalar>
void RDRDecomposition(  // TODO(phl):Unicode.
    UnboundedUpperTriangularMatrix<Scalar> const& a,
    UnboundedUpperTriangularMatrix<double>& r,
    UnboundedVector<Scalar>& d) {
  for (int i = 0; i < a.columns(); ++i) {
    Scalar Σrki2dk{};  //TODO(phl):Unicode.
    for (int k = 0; k < i; ++k) {
      Σrki2dk += Pow<2>(r[k][i]) * d[k];
    }
    d[i] = a[i][i] - Σrki2dk;
    for (int j = i + 1; j < a.columns(); ++j) {
      Scalar Σrkirkjdk{};  //TODO(phl):Unicode.
      for (int k = 0; k < i; ++k) {
        Σrkirkjdk += r[k][i] * r[k][j] * d[k];
      }
      r[i][j] = (a[i][j] - Σrkirkjdk) / d[i];
    }
    r[i][i] = 1;
  }
}

template<typename LScalar, typename RScalar>
UnboundedVector<Quotient<RScalar, LScalar>> BackSubstitution(
    UnboundedUpperTriangularMatrix<LScalar> const& u,
    UnboundedVector<RScalar> const& b) {
  UnboundedVector<Quotient<RScalar, LScalar>> x(b.size(), uninitialized);
  int const n = b.size() - 1;
  x[n] = b[n] / u[n][n];
  for (int i = n - 1; i >= 0; --i) {
    auto s = b[i];
    for (int j = i + 1; j <= n; ++j) {
      s -= u[i][j] * x[j];
    }
    x[i] = s / u[i][i];
  }
  return x;
}

//TODO(phl): https://en.wikipedia.org/wiki/Triangular_matrix#Forward_substitution
template<typename LScalar, typename RScalar>
UnboundedVector<Quotient<RScalar, LScalar>> ForwardSubstitution(
    UnboundedLowerTriangularMatrix<LScalar> const& l,
    UnboundedVector<RScalar> const& b) {
  UnboundedVector<Quotient<RScalar, LScalar>> x(b.size(), uninitialized);
  x[0] = b[0] / l[0][0];
  for (int i = 1; i < b.size(); ++i) {
    auto s = b[i];
    for (int j = 0; j < i; ++j) {
      s -= l[i][j] * x[j];
    }
    x[i] = s / l[i][i];
  }
  return x;
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
  using Norm² = typename Hilbert<Value>::Norm²Type;
  using Normalized = typename Hilbert<Value>::NormalizedType;
  using Series = PoissonSeries<Value,
                               aperiodic_degree, periodic_degree,
                               Evaluator>;

  // This code follows [Kud07], section 2.  Our indices start at 0, unlike those
  // of Кудрявцев which start at 1.

  Instant const& t0 = weight.origin();

  std::optional<AngularFrequency> ω = calculator(function);
  CHECK(ω.has_value());

  std::vector<Series> basis;
  // The Poisson series basis[k] belongs to the subspace basis_subspaces[k];
  // this remains true after orthonormalization, i.e., q[k] belongs to the
  // subspace basis_subspaces[k] below.
  std::vector<PoissonSeriesSubspace> basis_subspaces;

  int basis_size;
  // TODO(phl): This is replicated below.
  if (ω.value() == AngularFrequency{}) {
    auto const ω_basis =
        PoissonSeriesBasisGenerator<Series,
                                    aperiodic_degree>::Basis(t0);
    auto const ω_basis_subspaces =
        PoissonSeriesBasisGenerator<Series, aperiodic_degree>::Subspaces(t0);
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    std::move(ω_basis_subspaces.begin(),
              ω_basis_subspaces.end(),
              std::back_inserter(basis_subspaces));
  } else {
    auto const ω_basis =
        PoissonSeriesBasisGenerator<Series,
                                    periodic_degree>::Basis(ω.value(), t0);
    auto const ω_basis_subspaces =
        PoissonSeriesBasisGenerator<Series, periodic_degree>::Subspaces(
            ω.value(), t0);
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    std::move(ω_basis_subspaces.begin(),
              ω_basis_subspaces.end(),
              std::back_inserter(basis_subspaces));
  }

  UnboundedUpperTriangularMatrix<Norm²> C(basis_size);  // Zero-initialized.
  UnboundedVector<Norm²> c(basis_size, uninitialized);

  int m_begin = 0;
  for (;;) {
    for (int m = m_begin; m < basis_size; ++m) {
      auto const& aₘ = basis[m];
      for (int k = 0; k <= m; ++k) {
        if (!PoissonSeriesSubspace::orthogonal(basis_subspaces[k],
                                               basis_subspaces[m])) {
          C[k][m] = InnerProduct(basis[k], aₘ, weight, t_min, t_max);
          //C[k][m] = (PointwiseInnerProduct(basis[k], aₘ) * weight)
          //              .Integrate(t_min, t_max) /
          //          (t_max - t_min);
        }
      }
      c[m] = InnerProduct(function, aₘ, weight, t_min, t_max);
    }
#define PRINCIPIA_USE_CHOLESKY 0
#if PRINCIPIA_USE_CHOLESKY
    UnboundedUpperTriangularMatrix<Norm> const R = CholeskyFactorization(C);
    UnboundedVector<Norm> const y =
        ForwardSubstitution(R.Transpose(), c);  //TODO(phl): Costly?
    auto const x = BackSubstitution(R, y);
#else
    UnboundedUpperTriangularMatrix<double> R(basis_size, uninitialized);
    UnboundedVector<Norm²> D(basis_size, uninitialized);
    RDRDecomposition(C, R, D);
    UnboundedVector<Norm²> y =
        ForwardSubstitution(R.Transpose(), c);  //TODO(phl): Costly?
    UnboundedVector<double> yDminus1(basis_size,
                                     uninitialized);  // TODO(phl): Unicode.
    for (int m = 0; m < y.size(); ++m) {
      LOG_IF(ERROR, D[m] < Norm²{})<<"D["<<m<<"]: "<<D[m];
      yDminus1[m] = y[m] / D[m];
    }
    auto const x = BackSubstitution(R, yDminus1);
#endif
#undef PRINCIPIA_USE_CHOLESKY

    auto f = function - x[0] * basis[0];
    Series F = x[0] * basis[0];
    for (int m = 1; m < x.size(); ++m) {
      //TODO(phl):factor
      f -= x[m] * basis[m];
      F += x[m] * basis[m];
    }

    ω = calculator(f);
    if (!ω.has_value()) {
      return F;
    }

    int ω_basis_size;
    if (ω.value() == AngularFrequency{}) {
      auto const ω_basis =
          PoissonSeriesBasisGenerator<Series,
                                      aperiodic_degree>::Basis(t0);
      auto const ω_basis_subspaces =
          PoissonSeriesBasisGenerator<Series, aperiodic_degree>::Subspaces(t0);
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
      std::move(ω_basis_subspaces.begin(),
                ω_basis_subspaces.end(),
                std::back_inserter(basis_subspaces));
    } else {
      auto const ω_basis =
          PoissonSeriesBasisGenerator<Series,
                                      periodic_degree>::Basis(ω.value(), t0);
      auto const ω_basis_subspaces =
          PoissonSeriesBasisGenerator<Series, periodic_degree>::Subspaces(
              ω.value(), t0);
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
      std::move(ω_basis_subspaces.begin(),
                ω_basis_subspaces.end(),
                std::back_inserter(basis_subspaces));
    }
    m_begin = basis_size;
    basis_size += ω_basis_size;
    C.Extend(ω_basis_size);
    c.Extend(ω_basis_size);
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
