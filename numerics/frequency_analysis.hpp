﻿
#pragma once

#include <string_view>
#include <type_traits>

#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/poisson_series.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using geometry::Instant;
using geometry::Interval;
using quantities::AngularFrequency;
using quantities::Primitive;
using quantities::Product;
using quantities::Time;

// Computes the precise mode of a quasi-periodic |function|, assuming that the
// mode is over the interval |fft_mode| (so named because it has presumably been
// obtained using FFT).  The |Function| must have a member |FourierTransform|
// that returns its spectrum.  See [Cha95].
template<typename Function,
         int wdegree_,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight);

// In the projection functions the |Function| must have an |InnerProduct| with
// |PoissonSeries| or |PiecewisePoissonSeries|.

// Computes the Кудрявцев projection of |function| on a basis with angular
// frequency ω and maximum degree |degree_|.  See [Kud07].
template<int degree_,
         typename Function,
         int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
Projection(Function const& function,
           AngularFrequency const& ω,
           PoissonSeries<double, wdegree_, Evaluator> const& weight,
           Instant const& t_min,
           Instant const& t_max);

// AngularFrequencyCalculator is a templated functor that implements the
// extraction of the most relevant frequency out of a (mostly periodic)
// residual.  Its declaration must look like:
//
// class AngularFrequencyCalculator {
//  public:
//   ...
//   template<typename Residual>
//   std::optional<AngularFrequency>
//   operator()(Residual const& residual) const;
//   ...
//};
//
// Where Residual is a functor that takes an Instant and returns an element of a
// vector space.  The first call to the calculator is with the |function| passed
// to |IncrementalProjection|.
// If the calculator cannot find a suitable frequency, or if it wants to stop
// the algorithm, it does so by returning std::nullopt.
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
                      std::string_view celestial,
                      mathematica::Logger& logger);

}  // namespace internal_frequency_analysis

using internal_frequency_analysis::IncrementalProjection;
using internal_frequency_analysis::PreciseMode;
using internal_frequency_analysis::Projection;

}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia

#include "numerics/frequency_analysis_body.hpp"
