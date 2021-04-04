﻿#pragma once

#include "astronomy/orbital_elements.hpp"

#include <algorithm>
#include <vector>

#include "base/jthread.hpp"
#include "physics/kepler_orbit.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbital_elements {

using base::Error;
using base::Status;
using base::this_stoppable_thread;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using physics::KeplerOrbit;
using quantities::ArcTan;
using quantities::Cos;
using quantities::Mod;
using quantities::Pow;
using quantities::Product;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Square;
using quantities::Tan;
using quantities::UnwindFrom;
using quantities::si::Radian;

template<typename PrimaryCentred>
StatusOr<OrbitalElements> OrbitalElements::ForTrajectory(
    DiscreteTrajectory<PrimaryCentred> const& trajectory,
    MassiveBody const& primary,
    Body const& secondary) {
  OrbitalElements orbital_elements;
  if (trajectory.Size() < 2) {
    return Status(Error::INVALID_ARGUMENT,
                  "trajectory.Size() is " + std::to_string(trajectory.Size()));
  }
  orbital_elements.osculating_equinoctial_elements_ =
      OsculatingEquinoctialElements(trajectory, primary, secondary);
  orbital_elements.radial_distances_ = RadialDistances(trajectory);
  auto const sidereal_period =
      SiderealPeriod(orbital_elements.osculating_equinoctial_elements_);
  RETURN_IF_ERROR(sidereal_period);
  orbital_elements.sidereal_period_ = sidereal_period.ValueOrDie();
  if (!IsFinite(orbital_elements.sidereal_period_) ||
      orbital_elements.sidereal_period_ <= Time{}) {
    // Guard against NaN sidereal periods (from hyperbolic orbits) or negative
    // sidereal periods (from aberrant trajectories, see #2811).
    return Status(
        Error::OUT_OF_RANGE,
        "sidereal period is " + DebugString(orbital_elements.sidereal_period_));
  }
  auto mean_equinoctial_elements =
      MeanEquinoctialElements(orbital_elements.osculating_equinoctial_elements_,
                              orbital_elements.sidereal_period_);
  RETURN_IF_ERROR(mean_equinoctial_elements);
  orbital_elements.mean_equinoctial_elements_ =
      std::move(mean_equinoctial_elements).ValueOrDie();
  if (orbital_elements.mean_equinoctial_elements_.size() < 2) {
    return Status(
        Error::OUT_OF_RANGE,
        "trajectory does not span one sidereal period: sidereal period is " +
            DebugString(orbital_elements.sidereal_period_) +
            ", trajectory spans " +
            DebugString(trajectory.back().time -
                        trajectory.front().time));
  }
  auto mean_classical_elements =
      ToClassicalElements(orbital_elements.mean_equinoctial_elements_);
  RETURN_IF_ERROR(mean_classical_elements);
  orbital_elements.mean_classical_elements_ =
      std::move(mean_classical_elements).ValueOrDie();
  RETURN_IF_ERROR(orbital_elements.ComputePeriodsAndPrecession());
  RETURN_IF_ERROR(orbital_elements.ComputeIntervals());
  return orbital_elements;
}

inline std::vector<OrbitalElements::ClassicalElements> const&
OrbitalElements::mean_elements() const {
  return mean_classical_elements_;
}

inline Time OrbitalElements::sidereal_period() const {
  return sidereal_period_;
}

inline Time OrbitalElements::nodal_period() const {
  return nodal_period_;
}

inline Time OrbitalElements::anomalistic_period() const {
  return anomalistic_period_;
}

inline AngularFrequency OrbitalElements::nodal_precession() const {
  return nodal_precession_;
}

inline Interval<Length> OrbitalElements::mean_semimajor_axis_interval() const {
  return mean_semimajor_axis_interval_;
}

inline Interval<double> OrbitalElements::mean_eccentricity_interval() const {
  return mean_eccentricity_interval_;
}

inline Interval<Angle> OrbitalElements::mean_inclination_interval() const {
  return mean_inclination_interval_;
}

inline Interval<Angle>
OrbitalElements::mean_longitude_of_ascending_node_interval() const {
  return mean_longitude_of_ascending_node_interval_;
}

inline Interval<Angle> OrbitalElements::mean_argument_of_periapsis_interval()
    const {
  return mean_argument_of_periapsis_interval_;
}

inline Interval<Length> OrbitalElements::mean_periapsis_distance_interval()
    const {
  return mean_periapsis_distance_interval_;
}

inline Interval<Length> OrbitalElements::mean_apoapsis_distance_interval()
    const {
  return mean_apoapsis_distance_interval_;
}

inline Interval<Length> OrbitalElements::radial_distance_interval() const {
  return radial_distance_interval_;
}

template<typename PrimaryCentred>
std::vector<OrbitalElements::EquinoctialElements>
OrbitalElements::OsculatingEquinoctialElements(
    DiscreteTrajectory<PrimaryCentred> const& trajectory,
    MassiveBody const& primary,
    Body const& secondary) {
  DegreesOfFreedom<PrimaryCentred> const primary_dof{
      PrimaryCentred::origin, PrimaryCentred::unmoving};
  std::vector<EquinoctialElements> result;
  for (auto const& [time, degrees_of_freedom] : trajectory) {
    auto const osculating_elements =
        KeplerOrbit<PrimaryCentred>(primary,
                                    secondary,
                                    degrees_of_freedom - primary_dof,
                                    time)
            .elements_at_epoch();
    double const& e = *osculating_elements.eccentricity;
    Angle const& ϖ = *osculating_elements.longitude_of_periapsis;
    Angle const& Ω = osculating_elements.longitude_of_ascending_node;
    Angle const& M = *osculating_elements.mean_anomaly;
    Angle const& i = osculating_elements.inclination;
    double const tg_½i = Tan(i / 2);
    double const cotg_½i = 1 / tg_½i;
    result.push_back(
        {.t = time,
         .a = *osculating_elements.semimajor_axis,
         .h = e * Sin(ϖ),
         .k = e * Cos(ϖ),
         .λ = result.empty() ? ϖ + M : UnwindFrom(result.back().λ, ϖ + M),
         .p = tg_½i * Sin(Ω),
         .q = tg_½i * Cos(Ω),
         .pʹ = cotg_½i * Sin(Ω),
         .qʹ = cotg_½i * Cos(Ω)});
  }
  return result;
}

template<typename PrimaryCentred>
std::vector<Length> OrbitalElements::RadialDistances(
    DiscreteTrajectory<PrimaryCentred> const& trajectory) {
  std::vector<Length> radial_distances;
  radial_distances.reserve(trajectory.Size());
  DegreesOfFreedom<PrimaryCentred> const primary_dof{PrimaryCentred::origin,
                                                     PrimaryCentred::unmoving};
  for (auto const& [time, degrees_of_freedom] : trajectory) {
    radial_distances.push_back(
        (degrees_of_freedom.position() - primary_dof.position()).Norm());
  }
  return radial_distances;
}

inline std::vector<OrbitalElements::EquinoctialElements> const&
OrbitalElements::osculating_equinoctial_elements() const {
  return osculating_equinoctial_elements_;
}

inline std::vector<OrbitalElements::EquinoctialElements> const&
OrbitalElements::mean_equinoctial_elements() const {
  return mean_equinoctial_elements_;
}

inline StatusOr<Time> OrbitalElements::SiderealPeriod(
    std::vector<EquinoctialElements> const& equinoctial_elements) {
  Time const Δt =
      equinoctial_elements.back().t - equinoctial_elements.front().t;
  Instant const t0 = equinoctial_elements.front().t + Δt / 2;
  Product<Angle, Square<Time>> ʃ_λt_dt;

  auto const first = equinoctial_elements.begin();
  Time const first_t = first->t - t0;
  Product<Angle, Time> previous_λt = first->λ * first_t;
  for (auto previous = first, it = first + 1;
       it != equinoctial_elements.end();
       previous = it, ++it) {
    RETURN_IF_STOPPED;
    Time const t = it->t - t0;
    auto const λt = it->λ * t;
    Time const dt = it->t - previous->t;
    ʃ_λt_dt += (λt + previous_λt) / 2 * dt;
    previous_λt = λt;
  }
  return 2 * π * Radian * Pow<3>(Δt) / (12 * ʃ_λt_dt);
}

inline StatusOr<std::vector<OrbitalElements::EquinoctialElements>>
OrbitalElements::MeanEquinoctialElements(
    std::vector<EquinoctialElements> const& osculating,
    Time const& period) {
  Instant const& t_min = osculating.front().t;
  // This function averages the elements in |osculating| over |period|.
  // For each |EquinoctialElements osculating_elements = osculating[i]| in
  // |osculating| such that |osculating_elements.t <= osculating.back().t|, let
  // |tᵢ = osculating_elements.t|; the result contains an
  // |EquinoctialElements mean_elements| with
  // |mean_elements.t = tᵢ + period / 2|.
  // For all э in the set of equinoctial elements {a, h, k, λ, p, q, pʹ, qʹ},
  // |mean_elements.э| is the integral of the osculating э from |tᵢ| to
  // |tᵢ + period|, divided by |period|.

  // Instead of computing the integral from |tᵢ| to |tᵢ + period| directly, we
  // precompute the integrals from |t_min| to each of the |tᵢ|.
  // The integral from |tᵢ| to |tᵢ + period| is then computed as the integral
  // from |tᵢ| to the last |tⱼ₋₁| before |tᵢ + period| (obtained by subtracting
  // the integrals to |tᵢ| and to |tⱼ₋₁|), plus the remainder integral from
  // |tⱼ₋₁| to |tᵢ + period| (a partial trapezoid on [tⱼ₋₁, tⱼ]).

  // Precompute the integrals from |t_min| to each of the |tᵢ|.
  struct IntegratedEquinoctialElements {
    Instant t_max;
    // The integrals are from t_min to t_max.
    Product<Length, Time> ʃ_a_dt;
    Time ʃ_h_dt;
    Time ʃ_k_dt;
    Product<Angle, Time> ʃ_λ_dt;
    Time ʃ_p_dt;
    Time ʃ_q_dt;
    Time ʃ_pʹ_dt;
    Time ʃ_qʹ_dt;
  };
  std::vector<IntegratedEquinoctialElements> integrals;
  integrals.push_back({t_min});
  for (auto previous = osculating.begin(), it = osculating.begin() + 1;
       it != osculating.end();
       previous = it, ++it) {
    RETURN_IF_STOPPED;
    integrals.push_back(integrals.back());
    integrals.back().t_max = it->t;
    Time const dt = it->t - previous->t;
    integrals.back().ʃ_a_dt += (it->a + previous->a) / 2 * dt;
    integrals.back().ʃ_h_dt += (it->h + previous->h) / 2 * dt;
    integrals.back().ʃ_k_dt += (it->k + previous->k) / 2 * dt;
    integrals.back().ʃ_λ_dt += (it->λ + previous->λ) / 2 * dt;
    integrals.back().ʃ_p_dt += (it->p + previous->p) / 2 * dt;
    integrals.back().ʃ_q_dt += (it->q + previous->q) / 2 * dt;
    integrals.back().ʃ_pʹ_dt += (it->pʹ + previous->pʹ) / 2 * dt;
    integrals.back().ʃ_qʹ_dt += (it->qʹ + previous->qʹ) / 2 * dt;
  }

  // Now compute the averages.
  std::vector<EquinoctialElements> mean_elements;
  int j = 0;
  for (auto const& up_to_tᵢ : integrals) {
    RETURN_IF_STOPPED;
    // We are averaging the elements over the interval [tᵢ, tᵢ + period].
    Instant const tᵢ = up_to_tᵢ.t_max;
    while (integrals[j].t_max < tᵢ + period) {
      if (++j == integrals.size()) {
        return mean_elements;
      }
    }
    // We have tⱼ₋₁ < tᵢ + period ≤ tⱼ.
    auto const& tⱼ = osculating[j].t;
    auto const& tⱼ₋₁ = osculating[j - 1].t;

    auto const& up_to_tⱼ₋₁ = integrals[j - 1];
    // |element| should be a pointer to a member of |EquinoctialElements|;
    // Integrates that element on [tⱼ₋₁, tᵢ + period].
    auto ʃ = [j, &period, &tᵢ, &tⱼ, &tⱼ₋₁, &osculating](auto element) {
      Time const Δt = tⱼ - tⱼ₋₁;
      Time const dt = tᵢ + period - tⱼ₋₁;
      auto const element_at_end =
          osculating[j - 1].*element +
          (osculating[j].*element - osculating[j - 1].*element) * (dt / Δt);
      return (osculating[j - 1].*element + element_at_end) / 2 * dt;
    };
    mean_elements.emplace_back();
    mean_elements.back().t = tᵢ + period / 2;
    mean_elements.back().a =
        (up_to_tⱼ₋₁.ʃ_a_dt - up_to_tᵢ.ʃ_a_dt + ʃ(&EquinoctialElements::a)) /
        period;
    mean_elements.back().h =
        (up_to_tⱼ₋₁.ʃ_h_dt - up_to_tᵢ.ʃ_h_dt + ʃ(&EquinoctialElements::h)) /
        period;
    mean_elements.back().k =
        (up_to_tⱼ₋₁.ʃ_k_dt - up_to_tᵢ.ʃ_k_dt + ʃ(&EquinoctialElements::k)) /
        period;
    mean_elements.back().λ =
        (up_to_tⱼ₋₁.ʃ_λ_dt - up_to_tᵢ.ʃ_λ_dt + ʃ(&EquinoctialElements::λ)) /
        period;
    mean_elements.back().p =
        (up_to_tⱼ₋₁.ʃ_p_dt - up_to_tᵢ.ʃ_p_dt + ʃ(&EquinoctialElements::p)) /
        period;
    mean_elements.back().q =
        (up_to_tⱼ₋₁.ʃ_q_dt - up_to_tᵢ.ʃ_q_dt + ʃ(&EquinoctialElements::q)) /
        period;
    mean_elements.back().pʹ =
        (up_to_tⱼ₋₁.ʃ_pʹ_dt - up_to_tᵢ.ʃ_pʹ_dt + ʃ(&EquinoctialElements::pʹ)) /
        period;
    mean_elements.back().qʹ =
        (up_to_tⱼ₋₁.ʃ_qʹ_dt - up_to_tᵢ.ʃ_qʹ_dt + ʃ(&EquinoctialElements::qʹ)) /
        period;
  }
  return mean_elements;
}

inline StatusOr<std::vector<OrbitalElements::ClassicalElements>>
OrbitalElements::ToClassicalElements(
    std::vector<EquinoctialElements> const& equinoctial_elements) {
  std::vector<ClassicalElements> classical_elements;
  for (auto const& equinoctial : equinoctial_elements) {
    RETURN_IF_STOPPED;
    double const tg_½i = Sqrt(Pow<2>(equinoctial.p) + Pow<2>(equinoctial.q));
    double const cotg_½i =
        Sqrt(Pow<2>(equinoctial.pʹ) + Pow<2>(equinoctial.qʹ));
    Angle const i =
        cotg_½i > tg_½i ? 2 * ArcTan(tg_½i) : 2 * ArcTan(1 / cotg_½i);
    Angle const Ω = cotg_½i > tg_½i ? ArcTan(equinoctial.p, equinoctial.q)
                                    : ArcTan(equinoctial.pʹ, equinoctial.qʹ);
    double const e = Sqrt(Pow<2>(equinoctial.h) + Pow<2>(equinoctial.k));
    Angle const ϖ = ArcTan(equinoctial.h, equinoctial.k);
    Angle const ω = ϖ - Ω;
    Angle const M = equinoctial.λ - ϖ;
    classical_elements.push_back(
        {.time = equinoctial.t,
         .semimajor_axis = equinoctial.a,
         .eccentricity = e,
         .inclination = i,
         .longitude_of_ascending_node = classical_elements.empty()
             ? Mod(Ω, 2 * π * Radian)
             : UnwindFrom(classical_elements.back().longitude_of_ascending_node,
                          Ω),
         .argument_of_periapsis = classical_elements.empty()
             ? Mod(ω, 2 * π * Radian)
             : UnwindFrom(classical_elements.back().argument_of_periapsis, ω),
         .mean_anomaly = classical_elements.empty()
             ? Mod(M, 2 * π * Radian)
             : UnwindFrom(classical_elements.back().mean_anomaly, M),
         .periapsis_distance = (1 - e) * equinoctial.a,
         .apoapsis_distance = (1 + e) * equinoctial.a});
  }
  return classical_elements;
}

inline Status OrbitalElements::ComputePeriodsAndPrecession() {
  Time const Δt = mean_classical_elements_.back().time -
                  mean_classical_elements_.front().time;
  auto const Δt³ = Pow<3>(Δt);
  // We compute the mean rate (slope) of the mean anomaly M(t), the mean
  // argument of latitude u(t), and the longitude of the ascending node Ω(t).
  // On an interval [t̄ - Δt/2, t̄ + Δt/2], the slope of э is computed as
  //   ∫ (э(t) - э̄) (t - t̄) dt / ∫ (t - t̄)² dt;
  // this is the continuous analogue of a simple linear regression.
  // With ∫ (t - t̄)² dt = Δt³ / 12 and ∫ э̄ (t - t̄) = э̄ ∫ (t - t̄) dt = 0,
  // this simplifies to
  //   12 ∫ э(t) (t - t̄) dt / Δt³.
  // We first compute ∫ э(t) (t - t̄) dt for the three elements of interest.

  Product<Angle, Square<Time>> ʃ_Mt_dt;
  Product<Angle, Square<Time>> ʃ_ut_dt;
  Product<Angle, Square<Time>> ʃ_Ωt_dt;

  Instant const t̄ = mean_classical_elements_.front().time + Δt / 2;
  auto const first = mean_classical_elements_.begin();
  Time first_t = first->time - t̄;
  Product<Angle, Time> previous_Mt = first->mean_anomaly * first_t;
  Product<Angle, Time> previous_ut =
      (first->argument_of_periapsis + first->mean_anomaly) * first_t;
  Product<Angle, Time> previous_Ωt =
      first->longitude_of_ascending_node * first_t;
  for (auto previous = first, it = first + 1;
       it != mean_classical_elements_.end();
       previous = it, ++it) {
    RETURN_IF_STOPPED;
    Angle const u = it->argument_of_periapsis + it->mean_anomaly;
    Angle const& M = it->mean_anomaly;
    Angle const& Ω = it->longitude_of_ascending_node;
    Time const t = it->time - t̄;
    auto const Mt = M * t;
    auto const ut = u * t;
    auto const Ωt = Ω * t;
    Time const dt = it->time - previous->time;
    ʃ_Mt_dt += (Mt + previous_Mt) / 2 * dt;
    ʃ_ut_dt += (ut + previous_ut) / 2 * dt;
    ʃ_Ωt_dt += (Ωt + previous_Ωt) / 2 * dt;
    previous_Mt = Mt;
    previous_ut = ut;
    previous_Ωt = Ωt;
  }

  // The periods are 2π over the mean rate of the relevant element; the nodal
  // precession is the mean rate of Ω.

  anomalistic_period_ = 2 * π * Radian * Δt³ / (12 * ʃ_Mt_dt);
  nodal_period_ = 2 * π * Radian * Δt³ / (12 * ʃ_ut_dt);
  nodal_precession_ = 12 * ʃ_Ωt_dt / Δt³;
  return Status::OK;
}

inline Status OrbitalElements::ComputeIntervals() {
  for (auto const& r : radial_distances_) {
    RETURN_IF_STOPPED;
    radial_distance_interval_.Include(r);
  }
  for (auto const& elements : mean_classical_elements_) {
    RETURN_IF_STOPPED;
    mean_semimajor_axis_interval_.Include(elements.semimajor_axis);
    mean_eccentricity_interval_.Include(elements.eccentricity);
    mean_inclination_interval_.Include(elements.inclination);
    mean_longitude_of_ascending_node_interval_.Include(
        elements.longitude_of_ascending_node);
    mean_argument_of_periapsis_interval_.Include(
        elements.argument_of_periapsis);
    mean_periapsis_distance_interval_.Include(elements.periapsis_distance);
    mean_apoapsis_distance_interval_.Include(elements.apoapsis_distance);
  }
  return Status::OK;
}


}  // namespace internal_orbital_elements
}  // namespace astronomy
}  // namespace principia
