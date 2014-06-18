﻿#pragma once

#include <memory>
#include <vector>

#include "integrators/symplectic_integrator.hpp"
#include "physics/body.hpp"
#include "quantities/quantities.hpp"

using principia::integrators::SymplecticIntegrator;
using principia::quantities::Acceleration;
using principia::quantities::Length;
using principia::quantities::Speed;
using principia::quantities::Time;

namespace principia {
namespace physics {

template<typename InertialFrame>
class NBodySystem {
 public:
  typedef std::vector<std::unique_ptr<Body<InertialFrame>>> Bodies;
  typedef SymplecticIntegrator<Length, Speed> Integrator;

  // NOTE(phl): We would prefer to pass the unique_ptr<> by value, but that
  // confuses the compiler.  So for now, we'll use r-value references.
  NBodySystem(std::unique_ptr<Bodies>&& massive_bodies,
              std::unique_ptr<Bodies>&& massless_bodies);
  ~NBodySystem() = default;

  // No transfer of ownership.
  std::vector<Body<InertialFrame> const*> massive_bodies() const;
  std::vector<Body<InertialFrame> const*> massless_bodies() const;
  std::vector<Body<InertialFrame> const*> bodies() const;

  // The |integrator| must already have been initialized.  Performs an
  // integration for all the bodies passed at construction.
  //TODO(phl): From which time?
  void IntegrateFull(Integrator const& integrator,
                     Time const& tmax,
                     Time const& Δt,
                     int const sampling_period);

  // The |integrator| must already have been initialized.  Performs an
  // integration for only the given |massless_bodies| (which may be a subset of
  // the ones passed at construction).  All the other bodies are assumed to
  // remain stationary.
  void IntegrateSubset(std::vector<Body<InertialFrame>> const& massless_bodies,
                       Integrator const& integrator,
                       Time const& tmax,
                       Time const& Δt,
                       int const sampling_period);

 private:
  void ComputeGravitationalAccelerations(
      Time const& t,
      std::vector<Length> const& q,
      std::vector<Acceleration>* result) const;
  static void ComputeGravitationalVelocities(std::vector<Speed> const& p,
                                             std::vector<Speed>* result);

  std::unique_ptr<Bodies const> const massive_bodies_;  // Never null.
  std::unique_ptr<Bodies const> const massless_bodies_;  // Never null.

  // The pointers are not owned.  The massive bodies come first.
  std::vector<Body<InertialFrame>*> bodies_;
};

}  // namespace physics
}  // namespace principia

#include "physics/n_body_system_body.hpp"
