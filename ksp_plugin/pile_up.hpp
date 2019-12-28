﻿
#pragma once

#include <functional>
#include <future>
#include <list>
#include <map>

#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/grassmann.hpp"
#include "integrators/integrators.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/euler_solver.hpp"
#include "physics/inertia_tensor.hpp"
#include "physics/massless_body.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {

FORWARD_DECLARE_FROM(part, class, Part);

namespace internal_pile_up {

using base::not_null;
using base::Status;
using geometry::Bivector;
using geometry::Frame;
using geometry::Instant;
using geometry::NonInertial;
using geometry::Vector;
using integrators::Integrator;
using physics::DiscreteTrajectory;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::EulerSolver;
using physics::InertiaTensor;
using physics::MasslessBody;
using physics::RelativeDegreesOfFreedom;
using quantities::AngularMomentum;
using quantities::Force;
using quantities::Mass;

// A |PileUp| handles a connected component of the graph of |Parts| under
// physical contact.  It advances the history and psychohistory of its component
// |Parts|, modeling them as a massless body at their centre of mass.
class PileUp {
 public:
  PileUp(std::list<not_null<Part*>>&& parts,
         Instant const& t,
         Ephemeris<Barycentric>::AdaptiveStepParameters const&
             adaptive_step_parameters,
         Ephemeris<Barycentric>::FixedStepParameters const&
             fixed_step_parameters,
         not_null<Ephemeris<Barycentric>*> ephemeris,
         std::function<void()> deletion_callback);

  virtual ~PileUp();

  // This class is moveable.
  PileUp(PileUp&& pile_up) = default;
  PileUp& operator=(PileUp&& pile_up) = default;

  std::list<not_null<Part*>> const& parts() const;

  // Set the |degrees_of_freedom| for the given |part|.  These degrees of
  // freedom are *apparent* in the sense that they were reported by the game but
  // we know better since we are doing science.
  void SetPartApparentDegreesOfFreedom(
      not_null<Part*> part,
      DegreesOfFreedom<ApparentBubble> const& degrees_of_freedom);

  // Deforms the pile-up, advances the time, and nudges the parts, in sequence.
  // Does nothing if the psychohistory is already advanced beyond |t|.  Several
  // executions of this method may happen concurrently on multiple threads, but
  // not concurrently with any other method of this class.
  Status DeformAndAdvanceTime(Instant const& t);

  // Recomputes the state of motion of the pile-up based on that of its parts.
  void RecomputeFromParts(Instant const& t);

  // We'd like to return |not_null<std::shared_ptr<PileUp> const&|, but the
  // compiler gets confused when defining the corresponding lambda, and thinks
  // that we return a local variable even though we capture by reference.
  // TODO(phl): Try to fix in VS2017 or later.
  using PileUpForSerializationIndex =
      std::function<not_null<std::shared_ptr<PileUp>>(int)>;
  using SerializationIndexForPileUp =
      std::function<int(not_null<PileUp const*>)>;

  void WriteToMessage(not_null<serialization::PileUp*> message) const;
  static not_null<std::unique_ptr<PileUp>> ReadFromMessage(
      serialization::PileUp const& message,
      std::function<not_null<Part*>(PartId)> const& part_id_to_part,
      not_null<Ephemeris<Barycentric>*> ephemeris,
      std::function<void()> deletion_callback);

 private:
  // A pointer to a member function of |Part| used to append a point to either
  // trajectory (history or psychohistory).
  using AppendToPartTrajectory =
      void (Part::*)(Instant const&, DegreesOfFreedom<Barycentric> const&);

  // The principal axes of the pile-up.
  using RigidPileUpPrincipalAxes = Frame<enum class RigidPileUpPrincipalAxesTag,
                                         NonInertial,
                                         RigidPileUp::handedness>;

  // For deserialization.
  PileUp(std::list<not_null<Part*>>&& parts,
         Ephemeris<Barycentric>::AdaptiveStepParameters const&
             adaptive_step_parameters,
         Ephemeris<Barycentric>::FixedStepParameters const&
             fixed_step_parameters,
         not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> history,
         DiscreteTrajectory<Barycentric>* psychohistory,
         not_null<Ephemeris<Barycentric>*> ephemeris,
         std::function<void()> deletion_callback);

  // Update the degrees of freedom (in |RigidPileUp|) of all the parts by
  // translating the *apparent* degrees of freedom so that their centre of mass
  // matches that computed by integration.
  // |SetPartApparentDegreesOfFreedom| must have been called for each part in
  // the pile-up, or for none.
  // The degrees of freedom set by this method are used by |NudgeParts|.
  // NOTE(egg): Eventually, this will also change their velocities and angular
  // velocities so that the angular momentum matches that which has been
  // computed for |this| |PileUp|.
  void DeformPileUpIfNeeded();

  // Flows the history authoritatively as far as possible up to |t|, advances
  // the histories of the parts and updates the degrees of freedom of the parts
  // if the pile-up is in the bubble.  After this call, the tail (of |*this|)
  // and of its parts have a (possibly ahistorical) final point exactly at |t|.
  Status AdvanceTime(Instant const& t);

  // Adjusts the degrees of freedom of all parts in this pile up based on the
  // degrees of freedom of the pile-up computed by |AdvanceTime| and on the
  // |RigidPileUp| degrees of freedom of the parts, as set by
  // |DeformPileUpIfNeeded|.
  void NudgeParts() const;

  template<AppendToPartTrajectory append_to_part_trajectory>
  void AppendToPart(DiscreteTrajectory<Barycentric>::Iterator it) const;

  // Computes the angular momentum, inertia tensor and intrinsic force from the
  // list of parts.  Returns the barycentre of the parts.
  DegreesOfFreedom<Barycentric> RecomputeFromParts(
      Instant const& t,
      std::list<not_null<Part*>> const& parts);

  // Wrapped in a |unique_ptr| to be moveable.
  not_null<std::unique_ptr<absl::Mutex>> lock_;

  std::list<not_null<Part*>> parts_;
  not_null<Ephemeris<Barycentric>*> ephemeris_;
  Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters_;
  Ephemeris<Barycentric>::FixedStepParameters fixed_step_parameters_;

  // Recomputed by the parts subset on every change.  Not serialized.
  Bivector<AngularMomentum, RigidPileUp> angular_momentum_;
  InertiaTensor<RigidPileUp> inertia_tensor_;
  Vector<Force, Barycentric> intrinsic_force_;
  std::unique_ptr<EulerSolver<RigidPileUp, RigidPileUpPrincipalAxes>>
      euler_solver_;

  // The |history_| is the past trajectory of the pile-up.  It is normally
  // integrated with a fixed step using |fixed_instance_|, except in the
  // presence of intrinsic acceleration.  It is authoritative in the sense that
  // it is never going to change.
  not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> history_;

  // The |psychohistory_| is the recent past trajectory of the pile-up.  Since
  // we need to draw something between the last point of the |history_| and the
  // current time, we must have a bit of trajectory that may not cover an entire
  // fixed step.  This part is the |psychohistory_|, and it is forked at the end
  // of the |history_|.  It is not authoritative in the sense that it may not
  // match the |history_| that we'll ultimately compute.  The name comes from
  // the fact that we are trying to predict the future, but since we are not as
  // good as Hari Seldon we only do it over a short period of time.
  DiscreteTrajectory<Barycentric>* psychohistory_ = nullptr;

  // When present, this instance is used to integrate the trajectory of this
  // pile-up using a fixed-step integrator.  This instance is destroyed
  // if a variable-step integrator needs to be used because of an intrinsic
  // acceleration.
  std::unique_ptr<typename Integrator<
      Ephemeris<Barycentric>::NewtonianMotionEquation>::Instance>
      fixed_instance_;

  PartTo<DegreesOfFreedom<RigidPileUp>> actual_part_degrees_of_freedom_;
  PartTo<DegreesOfFreedom<ApparentBubble>> apparent_part_degrees_of_freedom_;

  // Called in the destructor.
  std::function<void()> deletion_callback_;

  friend class TestablePileUp;
};

// A convenient data object to track a pile-up and the result of integrating it.
struct PileUpFuture {
  PileUpFuture(not_null<PileUp const*> pile_up, std::future<Status> future);
  not_null<PileUp const*> pile_up;
  std::future<Status> future;
};

}  // namespace internal_pile_up

using internal_pile_up::PileUp;
using internal_pile_up::PileUpFuture;

}  // namespace ksp_plugin
}  // namespace principia
