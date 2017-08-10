#pragma once

#include "ksp_plugin/planetarium.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using geometry::RigidTransformation;
using quantities::si::Metre;
using quantities::si::Radian;
using testing_utilities::make_not_null;

class MockPlanetarium : public Planetarium {
 public:
  MockPlanetarium()
      : Planetarium(Planetarium::Parameters(1.0, 1.0 * Radian, 1.0 * Radian),
                    Perspective<Navigation, Camera>(
                        RigidTransformation<Navigation, Camera>::Identity(),
                        1 * Metre),
                    make_not_null<Ephemeris<Barycentric> const*>(),
                    make_not_null<NavigationFrame const*>()) {}
};

}  // namespace internal_planetarium

using internal_planetarium::MockPlanetarium;

}  // namespace ksp_plugin
}  // namespace principia
