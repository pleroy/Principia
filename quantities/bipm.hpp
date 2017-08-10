﻿
#pragma once

#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

// This namespace contains the other non-SI units listed in the BIPM's
// SI brochure 8, section 4.1, table 8,
// http://www.bipm.org/en/si/si_brochure/chapter4/table8.html.
namespace bipm {

constexpr Pressure Bar                 = 1e5 * si::Pascal;
constexpr Pressure MillimetreOfMercury = 133.322 * si::Pascal;
constexpr Length   Ångström            = 1e-10 * si::Metre;
constexpr Length   NauticalMile        = 1852 * si::Metre;
constexpr Speed    Knot                = 1 * NauticalMile / si::Hour;
constexpr Area     Barn                = 1e-28 * Pow<2>(si::Metre);

}  // namespace bipm
}  // namespace quantities
}  // namespace principia
