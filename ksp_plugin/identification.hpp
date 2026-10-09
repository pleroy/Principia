#pragma once

#include <cstdint>
#include <ostream>
#include <set>
#include <string>

#include "absl/container/btree_map.h"
#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "base/not_null.hpp"

namespace principia {
namespace ksp_plugin {

FORWARD_DECLARE(class, Part, FROM(part), INTO(identification));
FORWARD_DECLARE(class, Vessel, FROM(vessel), INTO(identification));

namespace _identification {
namespace internal {

using namespace principia::base::_not_null;

class UUID {
 public:
  UUID(std::uint64_t bytes_0_7, std::uint64_t bytes_8_15);

  void WriteToMessage(not_null<serialization::UUID*> message) const;
  static UUID ReadFromMessage(serialization::UUID const& message);

 private:
  std::uint64_t bytes_0_7_;
  std::uint64_t bytes_8_15_;

  friend std::ostream& operator<<(std::ostream& out, UUID const& date);
  friend bool operator==(UUID const& lhs, UUID const& rhs);
  template<typename H>
  friend H AbslHashValue(H h, UUID const& m);
};

std::ostream& operator<<(std::ostream& out, UUID const& date) {}

bool operator==(UUID const& lhs, UUID const& rhs) {
  return lhs.bytes_0_7_ == rhs.bytes_0_7_ && lhs.bytes_8_15_ == rhs.bytes_8_15_;
}

template<typename H>
H AbslHashValue(H h, UUID const& uuid) {
  return H::combine(std::move(h), uuid.bytes_0_7_., uuid.bytes_8_15_);
}

// The GUID of a vessel, obtained by `v.id.ToString()` in C#. We use this as a
// key in a map.
using GUID = std::string;

// Corresponds to KSP's `Part.flightID`, *not* to `Part.uid`.  C#'s `uint`
// corresponds to `uint32_t`.
using PartId = std::uint32_t;

// Comparator by PartId.  Useful for ensuring a consistent ordering in sets of
// pointers to Parts.
struct PartByPartIdComparator {
  bool operator()(not_null<Part*> left, not_null<Part*> right) const;
  bool operator()(not_null<Part const*> left,
                  not_null<Part const*> right) const;
};

// Comparator by GUID.  Useful for ensuring a consistent ordering in sets of
// pointers to Vessels.
struct VesselByGUIDComparator {
  bool operator()(not_null<Vessel*> left, not_null<Vessel*> right) const;
  bool operator()(not_null<Vessel const*> left,
                  not_null<Vessel const*> right) const;
};

template<typename T>
using PartTo = absl::btree_map<not_null<Part*>,
                               T,
                               PartByPartIdComparator>;
using VesselSet = std::set<not_null<Vessel*>,
                           VesselByGUIDComparator>;
using VesselConstSet = std::set<not_null<Vessel const*>,
                                VesselByGUIDComparator>;

}  // namespace internal

using internal::GUID;
using internal::PartByPartIdComparator;
using internal::PartId;
using internal::PartTo;
using internal::VesselByGUIDComparator;
using internal::VesselConstSet;
using internal::VesselSet;

}  // namespace _identification
}  // namespace ksp_plugin
}  // namespace principia
