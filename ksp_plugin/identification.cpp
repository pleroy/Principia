#include "ksp_plugin/part.hpp"
#include "ksp_plugin/vessel.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace _identification {
namespace internal {

using namespace principia::ksp_plugin::_part;
using namespace principia::ksp_plugin::_vessel;

UUID::UUID(std::uint64_t const bytes_0_7, std::uint64_t const bytes_8_15)
    : bytes_0_7_(bytes_0_7), bytes_8_15_(bytes_8_15) {}

void WriteToMessage(not_null<serialization::UUID*> const message) const {
  message->set_bytes_0_7(bytes_0_7_);
  message->set_bytes_8_15(bytes_8_15);
}

UUID UUID::ReadFromMessage(serialization::UUID const& message) {
  return UUID(message.bytes_0_7(), message.bytes_8_15());
}

bool operator==(UUID const& lhs, UUID const& rhs) {
  return lhs.bytes_0_7_ == rhs.bytes_0_7_ && lhs.bytes_8_15_ == rhs.bytes_8_15_;
}

std::ostream& operator<<(std::ostream& out, UUID const& uuid) {
  // TODO(phl): Output a proper UUID as it would be stringified by C#.
  out << uuid.bytes_0_7_ << "/" << uuid.bytes_8_15_;
  return out;
}

bool PartByPartIdComparator::operator()(not_null<Part*> const left,
                                        not_null<Part*> const right) const {
  return left->part_id() < right->part_id();
}

bool PartByPartIdComparator::operator()(
    not_null<Part const*> const left,
    not_null<Part const*> const right) const {
  return left->part_id() < right->part_id();
}

bool VesselByGUIDComparator::operator()(not_null<Vessel*> const left,
                                        not_null<Vessel*> const right) const {
  return left->guid() < right->guid();
}

bool VesselByGUIDComparator::operator()(
    not_null<Vessel const*> const left,
    not_null<Vessel const*> const right) const {
  return left->guid() < right->guid();
}

}  // namespace internal
}  // namespace _identification
}  // namespace ksp_plugin
}  // namespace principia
