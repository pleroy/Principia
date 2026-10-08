#include "ksp_plugin/part.hpp"
#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {
namespace _identification {
namespace internal {

using namespace principia::ksp_plugin::_part;
using namespace principia::ksp_plugin::_vessel;

Index::Index(int const index)
    : value_(index) {}

int Index::value() const {
  return value_;
}

bool operator==(Index const& lhs, Index const& rhs) {
  return lhs.value_ == rhs.value_;
}

std::ostream& operator<<(std::ostream& out, Index const& index) {
  // TODO: insert return statement here
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
