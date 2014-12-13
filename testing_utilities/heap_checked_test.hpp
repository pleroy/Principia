#include "base/heap_checker.hpp"
#include "gtest/gtest.h"

namespace principia {
namespace testing_utilities {

class HeapCheckedTest : public testing::Test {
 private:
  base::HeapChecker heap_checker_;
};

}  // namespace testing_utilities
}  // namespace principia
