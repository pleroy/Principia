#include "base/heap_checker.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace principia {
namespace testing_utilities {

class HeapCheckedTest : public testing::Test {
protected:
  HeapCheckedTest() : testing::Test() {
    FLAGS_minloglevel = google::FATAL;
  }
  
 private:
  base::HeapChecker heap_checker_;
};

}  // namespace testing_utilities
}  // namespace principia
