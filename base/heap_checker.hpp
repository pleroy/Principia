#pragma once

#ifdef _DEBUG
#include <crtdbg.h>
#endif

namespace principia {
namespace base {

// This is a RAII class that is intented to be used in test to detect storage
// leaks.  Just declare an object of this type in the test fixture and it will
// automatically compare memory usage before and after the test and fail the
// test in case of discrepancy:
//
//  class MyFunkyTest : public testing::Test {
//   private:
//    base::HeapChecker heap_checker_;
//  };
//
// The output is not so nice as it doesn't report where the leak comes from.
// Sorry about that.
// This class is mostly harmless in Release builds.
class HeapChecker {
 public:
  HeapChecker();
  ~HeapChecker();

 private:
#ifdef _DEBUG
  _CrtMemState memory_at_construction_;
#endif
};

}  // namespace base
}  // namespace principia

#include "base/heap_checker_body.hpp"
