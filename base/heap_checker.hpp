#pragma once

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#include "gtest/gtest.h"

namespace principia {
namespace base {

class HeapChecker {
 public:
  HeapChecker() {
    _CrtMemCheckpoint(&memory_at_construction_);
    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT );
    _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);
  }

  ~HeapChecker() {
    _CrtMemState memory_at_destruction;
    _CrtMemCheckpoint(&memory_at_destruction);
    _CrtMemState leaks;
    bool const has_leaks = _CrtMemDifference(&leaks,
                                             &memory_at_construction_,
                                             &memory_at_destruction) != 0;
    if (has_leaks) {
      _CrtMemDumpAllObjectsSince(&memory_at_construction_);
      EXPECT_FALSE(has_leaks) << "Storage leak detected:";
      _CrtMemDumpStatistics(&leaks);
    }
  }

  _CrtMemState memory_at_construction_;
};

}  // namespace base
}  // namespace principia
