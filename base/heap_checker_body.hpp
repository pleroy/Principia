#pragma once

#include "base/heap_checker.hpp"

#ifdef _DEBUG
#include "gtest/gtest.h"
#endif

// NOTE(phl): Strictly speaking this doesn't need to be a _body.hpp file, it
// might as well be a .cpp, but that's just more complicated in our setup.

namespace principia {
namespace base {

inline HeapChecker::HeapChecker() {
#ifdef _DEBUG
  _CrtMemCheckpoint(&memory_at_construction_);
  _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
  _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT );
  _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);
#endif
}

inline HeapChecker::~HeapChecker() {
#ifdef _DEBUG
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
#endif
}

}  // namespace base
}  // namespace principia
