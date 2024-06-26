#pragma once

#include <cstdint>
#include <memory>

#include "base/thread_pool.hpp"

namespace principia {
namespace base {
namespace _graveyard {
namespace internal {

using namespace principia::base::_thread_pool;

class Graveyard {
 public:
  explicit Graveyard(std::int64_t number_of_threads);

  template<typename T>
  void Bury(std::unique_ptr<T> t);

 private:
  ThreadPool<void> gravedigger_;
};

}  // namespace internal

using internal::Graveyard;

}  // namespace _graveyard
}  // namespace base
}  // namespace principia

#include "base/graveyard_body.hpp"
