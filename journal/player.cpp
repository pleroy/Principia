﻿
#include "journal/player.hpp"

#include <chrono>
#include <string>

#include "base/array.hpp"
#include "base/get_line.hpp"
#include "base/hexadecimal.hpp"
#include "journal/profiles.hpp"
#include "glog/logging.h"

namespace principia {

using base::GetLine;
using base::HexadecimalDecode;
using base::UniqueBytes;

namespace journal {

Player::Player(std::experimental::filesystem::path const& path)
    : stream_(path, std::ios::in) {
  CHECK(!stream_.fail());
}

bool Player::Play() {
  std::unique_ptr<serialization::Method> method_in = Read();
  if (method_in == nullptr) {
    // End of input file.
    return false;
  }
  std::unique_ptr<serialization::Method> method_out_return = Read();
  if (method_out_return == nullptr) {
    LOG(ERROR) << "Unpaired method:\n" << method_in->DebugString();
    return false;
  }

#if 0
  LOG(ERROR) << "\n" << method_in->ShortDebugString() << "\n"
             << method_out_return->ShortDebugString();
#endif

  auto const before = std::chrono::system_clock::now();

#include "journal/player.generated.cc"

  auto const after = std::chrono::system_clock::now();
  if (after - before > std::chrono::milliseconds(100)) {
    LOG(ERROR) << "Long method:\n" << method_in->DebugString();
  }

  last_method_in_.swap(method_in);
  last_method_out_return_.swap(method_out_return);

  return true;
}

serialization::Method const& Player::last_method_in() const {
  return *last_method_in_;
}

serialization::Method const& Player::last_method_out_return() const {
  return *last_method_out_return_;
}

std::unique_ptr<serialization::Method> Player::Read() {
  std::string const line = GetLine(stream_);
  if (line.empty()) {
    return nullptr;
  }

  std::uint8_t const* const hexadecimal =
      reinterpret_cast<std::uint8_t const*>(line.c_str());
  int const hexadecimal_size = strlen(line.c_str());
  UniqueBytes bytes(hexadecimal_size >> 1);
  HexadecimalDecode({hexadecimal, hexadecimal_size},
                    {bytes.data.get(), bytes.size});
  auto method = std::make_unique<serialization::Method>();
  CHECK(method->ParseFromArray(bytes.data.get(),
                               static_cast<int>(bytes.size)));

  return method;
}

}  // namespace journal
}  // namespace principia
