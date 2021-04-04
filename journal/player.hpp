﻿
#pragma once

#include <filesystem>
#include <fstream>
#include <map>
#include <memory>

#include "serialization/journal.pb.h"

#define PRINCIPIA_PLAYER_ALLOW_VERSION_MISMATCH 0

namespace principia {
namespace journal {

class Player final {
 public:
  using PointerMap = std::map<std::uint64_t, void*>;

  explicit Player(std::filesystem::path const& path);

  // Replays the next message in the journal.  Returns false at end of journal.
  // |index| is the 0-based index of the message in the journal.
  bool Play(int index);

  // Return the last replayed messages.
  serialization::Method const& last_method_in() const;
  serialization::Method const& last_method_out_return() const;

 private:
  // Reads one message from the stream.  Returns a |nullptr| at end of stream.
  std::unique_ptr<serialization::Method> Read();

  template<typename Profile>
  bool RunIfAppropriate(serialization::Method const& method_in,
                        serialization::Method const& method_out_return);

  PointerMap pointer_map_;
  std::ifstream stream_;

  std::unique_ptr<serialization::Method> last_method_in_;
  std::unique_ptr<serialization::Method> last_method_out_return_;

  friend class PlayerTest;
  friend class RecorderTest;
};

}  // namespace journal
}  // namespace principia

#include "journal/player_body.hpp"
