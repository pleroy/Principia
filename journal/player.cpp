#include "journal/player.hpp"

#include <chrono>
#include <filesystem>
#include <memory>
#include <string>

#include "absl/strings/match.h"
#include "base/array.hpp"
#include "base/get_line.hpp"
#include "base/hexadecimal.hpp"
#include "base/version.hpp"
#include "glog/logging.h"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.

#define PRINCIPIA_PLAYER_ALLOW_VERSION_MISMATCH 0

namespace principia {
namespace journal {
namespace _player {
namespace internal {

using interface::principia__ActivatePlayer;
using namespace principia::base::_array;
using namespace principia::base::_get_line;
using namespace principia::base::_hexadecimal;
using namespace principia::base::_version;

using namespace std::chrono_literals;

Player::Player(std::filesystem::path const& path)
    : stream_(path, std::ios::in) {
  principia__ActivatePlayer();
  CHECK(!stream_.fail());
}

bool Player::Play(int const index) {
  return Process(/*method_in=*/Read(), index, /*play=*/true);
}

bool Player::Scan(int const index) {
  return Process(/*method_in=*/Read(), index, /*play=*/false);
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

  static auto* const encoder = new HexadecimalEncoder</*null_terminated=*/true>;
  auto const bytes = encoder->Decode({line.c_str(), strlen(line.c_str())});
  auto method = std::make_unique<serialization::Method>();
  CHECK(method->ParseFromArray(bytes.data.get(), static_cast<int>(bytes.size)));

  return method;
}

bool Player::Process(std::unique_ptr<serialization::Method> method_in,
                     int const index, bool const play) {
  if (method_in == nullptr) {
    // End of input file.
    return false;
  }
  std::unique_ptr<serialization::Method> method_out_return = Read();
  if (method_out_return == nullptr) {
    LOG(ERROR) << "Unpaired method:\n" << method_in->DebugString();
    return false;
  }

  // Check that the version of the journal matches that of the binary.  Remember
  // that a GetVersion message is logged in Recorder::Activate, so it's always
  // present.  The |StartsWith| test below allows the "-dirty" suffix, which is
  // typically present when debugging.
  if (method_in->HasExtension(serialization::GetVersion::extension)) {
    auto const& get_version_out =
        method_out_return->GetExtension(serialization::GetVersion::extension)
            .out();
    LOG_IF(FATAL,
           !absl::StartsWith(Version, get_version_out.version())  &&
               (PRINCIPIA_PLAYER_ALLOW_VERSION_MISMATCH == 0))
        << "Journal version is " << get_version_out.version()
        << ", running with a binary built at version " << Version
        << "; set PRINCIPIA_PLAYER_ALLOW_VERSION_MISMATCH to 1 in player.cpp "
        << "if this is intended.";
  }

#if 0
  LOG_IF(
      ERROR,
      method_in->HasExtension(serialization::VesselGetAnalysis::extension) ||
      method_in->HasExtension(serialization::DeleteInterchange::extension))
      << "index: " << index << "\n"
      << method_in->ShortDebugString() << "\n"
      << method_out_return->ShortDebugString();
#endif
#if 1
  LOG_IF(ERROR, index > 58889025 - 50) << "index: " << index << "\n"
                                 << method_in->ShortDebugString() << "\n"
                                 << method_out_return->ShortDebugString();
#endif

  if (play) {
    auto const before = std::chrono::system_clock::now();

#include "journal/player.generated.cc"

    auto const after = std::chrono::system_clock::now();
    if (after - before > 100ms) {
      LOG(ERROR) << "Long method (" << (after - before) / 1ms << " ms):\n"
                 << method_in->DebugString();
    }
  }

  last_method_in_.swap(method_in);
  last_method_out_return_.swap(method_out_return);

  return true;
}

}  // namespace internal
}  // namespace _player
}  // namespace journal
}  // namespace principia
