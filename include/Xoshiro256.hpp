#pragma once
#include <cstdint>
#include <array>
#include <ctime>
#include <limits>
#include <type_traits>

namespace XoshiroCpp {
inline static const std::uint64_t kDefaultSeed = std::time(nullptr);
class SplitMix64 {
 public:
  using StateType = std::uint64_t;
  using ResultType = std::uint64_t;
  [[nodiscard]] explicit constexpr SplitMix64(
      StateType state = kDefaultSeed) noexcept;
  constexpr ResultType operator()() noexcept;
  template <std::size_t N>
  [[nodiscard]] constexpr std::array<std::uint64_t, N>
  generateSeedSequence() noexcept;

 private:
  StateType state_;
};

class Xoshiro256PlusPlus {
 public:
  using StateType = std::array<std::uint64_t, 4>;
  using ResultType = std::uint64_t;
  [[nodiscard]] explicit constexpr Xoshiro256PlusPlus(
      std::uint64_t seed = kDefaultSeed) noexcept;
  [[nodiscard]] explicit constexpr Xoshiro256PlusPlus(StateType state) noexcept;
  constexpr ResultType operator()() noexcept;

 private:
  constexpr void jump() noexcept;
  constexpr void longJump() noexcept;
  StateType state_;
};

}  // namespace XoshiroCpp

namespace detail {
[[nodiscard]] static constexpr std::uint64_t RotL(const std::uint64_t x,
                                                  const int s) noexcept {
  return (x << s) | (x >> (64 - s));
}

[[nodiscard]] static constexpr std::uint32_t RotL(const std::uint32_t x,
                                                  const int s) noexcept {
  return (x << s) | (x >> (32 - s));
}
}  // namespace detail

inline constexpr XoshiroCpp::SplitMix64::SplitMix64(
    const StateType state) noexcept
    : state_(state) {}

inline constexpr XoshiroCpp::SplitMix64::ResultType
XoshiroCpp::SplitMix64::operator()() noexcept {
  std::uint64_t z = (state_ += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

template <std::size_t N>
inline constexpr std::array<std::uint64_t, N>
XoshiroCpp::SplitMix64::generateSeedSequence() noexcept {
  std::array<std::uint64_t, N> seeds = {};
  for (auto& seed : seeds) {
    seed = operator()();
  }
  return seeds;
}

inline constexpr XoshiroCpp::Xoshiro256PlusPlus::Xoshiro256PlusPlus(
    const std::uint64_t seed) noexcept
    : state_(SplitMix64{seed}.generateSeedSequence<4>()) {}

inline constexpr XoshiroCpp::Xoshiro256PlusPlus::Xoshiro256PlusPlus(
    const StateType state) noexcept
    : state_(state) {}

inline constexpr XoshiroCpp::Xoshiro256PlusPlus::ResultType
XoshiroCpp::Xoshiro256PlusPlus::operator()() noexcept {
  const std::uint64_t result =
      detail::RotL(state_[0] + state_[3], 23) + state_[0];
  const std::uint64_t t = state_[1] << 17;
  state_[2] ^= state_[0];
  state_[3] ^= state_[1];
  state_[1] ^= state_[2];
  state_[0] ^= state_[3];
  state_[2] ^= t;
  state_[3] = detail::RotL(state_[3], 45);
  return result;
}

inline constexpr void XoshiroCpp::Xoshiro256PlusPlus::jump() noexcept {
  constexpr std::uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c,
                                    0xa9582618e03fc9aa, 0x39abdc4529b1661c};

  std::uint64_t s0 = 0;
  std::uint64_t s1 = 0;
  std::uint64_t s2 = 0;
  std::uint64_t s3 = 0;

  for (std::uint64_t jump : JUMP) {
    for (int b = 0; b < 64; ++b) {
      if (jump & 1ULL << b) {
        s0 ^= state_[0];
        s1 ^= state_[1];
        s2 ^= state_[2];
        s3 ^= state_[3];
      }
      operator()();
    }
  }

  state_[0] = s0;
  state_[1] = s1;
  state_[2] = s2;
  state_[3] = s3;
}

inline constexpr void XoshiroCpp::Xoshiro256PlusPlus::longJump() noexcept {
  constexpr std::uint64_t LONG_JUMP[] = {0x76e15d3efefdcbbf, 0xc5004e441c522fb3,
                                         0x77710069854ee241,
                                         0x39109bb02acbe635};

  std::uint64_t s0 = 0;
  std::uint64_t s1 = 0;
  std::uint64_t s2 = 0;
  std::uint64_t s3 = 0;

  for (std::uint64_t jump : LONG_JUMP) {
    for (int b = 0; b < 64; ++b) {
      if (jump & 1ULL << b) {
        s0 ^= state_[0];
        s1 ^= state_[1];
        s2 ^= state_[2];
        s3 ^= state_[3];
      }
      operator()();
    }
  }

  state_[0] = s0;
  state_[1] = s1;
  state_[2] = s2;
  state_[3] = s3;
}
