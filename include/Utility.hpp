#pragma once

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <tsl/robin_map.h>
#include <tsl/robin_set.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>

#include "Xoshiro256.hpp"

/* ---------------------------- */

#if defined(LOGGING) || defined(DEBUG)
#define TRACE_ARGLESS(text) fmt::print(text);
#define TRACE(text, ...) fmt::print(text, __VA_ARGS__);
#else
#define TRACE_ARGLESS(text)
#define TRACE(text, ...)
#endif

/* ---------------------------- */

namespace Utility {
static XoshiroCpp::Xoshiro256PlusPlus generator;
static constexpr auto EPSILON = 10e-6;
static Eigen::IOFormat FORMATTER(4, 0, ", ", "\n", "[", "]");
namespace Timer {
template <typename C = std::chrono::high_resolution_clock>
class TimerClass {
  const typename C::time_point start_point;
 public:
  TimerClass() : start_point(C::now()) {}
  template <typename R = typename C::duration::rep,
            typename U = typename C::duration>
  R elapsed() const {
    std::atomic_thread_fence(std::memory_order_relaxed);
    auto counted_time = std::chrono::duration_cast<U>(C::now() - start_point).count();
    std::atomic_thread_fence(std::memory_order_relaxed);
    return static_cast<R>(counted_time);
  }
};
using precise_stopwatch = TimerClass<>;
using system_stopwatch = TimerClass<std::chrono::system_clock>;
using monotonic_stopwatch = TimerClass<std::chrono::steady_clock>;
}  // namespace Timer

template <typename T>
inline T randomScaled(T min, T max) {
  thread_local static std::mt19937 eng{std::random_device{}()};
  thread_local static std::uniform_int_distribution<T> dist;
  return dist(eng, typename std::uniform_int_distribution<T>::param_type{min, max});
}

namespace FloatingPoint {
  inline bool IsInteger(const float x) {
    return std::abs(std::round(x) - x) < EPSILON;
  }

  inline bool AlmostZero(const float x) {
    return std::abs(x) < EPSILON;
  }

  inline bool CloseToIntegerByAbs(const float x, const std::int32_t integer) {
    return std::abs(std::round(x)) == integer;
  }

  inline float RoundToNearestIfInteger(const float x) {
    if (Utility::FloatingPoint::IsInteger(x)) {
      return std::lroundf(x);
    }
    return x;
  }

  template <std::ranges::range Iterable, typename = std::enable_if_t<std::is_same_v<typename Iterable::value_type, float>>>
  inline Iterable FixFloatingPointError(Iterable iterable) {
    for (float& f: iterable) {
      f = Utility::FloatingPoint::RoundToNearestIfInteger(f);
      if (Utility::FloatingPoint::AlmostZero(f)) {
        f = 0;
      }
    }
    return iterable;
  }

}

inline auto rand01Vec(const Eigen::Index size) {
  Eigen::VectorXf vector;
  vector.resize(size);
  for (Eigen::Index i = 0; i < size; ++i) {
    vector[i] = static_cast<float>(Utility::generator() % 2);
  }
  return vector;
}

inline Eigen::MatrixXf parseMatrix(const std::string& filename) {
  std::ifstream stream(filename);
  std::vector<std::vector<float>> stlMatrix;
  Eigen::MatrixXf matrix;
  std::string parsed;
  while (std::getline(stream, parsed)) {
    std::istringstream iss(parsed);
    std::int32_t token;
    std::vector<float> row;
    while (iss >> token) {
      row.emplace_back(token);
    }
    stlMatrix.emplace_back(std::move(row));
  }
  const auto size = static_cast<Eigen::Index>(stlMatrix.size());
  matrix.resize(size, size);
  for (std::uint32_t i = 0; i < stlMatrix.size(); ++i) {
    auto& row = stlMatrix[i];
    Eigen::VectorXf vector = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(row.data(), static_cast<Eigen::Index>(row.size()));
    matrix.row(i) = vector;
  }
  return matrix;
}

inline std::vector<Eigen::MatrixXf> parseMatrices(const std::string& filename) {
  static const auto toMatrix = [](std::vector<std::vector<float>>& stlMatrix) {
    Eigen::MatrixXf matrix;
    const auto size = static_cast<Eigen::Index>(stlMatrix.size());
    matrix.resize(size, size);
    for (std::uint32_t i = 0; i < stlMatrix.size(); ++i) {
      auto& row = stlMatrix[i];
      Eigen::VectorXf vector = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(row.data(), static_cast<Eigen::Index>(row.size()));
      matrix.row(i) = vector;
    }
    return matrix;
  };
  static const auto split = [](const std::string& string){
    static const auto splitter = [](const std::string& str) {
        std::vector<float> vector;
        std::stringstream ss(str);
        for (int i; ss >> i;) {
            vector.push_back(i);
            if (ss.peek() == ',' || ss.peek() == ' ') ss.ignore();
        }
        return vector;
    };
    const std::string_view digits = "0123456789";
    const auto first = std::find_first_of(string.cbegin(), string.cend(), digits.cbegin(), digits.cend());
    const auto last = std::find_first_of(string.crbegin(), string.crend(), digits.cbegin(), digits.cend());
    std::string toParse(first, last.base());
    auto parsed = splitter(toParse);
    return parsed;
  };
  std::ifstream stream(filename);
  std::vector<std::vector<float>> stlMatrix;
  std::vector<Eigen::MatrixXf> matrices;
  std::string parsed;
  while (std::getline(stream, parsed)) {
    if (parsed.empty() || parsed == "\r") {
      matrices.emplace_back(toMatrix(stlMatrix));
      stlMatrix.clear();
      continue;
    }
    stlMatrix.emplace_back(split(parsed));
  }
  return matrices;
}

}  // namespace Utility
