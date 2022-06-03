#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <sstream>
#include <filesystem>

#include "Solver.hpp"
#include "Utility.hpp"

/*
  auto t1 = Utility::Timer::system_stopwatch();
  auto matrix = Utility::parseMatrix("/home/user/matrix.txt");
  auto e1 = t1.elapsed();
  std::cout << e1 << std::endl;
  auto t2 = Utility::Timer::system_stopwatch();
  auto sol = Solver().generateUnimodularMatrices(matrix.transpose());
  auto e2 = t1.elapsed();
  std::cout << e2 << std::endl;
  std::cout << sol.size() << std::endl;
*/

std::vector<std::int32_t> split(const std::string& string) {
  static const auto splitter = [](const std::string& str) {
      std::vector<std::int32_t> vector;
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
}

int main() {
    const auto matrix = Utility::parseMatrix("/home/user/matrix.txt");
    Solver solver;
    Eigen::VectorXf vector;
    vector.resize(2);
    vector << 1, 2;
    solver.solveForFixedVector(matrix, );
}
