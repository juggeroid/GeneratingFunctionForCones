#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <tl/expected.hpp>

#include "Utility.hpp"

// No time to refactor, maybe later.

class Solver {
 public:
    enum CONE_TYPE {
      BASIC, // cone(A)
      CONSTRAINTED // Ax <= b
    };
 public:
  const Eigen::MatrixXf& getMatrixRef() const noexcept;
  Eigen::MatrixXf getMatrix() const noexcept;

  void setMatrix(const Eigen::MatrixXf& matrix) noexcept;
  void setMatrix(Eigen::MatrixXf&& matrix) noexcept;

  void produceGeneratingFunction(CONE_TYPE type) noexcept;
  void produceGeneratingFunction(CONE_TYPE type, const Eigen::VectorXf& b) noexcept;

 public:
  Eigen::MatrixXf matrix_;
  std::vector<std::pair<std::int32_t, Eigen::MatrixXf>> unimodular_;
  std::vector<Eigen::VectorXf> shifts_;

  using FloatType = float;

 public:
  tsl::robin_map<Eigen::Index, FloatType>
  findIndicesOfRational(const Eigen::VectorXf& vector) const noexcept;

  bool
  isRational(const Eigen::VectorXf& vector) const noexcept;

  [[nodiscard]]
  tl::expected<Eigen::VectorXf, std::string>
  solveForRandZ2AndAdjust(const Eigen::MatrixXf& matrix) const noexcept;

  [[nodiscard]]
  Eigen::VectorXf
  adjustComponent(Eigen::VectorXf&& vector) const noexcept;

  // Generate unimodular for an arbitrary cone(A) or Ax <= 0.
  [[nodiscard]]
  std::vector<std::pair<std::int32_t, Eigen::MatrixXf>>
  generateUnimodularMatrices(const Eigen::MatrixXf& matrix) const noexcept;

  // Solve for Ax <= b.
  [[nodiscard]]
  Eigen::VectorXf
  solveForFixedVector(const Eigen::MatrixXf& matrix, const Eigen::VectorXf& vector) noexcept;

  // Solve for v + cone(A).
  [[nodiscard]]
  Eigen::VectorXf
  getCoefficientsForShiftedCone(const Eigen::MatrixXf& matrix, const Eigen::VectorXf& vector) noexcept;

};
