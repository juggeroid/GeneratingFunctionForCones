#include <fmt/ostream.h>

#include "Solver.hpp"

namespace {
constexpr auto kSolutionAttempts = 10ULL;
constexpr auto kAvgDetResizeFactor = 10ULL;

}  // namespace

Eigen::MatrixXf Solver::getMatrix() const noexcept { return matrix_; }
const Eigen::MatrixXf &Solver::getMatrixRef() const noexcept { return matrix_; }

void Solver::setMatrix(const Eigen::MatrixXf &matrix) noexcept { matrix_ = matrix; };

void Solver::setMatrix(Eigen::MatrixXf &&matrix) noexcept { matrix_ = std::move(matrix); };

void Solver::produceGeneratingFunction(CONE_TYPE type) noexcept {
  // TO DO: proper interface and filewriter
  std::ofstream fstream("unimodular.txt");
  unimodular_ = generateUnimodularMatrices(matrix_);
  fstream << "For matrix Ax <= 0, the short rational generating function is defined as: ";
  fstream << "\\sum_{i} e_i * 1 / (1 - x^{U_i1}) * ... * (1 - x^{U_in}).\n";
  fstream << "Unimodular (U):\n";
  if (type == CONE_TYPE::CONSTRAINTED) {
    for (auto& matrix: unimodular_) {
      auto& ref = matrix.second;
      ref = ref.inverse().transpose() * -1;
    }
  }
  for (std::size_t i = 0; i < unimodular_.size(); ++i) {
    fstream << fmt::format("U_{} = \n{}\n", i, unimodular_[i].second);
  }
  for (std::size_t i = 0; i < unimodular_.size(); ++i) {
    fstream << fmt::format("e_{} = {}\n", i, unimodular_[i].first);
  }
};

void Solver::produceGeneratingFunction(CONE_TYPE type, const Eigen::VectorXf &b) noexcept {
  // Re-factor later.
  std::ofstream fstream("unimodular_with_b.txt");
  Eigen::VectorXf numeratorCoeff;
  switch (type) {
    case CONE_TYPE::BASIC:
      numeratorCoeff = getCoefficientsForShiftedCone(matrix_, b);
      break;
    case CONE_TYPE::CONSTRAINTED:
      numeratorCoeff = solveForFixedVector(matrix_, b);
      break;
    default:
      break;
  }
  unimodular_ = generateUnimodularMatrices(matrix_);
  fstream << "For matrix Ax <= b, the short rational generating function is defined as: ";
  fstream << "\\sum_{i} e_i * x^{w_i = U*ceil(x)} / (1 - x^{U_i1}) * ... * (1 - x^{U_in}).\n";
  fstream << "Unimodular (U):\n";
  for (std::size_t i = 0; i < unimodular_.size(); ++i) {
    fstream << fmt::format("U_{} = \n{}\n", i, unimodular_[i].second);
  }
  for (std::size_t i = 0; i < unimodular_.size(); ++i) {
    fstream << fmt::format("e_{} = {}\n", i, unimodular_[i].first);
  }
  for (std::size_t i = 0; i < numeratorCoeff.size(); ++i) {
    fstream << fmt::format("w_{} = {}\n", i, numeratorCoeff[i]);
  }
};

tsl::robin_map<Eigen::Index, Solver::FloatType> Solver::findIndicesOfRational(
    const Eigen::VectorXf &vector) const noexcept {
  tsl::robin_map<Eigen::Index, FloatType> floatToIndex;
  TRACE_ARGLESS("\033[1mFloats:\033[0m \n");
  for (Eigen::Index i = 0; i < vector.size(); ++i) {
    if (!Utility::FloatingPoint::IsInteger(vector[i])) {
      TRACE("[{} (index), {} (value)]\n", i, vector[i]);
      floatToIndex.emplace(i, vector[i]);
    }
  }
  return floatToIndex;
};

tl::expected<Eigen::VectorXf, std::string> Solver::solveForRandZ2AndAdjust(
    const Eigen::MatrixXf &matrix) const noexcept {
  // Check the dimension validity.
  if ((matrix.rows() != matrix.cols()) || matrix.rows() == 0) {
    return tl::unexpected<std::string>("Matrix is non-square or empty, abort.");
  }
  // The probability of getting a complete integer solution is 1 / (2 ^ kSolutionAttempts).
  for (std::size_t i = 0ULL; i < kSolutionAttempts; ++i) {
    TRACE("\033[1mGenerating... (#{})\033[0m\n", i);
    // Generate random v \in Z^2.
    const Eigen::VectorXf randomVector = Utility::rand01Vec(matrix.rows());
    TRACE("\033[1mRandom [0/1]:\033[0m [{}]\n",
          randomVector.transpose().format(Utility::FORMATTER));
    // Solve the system of equations using Eigen (LLT/LLDT works faster, but the solutions are
    // often inaccurate).
    Eigen::VectorXf solution = matrix.colPivHouseholderQr().solve(randomVector);
    // If we've got floats that are close to 0, adjust them beforehand.
    std::ranges::for_each(solution, [](float &x) {
      if (Utility::FloatingPoint::AlmostZero(x)) {
        x = 0;
      }
    });
    TRACE("\033[1mFixed FP Solution:\033[0m [{}]\n",
          solution.transpose().format(Utility::FORMATTER));
    // Find the iterator to first rational component, rounding to integer if sol < EPS = 10e-7.
    const Eigen::VectorXf::const_iterator iterator = std::ranges::find_if(
        solution, [](FloatType x) { return !Utility::FloatingPoint::IsInteger(x); });
    // Try once more if not found.
    if (iterator == std::cend(solution)) {
      continue;
    }
    return tl::expected<Eigen::VectorXf, std::string>(adjustComponent(std::move(solution)));
  }
  // Something went wrong, abort.
  return tl::unexpected<std::string>("No rational solution has been found.");
}

Eigen::VectorXf Solver::adjustComponent(Eigen::VectorXf &&vector) const noexcept {
  TRACE("\033[1mSolution before adjustment:\033[0m [{}]\n",
        vector.transpose().format(Utility::FORMATTER));
  // Take the float reminder and save the sign if x = -0.5.
  // e.g. -5.8 -> 0.2, -5.4 -> -0.4 ...
  static constexpr auto adjust = [](float x) {
    const bool neg = std::signbit(x);
    x -= std::floor(x + 0.5);
    return !neg && x == -0.5 ? 0.5 : x;
  };
  // Adjust the floats and set the integer values as 0.
  auto adjusted = std::move(vector);
  for (Eigen::Index i = 0; i < adjusted.size(); ++i) {
    auto &component = adjusted[i];
    if (Utility::FloatingPoint::IsInteger(component)) {
      component = 0;
      continue;
    }
    component = adjust(component);
  }
  TRACE("\033[1mSolution after adjustment:\033[0m [{}]\n",
        adjusted.transpose().format(Utility::FORMATTER));
  return adjusted;
}

std::vector<std::pair<std::int32_t, Eigen::MatrixXf>> Solver::generateUnimodularMatrices(
    const Eigen::MatrixXf &matrix) const noexcept {
  // Allocate two queues for generated unimodular matrices and intermediate matrices.
  std::vector<std::pair<std::int32_t, Eigen::MatrixXf>> unimodular;
  std::vector<std::pair<std::int32_t, Eigen::MatrixXf>> intermediate;
  // There are no more than dim^log(det(M)) unimodular matrices.
  unimodular.reserve(
      static_cast<std::size_t>(std::pow(matrix.cols(), std::log(kAvgDetResizeFactor)) + 1ULL));
  // Push the initial matrix onto the queue.
  intermediate.emplace_back(1, matrix);

  while (!intermediate.empty()) {
    std::pair<std::int32_t, Eigen::MatrixXf> current = std::move(intermediate.back());
    intermediate.pop_back();
    const FloatType determinant = current.second.determinant();
    TRACE("Matrix: \n{}, det(Matrix): {}\n", current.second, determinant);
    // If det == 0, discard the branch.
    if (Utility::FloatingPoint::AlmostZero(determinant)) {
      continue;
    }
    // If det == 1, the unimodular matrix is found. Remove from the queue and continue.
    if (Utility::FloatingPoint::CloseToIntegerByAbs(determinant, 1)) {
      TRACE_ARGLESS("Unimodular, remove from the queue.\n");
      unimodular.emplace_back(std::move(current));
      continue;
    }
    // Attempt to solve for random v \in Z^2.
    auto solution = solveForRandZ2AndAdjust(current.second);
    if (!solution) {
      TRACE("{}", solution.error());
      continue;
    }
    auto adjustedSolution = solution.value();
    TRACE("\033[1mSolution:\033[0m [{}]\n",
          adjustedSolution.transpose().format(Utility::FORMATTER));
    // Fix the possible float value to prevent amassing errors
    adjustedSolution = Utility::FloatingPoint::FixFloatingPointError(adjustedSolution);
    // Find the indices of rational values to branch.
    auto indices = findIndicesOfRational(adjustedSolution);
    // Find the product A x T.
    Eigen::VectorXf product = Utility::FloatingPoint::FixFloatingPointError(current.second * adjustedSolution);
    TRACE("\033[1mA x T:\033[0m [{}]\n", product.transpose().format(Utility::FORMATTER));
    // Assign a new matrix, replace its i-th row by A x T, do so for every float found in the
    // adjusted solution. In total, there should be no more than dim(M) additional branches.
    for (const auto &[index, value] : indices) {
      std::pair<std::int32_t, Eigen::MatrixXf> modified = current;
      modified.second.col(index) = product.transpose();
      const float updatedDet = value * determinant;
      modified.first *= (updatedDet > 0) - (updatedDet < 0);
      TRACE("\033[1m+ Decomposition Queue:\033[0m \n[{}]\n", modified.second);
      intermediate.emplace_back(std::move(modified));
    }
  }
  return unimodular;
}

Eigen::VectorXf Solver::solveForFixedVector(
    const Eigen::MatrixXf &matrix, const Eigen::VectorXf &vector) noexcept {
  // Find the solution of Av = b system.
  Eigen::VectorXf solution = Utility::FloatingPoint::FixFloatingPointError(matrix.colPivHouseholderQr().solve(vector));
  // Find the transpose A = A^T.
  const Eigen::MatrixXf transposed = matrix.transpose();
  // Pass the vector v along with the transposed matrix to v + cone(A) algorithm.
  auto coefficients = getCoefficientsForShiftedCone(transposed, solution);
  return coefficients;
}

// v + cone(A).
Eigen::VectorXf
Solver::getCoefficientsForShiftedCone(const Eigen::MatrixXf& matrix, const Eigen::VectorXf& vector) noexcept {
  // Solve for the Ax = v; then set w = A * ceil(x).
  Eigen::VectorXf solution = Utility::FloatingPoint::FixFloatingPointError(matrix.colPivHouseholderQr().solve(vector));
  std::for_each(solution.begin(), solution.end(), [&](float& component) {
    component = std::ceil(component);
  });
  Eigen::VectorXf w = matrix * solution;
  // As a result, the generating function is x^w f(cone(A); x).
  return w;
}

