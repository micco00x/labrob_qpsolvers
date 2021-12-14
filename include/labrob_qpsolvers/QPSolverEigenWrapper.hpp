#pragma once

#include "QPSolver.hpp"

// STL
#include <memory>
#include <iostream>

// Eigen
#include <Eigen/Core>

namespace labrob {
namespace qpsolvers {

template <typename Scalar>
class QPSolverEigenWrapper {
 public:
  QPSolverEigenWrapper(std::shared_ptr<QPSolver<Scalar>> qp_solver_ptr) :
      qp_solver_ptr_(qp_solver_ptr) {

  }

  void
  solve(
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& H,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& g,
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& A,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& b,
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& C,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& lg,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& ug) {
    qp_solver_ptr_->solve(
        H.data(),
        g.data(),
        A.data(),
        b.data(),
        C.data(),
        lg.data(),
        ug.data()
    );
  }

  Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
  get_solution() {
    Scalar* solution_ptr = qp_solver_ptr_->get_solution();
    const int num_variables = qp_solver_ptr_->num_variables_;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> solution(num_variables);
    for (int idx = 0; idx < num_variables; ++idx) {
      solution(idx) = *(solution_ptr + idx);
    }
    return solution;
  }

 private:
  std::shared_ptr<QPSolver<Scalar>> qp_solver_ptr_;
}; // end class QPSolverEigenWrapper

} // end namespace qpsolvers
} // end namespace labrob