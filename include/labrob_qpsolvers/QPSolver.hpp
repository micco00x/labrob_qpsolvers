#pragma once

namespace labrob {
namespace qpsolvers {

template <typename Scalar>
class QPSolver {
 public:
  virtual void solve(
      Scalar* H,
      Scalar* f,
      Scalar* A,
      Scalar* b,
      Scalar* C,
      Scalar* d_min,
      Scalar* d_max) = 0;
  virtual Scalar* get_solution() const = 0;

  const int num_variables_;
  const int num_equality_constraints_;
  const int num_inequality_constraints_;
}; // end class QPSolver

} // end namespace labrob::qpsolvers
} // end namespace labrob