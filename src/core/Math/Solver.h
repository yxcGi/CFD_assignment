#ifndef SOLVER_H_
#define SOLVER_H_

#include "SparseMatrix.hpp"


class Solver
{
public:
    Solver() = delete;
    explicit Solver(const SparseMatrix<Scalar>& matrix);

private:
    const SparseMatrix<Scalar>& equation_;
};





#endif // SOLVER_H_