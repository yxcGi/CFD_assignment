#include "Solver.h"

Solver::Solver(const SparseMatrix<Scalar>& matrix)
    : equation_(matrix)
{
    // 检查matrix 是否有效
    if (!equation_.isValid())
    {
        std::cerr << "Solver::Solver(const SparseMatrix<Scalar>& matrix) Error: This matrix is not valid." << std::endl;
        throw std::invalid_argument("this matrix is not valid");
    }
    
}