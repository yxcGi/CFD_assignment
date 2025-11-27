#ifndef SOLVER_H_
#define SOLVER_H_

#include "SparseMatrix.hpp"


class Solver
{
public:
    // 求解方法
    enum class Method
    {
        Jacobi,             // 雅可比迭代
        GaussSeidel,        // 高斯赛德尔迭代
        AMG,                // 代数多重网格法
    };

public:
    Solver() = delete;
    explicit Solver(const SparseMatrix<Scalar>& matrix);
    Solver(const Solver&) = delete;
    Solver(Solver&&) = delete;
    ~Solver() = default;

public:
    
private:
    const SparseMatrix<Scalar>& equation_;  // 方程矩阵（稀疏矩阵，压不压缩都行）
    std::vector<Scalar> X0_;        // 上一步的解
    std::vector<Scalar> X_;         // 本部的解
    Scalar relaxationFactor_ = 1.0; // 松弛因子
};





#endif // SOLVER_H_