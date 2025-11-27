#ifndef SOLVER_H_
#define SOLVER_H_

#include "SparseMatrix.hpp"
#include "Field.hpp"


template <typename Tp>
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
    explicit Solver(SparseMatrix<Tp>& matrix, Method method, int maxmaxIterationNum);
    Solver(const Solver&) = delete;
    Solver(Solver&&) = delete;
    ~Solver() = default;

public:
    // 初始化
    void init(const std::vector<Tp>& x0);
    void init(Tp value0);
    void init(const Field<Tp>& initField);

    // 求解方程


private:
    SparseMatrix<Tp>& equation_;  // 方程矩阵（稀疏矩阵，压不压缩都行）
    std::vector<Tp> x0_;        // 上一步的解
    std::vector<Tp> x_;         // 本步的解
    Scalar relaxationFactor_{ 1.0 };       // 松弛因子
    Method method_;             // 求解方法
    int maxIterationNum_;       // 最大迭代次数
    bool isInitialized_{ false };           // 是否赋初始值
};






template<typename Tp>
inline Solver<Tp>::Solver(SparseMatrix<Tp>& matrix, Method method, int maxmaxIterationNum)
    : equation_(matrix)
    , method_(method)
    , maxIterationNum_(maxmaxIterationNum)
{
    // 检查矩阵是否有效
    if (!equation_.isValid())
    {
        std::cerr << "Solver<Tp>::Solver(const SparseMatrix<Tp>& matrix) Error: matrix is not valid" << std::endl;
        throw std::invalid_argument("matrix is not valid");
    }

    // 矩阵有效，判断是否压缩，没压缩，则进行压缩
    if (!equation_.isCompressed())
    {
        equation_.compress();
    }
}

template<typename Tp>
inline void Solver<Tp>::init(const std::vector<Tp>& x0)
{
    // 不能重复初始化
    if (isInitialized_)
    {
        std::cerr << "Solver<Tp>::init(const std::vector<Tp>& x0) Error: cannot initialize again" << std::endl;
        throw std::runtime_error("cannot initialize again");
    }

    // 检查长度是否合法
    if (x0.size() != equation_.size())
    {
        std::cerr << "Solver<Tp>::init(const std::vector<Tp>& x0) Error: length of x0 is not equal to matrix size, length of x0i is illegal" << std::endl;
        throw std::invalid_argument("length of x0 is not equal to matrix size, length of x0i is illegal");
    }
    // 赋初值
    x0_ = x0;
    isInitialized_ = true;
}

template<typename Tp>
inline void Solver<Tp>::init(Tp value0)
{
    // 不能重复初始化
    if (isInitialized_)
    {
        std::cerr << "Solver<Tp>::init(Tp value0) Error: cannot initialize again" << std::endl;
        throw std::runtime_error("cannot initialize again");
    }

    x0_.resize(equation_.size(), value0);
    isInitialized_ = true;
}

template<typename Tp>
inline void Solver<Tp>::init(const Field<Tp>& initField)
{
    if (isInitialized_)
    {
        std::cerr << "Solver<Tp>::init(const Field<Tp>& initField) Error: cannot initialize again" << std::endl;
        throw std::runtime_error("cannot initialize again");
    }

    // 判断场是否有效
    if (!initField.isValid())
    {
        std::cerr << "Solver<Tp>::init(const Field<Tp>& initField) Error: initField is not valid" << std::endl;
        throw std::invalid_argument("initField is not valid");
    }

    Mesh* equationMesh = equation_.getMesh();
    Mesh* initFieldMesh = initField.getMesh();

    // 判断是否为同一个场
    // 如果稀疏矩阵不是由网格构造，则可忽略
    if (equationMesh != nullptr &&
        initFieldMesh != equationMesh)
    {
        std::cerr << "The mesh of the field is not equal to the mesh of the matrix" << std::endl;
        throw std::invalid_argument("The mesh of the field is not equal to the mesh of the matrix");
    }

    // 判断长度是否合理
    if (equation_.size() != initFieldMesh->getCellNumber())
    {
        std::cerr << "Solver<Tp>::init(const Field<Tp>& initField) Error: The size of equation matrix is not equal to the number of cells in the mesh" << std::endl;
        throw std::invalid_argument("The size of equation matrix is not equal to the number of cells in the mesh");
    }

    // 从field中的cellField_0_获取初始值
    x0_ = initField.getCellField_0().getData();
    isInitialized_ = true;
}






#endif // SOLVER_H_