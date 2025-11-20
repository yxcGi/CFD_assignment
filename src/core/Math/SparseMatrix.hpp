#pragma once

#include <vector>
#include "Mesh.h"

/**
 * @brief 采用CSR方法存储稀疏矩阵，仅为方阵
 * @tparam Tp 一般为double，int 不会是Vector
 */
template <typename Tp>
class SparseMatrix
{
    using ULL = unsigned long long;
    constexpr static double EPSILON = 1e-10;
public:
    SparseMatrix() = default;
    SparseMatrix(ULL size);
    SparseMatrix(const std::vector<std::vector<Tp>>& matrix);   // 直接用二维数组初始化

public:
    // 打印
    void printMatrix() const;
    // 设置矩阵大小
    void setSize(ULL size);
    // 设置系数矩阵i, j元素值
    void setElement(ULL i, ULL j, Tp value);

private:
    std::vector<Tp> values_;        // 存储数据（行优先），大小为矩阵非0元素个数
    std::vector<ULL> colIndex_;     // 每个元素列索引, 与values一一对应
    std::vector<ULL> rowPointer_;   // 行指针，每一行起始元素的索引，大小为矩阵行数
    ULL size_;                      // 矩阵大小
    bool isValid_{ false };         // 是否有效
};

template<typename Tp>
inline SparseMatrix<Tp>::SparseMatrix(ULL size)
    : size_(size)
{}


template<typename Tp>
inline SparseMatrix<Tp>::SparseMatrix(const std::vector<std::vector<Tp>>& matrix)
    : size_(matrix.size())
{
    // 判断是否为方阵
    for (const std::vector<Tp>& row : matrix)
    {
        if (row.size() != size_)
        {
            std::cerr << "SparseMatrix<Tp>::SparseMatrix(const std::vector<std::vector<Tp>>& matrix) Error: matrix is not square" << std::endl;
            throw std::invalid_argument("matrix is not square");
        }
    }

    // 开始构造
    rowPointer_.push_back(0);
    for (ULL i = 0; i < size_; ++i)
    {
        for (ULL j = 0; j < size_; ++j)
        {
            const Tp ele = matrix[i][j];
            if (std::abs(ele) > EPSILON)   // 非0元素
            {
                values_.push_back(ele);
                colIndex_.push_back(j);
            }
        }
        // 每行结束后记录下一行的起始值
        rowPointer_.push_back(colIndex_.size());
    }
    isValid_ = true;
}

template<typename Tp>
inline void SparseMatrix<Tp>::printMatrix() const
{
    for (ULL i = 0; i < size_; ++i)
    {
        // 打印本行第一个非0元素前的0，保留四位小数
        ULL j = 0;
        std::cout << std::setprecision(4) << std::fixed;
        for (; j < colIndex_[rowPointer_[i]]; j++)
        {
            std::cout << std::setw(10) << 0.0;
        }

        // 开始依次判断并打印本行非0元素
        for (ULL index = rowPointer_[i]; index < rowPointer_[i + 1]; )
        {
            if (j == colIndex_[index])  // 该列为非0元素
            {
                std::cout << std::setw(10) << values_[index];
                ++index;
            }
            else
            {
                std::cout << std::setw(10) << 0.0;
            }
            ++j;
        }

        // 填充本行后面的0
        for (; j < size_; ++j)
        {
            std::cout << std::setw(10) << 0.0;
        }
        std::cout << "\n";
    }
}

template<typename Tp>
inline void SparseMatrix<Tp>::setSize(ULL size)
{
    if (isValid_)   // 如果已经有效，则不能再设置
    {
        std::cerr << "SparseMatrix<Tp>::setSize(ULL size) Error: matrix is valid" << std::endl;
        throw std::invalid_argument("matrix is valid");
    }
    size_ = size;
    isValid_ = true;
}
