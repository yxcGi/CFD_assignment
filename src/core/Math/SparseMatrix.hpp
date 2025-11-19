#pragma once

#include <vector>
#include "Mesh.h"

/**
 * @brief 采用CSR方法存储稀疏矩阵
 * @tparam Tp 
 */
template <typename Tp>
class SparseMatrix
{
    using ULL = unsigned long long;
public:
    SparseMatrix() = default;
    SparseMatrix(ULL size);
    

private:
    std::vector<Tp> values_;        // 存储数据（行优先），大小为矩阵非0元素个数
    std::vector<ULL> colIndex_;     // 每个元素列索引, 与values一一对应
    std::vector<ULL> rowPointer_;   // 行指针，每一行起始元素的索引，大小为矩阵行数
};
