#pragma once

#include <vector>
#include "Mesh.h"
#include <queue>


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
    SparseMatrix(const Mesh* mesh);     // 通过读取网格信息初始化

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
    rowPointer_.reserve(size_ + 1);

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
inline SparseMatrix<Tp>::SparseMatrix(const Mesh* mesh)
{
    if (!mesh->isValid())
    {
        std::cerr << "SparseMatrix<Tp>::SparseMatrix(const Mesh* mesh) Error: mesh is not valid" << std::endl;
        throw std::invalid_argument("mesh is not valid");
    }


    size_ = mesh->getCellNumber();

    // 先取出必要的变量
    const std::vector<Face>& faces = mesh->getFaces();

    // 给values_和colIndex_预分配内存,二者一一对应
    ULL elementNumber = 0;
    for (const Cell& cell : mesh->getCells())
    {
        elementNumber += cell.getFaceNum() + 1; // 若为内部纯内部单元则elementNumber为该行的非0元素个数

        // 去除边界面的数量（无邻单元需要减去）
        for (ULL cellFaceId : cell.getFaceIndexes())
        {
            if (faces[cellFaceId].getNeighborIndex() == -1)
            {
                --elementNumber;
            }
        }
    }
    values_.resize(elementNumber);
    // colIndex_.resize(elementNumber);
    colIndex_.reserve(elementNumber);

    // 给rowPointer_预分配内存
    // rowPointer_.resize(size_ + 1);
    // rowPointer_[0] = 0;
    rowPointer_.reserve(size_ + 1);
    rowPointer_.emplace_back(0);


    // 构造colIndex_与rowPointer_
    ULL currentCellIndex = 0;  // 记录当前遍历到第几个单元（第几行）
    std::priority_queue<ULL, std::vector<ULL>, std::greater<>> neighborCellQueue;       // 小根堆
    ULL neighborCellIndex = -1; // 存储当前单元的邻居单元号
    ULL currentColIndex = 0;    // 下次添加colIndex_的索引
    for (const Cell& cell : mesh->getCells())
    {
        // 通过便利cell的所有面找到与当前cell邻居单元号（需要按顺序）
        std::vector<ULL> FaceIndexes = cell.getFaceIndexes();
        for (const ULL& neighborFaceId : FaceIndexes)
        {
            if (faces[neighborFaceId].getNeighborIndex() == -1)  // 边界面没有邻单元，直接跳过。
            {
                continue;
            }

            const Face& neighborFace = faces[neighborFaceId];

            // 先确定邻单元号
            neighborCellIndex = (
                currentCellIndex == neighborFace.getOwnerIndex() ?
                neighborFace.getNeighborIndex() :
                neighborFace.getOwnerIndex()
                );
            neighborCellQueue.emplace(neighborCellIndex);
        }
        neighborCellQueue.emplace(currentCellIndex);   // 再将当前自己单元号入队

        // 此时队列中的数字为该行矩阵的非0元素的列索引
        // 按从小到大顺序放入colIndex_，并构造rowPointer_
        while (!neighborCellQueue.empty())  // 从小到大依次添加
        {
            // colIndex_[currentColIndex++] = neighborCellQueue.top();
            colIndex_.emplace_back(neighborCellQueue.top());
            neighborCellQueue.pop();
        }
        // rowPointer_[currentCellIndex + 1] = currentColIndex;
        rowPointer_.emplace_back(colIndex_.size());
        ++currentCellIndex;
    }
    // isValid_ = true;
}

template<typename Tp>
inline void SparseMatrix<Tp>::printMatrix() const
{
    if (!isValid_)      // 如果无效，则不能打印
    {
        std::cerr << "SparseMatrix<Tp>::printMatrix() Error: matrix is not valid" << std::endl;
    }

    for (ULL i = 0; i < size_; ++i)
    {
        // 打印本行第一个非0元素前的0，保留四位小数
        ULL j = 0;
        std::cout << std::setprecision(4) << std::fixed;
        for (; j < colIndex_[rowPointer_[i]]; j++)
        {
            std::cout << std::setw(10) << 0.0;
        }

        // 开始依次判断并打印本行非0元素和0元素
        for (ULL index = rowPointer_[i];
            index < rowPointer_[i + 1]; ++j)
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
