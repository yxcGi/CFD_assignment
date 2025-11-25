#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_



#include <vector>
#include "Mesh.h"
#include <queue>
#include <fstream>

/**
 * @brief 采用CSR方法存储稀疏矩阵，仅为方阵
 * @tparam Tp 一般为double，int 不会是Vector
 */
template <typename Tp>
class SparseMatrix
{
    using ULL = unsigned long long;
    constexpr static double EPSILON = 1e-10;
    inline static int WIDTH = 10;       // 输出宽度 
    // 设置小数位数
    inline static int PRECISION = 4;    // 保留小数位数
public:
    SparseMatrix() = default;
    SparseMatrix(const std::vector<std::vector<Tp>>& matrix);   // 直接用二维数组初始化(拷贝一份二维数组，还未压缩)
    SparseMatrix(Mesh* mesh);     // 通过读取网格信息初始化(给相应位置留空位)

public:
    // 设置输出宽度
    static void setWidth(int width = 10);
    // 设置小数保留位数
    static void setPrecision(int precision = 4);

public:
    // 初始化函数
    void compress();
    void init(const Mesh* mesh);
    void init(const std::vector<std::vector<Tp>>& matrix);
    // 打印
    void printMatrix() const;
    // 设置系数矩阵i, j元素值
    void setValue(ULL i, ULL j, Tp value);
    // 添加矩阵元素
    void addValue(ULL i, ULL j, Tp value);
    // 查看特定位置元素，只读
    Tp at(ULL i, ULL j) const;


    // 获取矩阵元素，括号重载
    Tp& operator()(ULL i, ULL j);
    const Tp& operator()(ULL i, ULL j) const;

private:


private:
    std::vector<Tp> values_;        // 存储数据（行优先），大小为矩阵非0元素个数
    std::vector<ULL> colIndexs_;     // 每个元素列索引, 与values一一对应
    std::vector<ULL> rowPointer_;   // 行指针，每一行起始元素的索引，大小为矩阵行数
    ULL size_;                      // 矩阵大小
    std::vector<std::vector<Tp>> unCompressedMatrix_;   // 未压缩的矩阵
    Mesh* mesh_{ nullptr };
    bool isValid_{ false };         // 是否有效
    bool isCompressed_{ false };    // 是否压缩
};


#pragma region 函数实现


template<typename Tp>
inline SparseMatrix<Tp>::SparseMatrix(const std::vector<std::vector<Tp>>& matrix)
    : size_(matrix.size())
    , unCompressedMatrix_(matrix)
{
    // 检查是否为方阵，非方阵抛出异常
    for (const std::vector<Tp>& row : matrix)
    {
        if (row.size() != size_)
        {
            std::cerr << "SparseMatrix<Tp>::SparseMatrix(const std::vector<std::vector<Tp>>& matrix) Error: matrix is not square" << std::endl;
            throw std::invalid_argument("matrix is not square");
        }
    }
    isValid_ = true;
}


template<typename Tp>
inline SparseMatrix<Tp>::SparseMatrix(Mesh* mesh)
    : mesh_(mesh)
{
    init(mesh);
}

template<typename Tp>
inline void SparseMatrix<Tp>::setWidth(int width)
{
    if (width <= 0)
    {
        std::cerr << "SparseMatrix<Tp>::setWidth(int width) Error: width must be greater than 0" << std::endl;
        throw std::invalid_argument("width must be greater than 0");
    }

    if (width >= 30)
    {
        std::cerr << "SparseMatrix<Tp>::setWidth(int width) Error: width must be less than 30" << std::endl;
        throw std::invalid_argument("width must be less than 30");
    }
    WIDTH = width;
}

template<typename Tp>
inline void SparseMatrix<Tp>::setPrecision(int precision)
{
    if (precision <= 0)
    {
        std::cerr << "SparseMatrix<Tp>::setPrecision(int precision) Error: precision must be greater than 0" << std::endl;
        throw std::invalid_argument("precision must be greater than 0");
    }

    if (precision >= 10)
    {
        std::cerr << "SparseMatrix<Tp>::setPrecision(int precision) Error: precision must be less than 10" << std::endl;
        throw std::invalid_argument("precision must be less than 10");
    }

    PRECISION = precision;
}

template<typename Tp>
inline void SparseMatrix<Tp>::compress()
{
    if (!isValid_)  // 矩阵无效
    {
        std::cerr << "SparseMatrix<Tp>::init(const std::vector<std::vector<Tp>>& matrix) Error: matrix is not valid" << std::endl;
        throw std::invalid_argument("matrix is not valid");
    }

    if (isCompressed_)  // 矩阵已压缩则不能压缩
    {
        std::cerr << "SparseMatrix<Tp>::init(const std::vector<std::vector<Tp>>& matrix) Error: matrix is already compressed" << std::endl;
        throw std::invalid_argument("matrix is already compressed");
    }


    // 开始对二维数数组进行压缩
    rowPointer_.reserve(size_ + 1);

    // 开始构造
    rowPointer_.push_back(0);
    for (ULL i = 0; i < size_; ++i)
    {
        for (ULL j = 0; j < size_; ++j)
        {
            const Tp ele = unCompressedMatrix_[i][j];
            if (std::abs(ele) > EPSILON)   // 非0元素
            {
                values_.emplace_back(ele);
                colIndexs_.emplace_back(j);
            }
        }
        // 每行结束后记录下一行的起始值（对于空行依然，成立）
        rowPointer_.emplace_back(colIndexs_.size());
    }
    isCompressed_ = true;
    std::cout << "The matrix has been compressed" << std::endl;
}

template<typename Tp>
inline void SparseMatrix<Tp>::init(const Mesh* mesh)
{
    if (isValid_)
    {
        std::cerr << "SparseMatrix<Tp>::init(const Mesh* mesh) Error: matrix is already initialized" << std::endl;
        throw std::invalid_argument("matrix is already initialized");
    }

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
    values_.resize(elementNumber);      // 便于之后访问
    // colIndexs_.resize(elementNumber);
    colIndexs_.reserve(elementNumber);

    // 给rowPointer_预分配内存
    // rowPointer_.resize(size_ + 1);
    // rowPointer_[0] = 0;
    rowPointer_.reserve(size_ + 1);
    rowPointer_.emplace_back(0);


    // 构造colIndex_与rowPointer_
    ULL currentCellIndex = 0;  // 记录当前遍历到第几个单元（第几行）
    std::priority_queue<ULL, std::vector<ULL>, std::greater<>> neighborCellQueue;       // 小根堆
    ULL neighborCellIndex = -1; // 存储当前单元的邻居单元号
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
            // colIndexs_[currentColIndex++] = neighborCellQueue.top();
            colIndexs_.emplace_back(neighborCellQueue.top());
            neighborCellQueue.pop();
        }
        // rowPointer_[currentCellIndex + 1] = currentColIndex;
        rowPointer_.emplace_back(colIndexs_.size());
        ++currentCellIndex;
    }
    isValid_ = true;
    isCompressed_ = true;
}

template<typename Tp>
inline void SparseMatrix<Tp>::init(const std::vector<std::vector<Tp>>& matrix)
{
    size_ = matrix.size();
    unCompressedMatrix_ = matrix;

    // 检查是否为方阵，非方阵抛出异常
    for (const std::vector<Tp>& row : matrix)
    {
        if (row.size() != size_)
        {
            std::cerr << "SparseMatrix<Tp>::SparseMatrix(const std::vector<std::vector<Tp>>& matrix) Error: matrix is not square" << std::endl;
            throw std::invalid_argument("matrix is not square");
        }
    }
    isValid_ = true;
}

template<typename Tp>
inline void SparseMatrix<Tp>::printMatrix() const
{
    if (!isValid_)      // 如果无效，则不能打印
    {
        std::cerr << "SparseMatrix<Tp>::printMatrix() Error: matrix is not valid" << std::endl;
        throw std::invalid_argument("matrix is not valid");
    }

    if (size_ > 30) // 矩阵过大则不支持输出
    {
        std::cerr << "SparseMatrix<Tp>::printMatrix() Error: matrix is too large" << std::endl;
        return;
    }

    // 存储当前输出格式，便于打印完后恢复
    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);


    std::cout << std::fixed << std::setprecision(PRECISION);    // 设置保留小数位数
    if (!isCompressed_) // 如果未压缩直接通过二维数组进行打印
    {
        for (const std::vector<Tp>& row : unCompressedMatrix_)
        {
            for (Tp ele : row)
            {
                std::cout << std::setw(WIDTH) << ele;
            }
            std::cout << "\n";
        }
    }
    else    // 压缩矩阵打印
    {
        for (ULL i = 0; i < size_; ++i)
        {
            // 判断本行是否为空行(全0), 若为空行直接打印本行
            if (rowPointer_[i] == rowPointer_[i + 1])   // 空行
            {
                for (ULL j = 0; j < size_; ++j)
                {
                    std::cout << std::setw(WIDTH) << 0.0;
                }
                std::cout << "\n";
            }
            else        // 非空行f
            {
                // 打印本行第一个非0元素前的0，保留四位小数
                ULL j = 0;
                for (; j < colIndexs_[rowPointer_[i]]; j++)
                {
                    std::cout << std::setw(WIDTH) << 0.0;
                }

                // 开始依次判断并打印本行非0元素和0元素
                for (ULL index = rowPointer_[i];
                    index < rowPointer_[i + 1]; ++j)
                {
                    if (j == colIndexs_[index])  // 该列为非0元素
                    {
                        std::cout << std::setw(WIDTH) << values_[index];
                        ++index;
                    }
                    else
                    {
                        std::cout << std::setw(WIDTH) << 0.0;
                    }
                }

                // 填充本行后面的0
                for (; j < size_; ++j)
                {
                    std::cout << std::setw(WIDTH) << 0.0;
                }
                std::cout << "\n";
            }
        }
    }
    std::cout.copyfmt(oldState);    // 恢复输出格式
}

template<typename Tp>
inline void SparseMatrix<Tp>::setValue(ULL i, ULL j, Tp value)
{
    (*this)(i, j) = value;
}

template<typename Tp>
inline void SparseMatrix<Tp>::addValue(ULL i, ULL j, Tp value)
{
    (*this)(i, j) += value;
}

template<typename Tp>
inline Tp SparseMatrix<Tp>::at(ULL i, ULL j) const
{
    if (!isValid_)
    {
        std::cerr << "SparseMatrix<Tp>::at(ULL i, ULL j) Error: matrix is not valid" << std::endl;
        throw std::invalid_argument("matrix is not valid");
    }

    if (i >= size_ || j >= size_)
    {
        std::cerr << "SparseMatrix<Tp>::at(ULL i, ULL j) Error: index out of range" << std::endl;
        throw std::out_of_range("index out of range");
    }

    if (!isCompressed_)
    {
        return unCompressedMatrix_[i][j];
    }
    else
    {
        // 判断是否存在第i行元素
        if (rowPointer_[i] == rowPointer_[i + 1])   // 该行不存在元素，无法设置
        {
            std::cerr << "SparseMatrix<Tp>::at(ULL i, ULL j) Error: row " << i << " is empty" << std::endl;
            throw std::invalid_argument("row " + std::to_string(i) + " is empty");
        }

        // 寻找第j列元素
        for (ULL index = rowPointer_[i]; index < rowPointer_[i + 1]; ++index)
        {
            if (colIndexs_[index] == j)
            {
                return values_[index];
            }
        }

        // 若未找到则返回0
        return Tp(0);
    }
}

template<typename Tp>
inline Tp& SparseMatrix<Tp>::operator()(ULL i, ULL j)
{
    // 判断矩阵是否有效
    if (!isValid_)
    {
        std::cerr << "SparseMatrix<Tp>::operator()(ULL i, ULL j) Error: matrix is not valid" << std::endl;
        throw std::invalid_argument("matrix is not valid");
    }

    // 判断索引是否越界
    if (i >= size_ || j >= size_)
    {
        std::cerr << "SparseMatrix<Tp>::operator()(ULL i, ULL j) Error: index out of range" << std::endl;
        throw std::out_of_range("index out of range");
    }

    if (!isCompressed_) // 未被压缩直接修改底层二维数组
    {
        return unCompressedMatrix_[i][j];
    }
    else
    {
        // 判断是否存在第i行元素
        if (rowPointer_[i] == rowPointer_[i + 1])   // 该行不存在元素，无法设置
        {
            std::cerr << "SparseMatrix<Tp>::setValue(ULL i, ULL j, Tp value) Error: row " << i << " is empty" << std::endl;
            throw std::invalid_argument("row " + std::to_string(i) + " is empty");
        }

        // 寻找第j列元素
        for (ULL index = rowPointer_[i]; index < rowPointer_[i + 1]; ++index)
        {
            if (colIndexs_[index] == j)
            {
                return values_[index];
            }
        }

        // 走到这说明不存在j列元素
        std::cerr << "SparseMatrix<Tp>::operator()(ULL i, ULL j) Error: column " << j << " is not exist" << std::endl;
        throw std::invalid_argument("column " + std::to_string(j) + " is not exist");
    }
}

template<typename Tp>
inline const Tp& SparseMatrix<Tp>::operator()(ULL i, ULL j) const
{
    if (!isValid_)
    {
        std::cerr << "SparseMatrix<Tp>::operator()(ULL i, ULL j) Error: matrix is not valid" << std::endl;
        throw std::invalid_argument("matrix is not valid");
    }

    if (i >= size_ || j >= size_)
    {
        std::cerr << "SparseMatrix<Tp>::operator()(ULL i, ULL j) Error: index out of range" << std::endl;
        throw std::out_of_range("index out of range");
    }

    if (!isCompressed_) // 未被压缩直接返回底层二维数组元素
    {
        return unCompressedMatrix_[i][j];
    }
    else    // 已被压缩只能返回非0元素，否则抛出异常
    {
        // 是否存在该行元素
        if (rowPointer_[i] == rowPointer_[i + 1])
        {
            std::cerr << "SparseMatrix<Tp>::operator()(ULL i, ULL j) Error: row " << i << " is empty" << std::endl;
            throw std::invalid_argument("row " + std::to_string(i) + " is empty");
        }

        // 遍历该行元素
        for (ULL index = rowPointer_[i]; index < rowPointer_[i + 1]; ++index)
        {
            if (colIndexs_[index] == j)
            {
                return values_[index];
            }
        }

        // 不存在该j列元素
        std::cerr << "SparseMatrix<Tp>::operator()(ULL i, ULL j) Error: column " << j << " is not exist" << std::endl;
        throw std::invalid_argument("column " + std::to_string(j) + " is not exist");
    }
}




#endif // SPARSEMATRIX_H_
