#ifndef CELL_H_
#define CELL_H_


#include <vector>
#include "Vector.hpp"
#include <Face.h>



class Cell
{
    using ULL = unsigned long long;
    using Scalar = double;
public:

    Cell();
    Cell(const Cell&) = default;
    Cell(Cell&&) noexcept = default;
    Cell& operator=(const Cell& rhs);
    Cell& operator=(Cell&&) noexcept = default;

    ~Cell() = default;

public:
    // 获取单元的面索引列表
    const std::vector<ULL>& getFaceIndexes() const;
    const std::vector<ULL>& getPointIndexes() const;
    // 获取体积
    Scalar getVolume() const;
    // 获取单元中心
    const Vector<Scalar>& getCenter() const;
    // 获取面的数量
    int getFaceNum() const;
    // 添加face索引
    void addFaceIndex(ULL faceIndex);

    // 计算单元的几何属性
    void calculateCellInfo(
        const std::vector<Face>& faces,     // mesh里的面列表
        const std::vector<Point>& points    // mesh里的点列表
    );


    // 打印测试用
    void printCellInfo() const
    {
        std::cout << "Cell Info:" << std::endl;
        std::cout << "  Face Num: " << faceIndexes_.size() << std::endl;
        std::cout << "  Face Indices: ";
        for (const auto& idx : faceIndexes_)
        {
            std::cout << idx << " ";
        }
        std::cout << std::endl;
        std::cout << "  Volume: " << volume_ << std::endl;
        std::cout << "  Center: " << center_ << std::endl;
    }


private:
    // 计算Cell的几何属性

private:
    ULL id;                             // 单元ID
    std::vector<ULL> faceIndexes_;      // 面的索引
    std::vector<ULL> pointIndexes_;     // 点的索引
    Scalar volume_;                     // 单元体积
    Vector<Scalar> center_;             // 单元中心
};

#endif // CELL_H_