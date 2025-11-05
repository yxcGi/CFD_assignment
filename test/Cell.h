#ifndef CELL_H_
#define CELL_H_


#include <vector>
#include "Vector.hpp"
class Cell
{
    using ULL = unsigned long long;
    using Scalar = double;
public:
    Cell(int faceNum = 0);
    Cell(const Cell&) = default;
    Cell(Cell&&) noexcept = default;
    Cell& operator=(const Cell& rhs);
    Cell& operator=(Cell&&) noexcept = default;

    ~Cell() = default;

public:
    // 获取单元的面索引列表
    const std::vector<ULL>& getFaceIndices() const;
    // 获取体积
    Scalar getVolume() const;
    // 获取单元中心
    const Vector<Scalar>& getCenter() const;
    // 获取面的数量
    int getFaceNum() const;
private:
    // 计算Cell的几何属性
    void calculateGeometry();

private:
    std::vector<ULL> faceIndexes_;      // 面的索引
    Scalar volume_;                     // 单元体积
    Vector<Scalar> center_;             // 单元中心
};

#endif // CELL_H_