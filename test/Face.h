#ifndef FACE_H_
#define FACE_H_


#include <vector>
#include "Vector.hpp"

class Face
{
    using Scalar = double;
    using ULL = unsigned long long;
public:
    Face(int pointNum = 0);
    Face(const Face& face) = default;
    Face(Face&&) noexcept = default;
    Face& operator=(const Face& rhs);
    Face& operator=(Face&& rhs) noexcept = default;
    ~Face() = default;

public:
    // 获取面的点索引列表
    const std::vector<ULL>& getPointIndices() const;
    // 获取面的法向量
    const Vector<Scalar>& getNormal() const;
    // 获取面积
    Scalar getArea() const;
    // 获取面心
    const Vector<Scalar>& getCenter() const;
    // 获取点的数量
    int getPointNum() const;

private:
    // 计算Face的几何属性
    void calculateGeometry();

private:
    ULL id;
    std::vector<ULL> pointIndexes_;     // 构成面的点的索引
    ULL owner_;                         // 主单元
    ULL neighbor_;                      // 邻接单元
    Vector<Scalar> normal_;             // 法向量
    Scalar area_;                       // 面积
    Vector<Scalar> center_;             // 面心
};


#endif // FACE_H_

