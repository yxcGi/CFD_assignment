#include "Face.h"

using ULL = unsigned long long;
using LL = long long;
using Scalar = double;

Face::Face(
    std::vector<ULL>& pointIndexs,
    ULL ownerIndex,
    ULL neighborIndex
)
    : pointNum_(pointIndexs.size())
    , owner_(ownerIndex)
    , neighbor_(neighborIndex)
{
    pointIndexes_ = std::move(pointIndexs);
}

std::ostream& operator<<(std::ostream& out, const Face& face)
{
    out << face.pointNum_ << "(";
    for (auto ele : face.pointIndexes_)
    {
        std::cout << ele << " ";
    }
    std::cout << ")";
    return out;
}

const std::vector<ULL>& Face::getPointIndexes() const
{
    return pointIndexes_;
}

int Face::getPointNum() const
{
    return pointNum_;
}

ULL Face::getOwnerIndex() const
{
    return owner_;
}

LL Face::getNeighborIndex() const
{
    return neighbor_;
}

void Face::setNeighbor(ULL neighborIndex)
{
    neighbor_ = neighborIndex;
}

void Face::calculateFaceInfo(const std::vector<Vector<Scalar>>& points)
{

    // 计算所有点的坐标平均值
    Point center(0.0, 0.0, 0.0);
    for (size_t i = 0; i < pointNum_; ++i)
    {
        center += points[pointIndexes_[i]];
    }
    center /= pointNum_;

    Scalar allArea = 0.0;   // 总面积


    // 计算面积
    area_ = calculateArea(points);


    // 计算面心
    if (pointNum_ == 3)     // 三角形直接算平均
    {
        // 单元中心
        center_ = center;
    }
    else                    // 四边形以上需要加权
    {
        Point weightedCenter(0.0, 0.0, 0.0);
        for (size_t i = 0; i < pointNum_; ++i)
        {
            Point center_i = calculateCenter(
                points[pointIndexes_[i]],
                points[pointIndexes_[(i + 1) % pointNum_]],
                center
            );
            Vector<Scalar> vector1 = points[pointIndexes_[i]] - center;
            Vector<Scalar> vector2 = points[pointIndexes_[(i + 1) % pointNum_]] - center;
            Scalar area_i = 0.5 * (vector1 ^ vector2).magnitude();
            allArea += area_i;
            weightedCenter += area_i * center_i;
        }
        center_ = weightedCenter / allArea;
    }


    // 计算法向量，采用面积加权法
    Vector<Scalar> weightedNormalVector(0.0, 0.0, 0.0);
    if (pointNum_ == 3)     // 三角形不用加权
    {
        const Point& point0 = points[pointIndexes_[0]];
        const Point& point1 = points[pointIndexes_[1]];
        const Point& point2 = points[pointIndexes_[2]];

        // 两条边的向量
        Vector<Scalar> vector1 = point1 - point0;
        Vector<Scalar> vector2 = point2 - point0;

        normal_ = (vector1 ^ vector2).unitVector();
    }
    else                    // 四边形以上需要加权
    {
        for (size_t i = 0; i < pointNum_; ++i)
        {
            const Point& pointA = points[pointIndexes_[i]];
            const Point& pointB = points[pointIndexes_[(i + 1) % pointNum_]];

            Vector<Scalar> vector1 = pointA - center;
            Vector<Scalar> vector2 = pointB - center;

            Scalar area_i = 0.5 * (vector1 ^ vector2).magnitude();
            Vector<Scalar> normal_i = (vector1 ^ vector2).unitVector();
            weightedNormalVector += area_i * normal_i;
        }
        normal_ = weightedNormalVector.unitVector();
    }

    // 测试
    // printFaceInfo();
}

Scalar Face::calculateArea(const std::vector<Vector<Scalar>>& points) const
{
    Scalar area = 0.0;
    // 计算面积
    if (points.size() == 3)     // 三角形随便算
    {
        const Point& point0 = points[pointIndexes_[0]];
        const Point& point1 = points[pointIndexes_[1]];
        const Point& point2 = points[pointIndexes_[2]];

        // 两条边的向量
        Vector<Scalar> vector1 = point1 - point0;
        Vector<Scalar> vector2 = point2 - point0;

        // 三角形面积公式
        area = 0.5 * (vector1 ^ vector2).magnitude();
        return area;
    }

    // 计算单元中心
    Point center(0.0, 0.0, 0.0);
    for (size_t i = 0; i < pointNum_; ++i)
    {
        center += points[pointIndexes_[i]];
    }
    center /= pointNum_;

    // 计算面积，以中心向周围点连线构成三角形，分别求面积，再求和
    for (size_t i = 0; i < pointNum_; ++i)
    {
        const Point& pointA = points[pointIndexes_[i]];
        const Point& pointB = points[pointIndexes_[(i + 1) % pointNum_]];

        Vector<Scalar> vector1 = pointA - center;
        Vector<Scalar> vector2 = pointB - center;

        area += 0.5 * (vector1 ^ vector2).magnitude();
    }
    return area;
}

Point Face::calculateCenter(const Point& p1, const Point& p2, const Point& p3) const
{
    return Point(
        (p1.x() + p2.x() + p3.x()) / 3.0,
        (p1.y() + p2.y() + p3.y()) / 3.0,
        (p1.z() + p2.z() + p3.z()) / 3.0
    );
}

