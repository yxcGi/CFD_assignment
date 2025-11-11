#include "Cell.h"
#include <unordered_set>

Cell::Cell()
{
    faceIndexes_.reserve(4);
}

const std::vector<Cell::ULL>& Cell::getFaceIndices() const
{
    return faceIndexes_;
}

Scalar Cell::getVolume() const
{
    return volume_;
}

const Vector<Scalar>& Cell::getCenter() const
{
    return center_;
}

int Cell::getFaceNum() const
{
    return faceIndexes_.size();
}

void Cell::addFaceIndex(ULL faceIndex)
{
    faceIndexes_.emplace_back(faceIndex);
}

void Cell::calculateCellInfo(
    const std::vector<Face>& faces,
    const std::vector<Point>& points
)
{
    std::unordered_set<ULL> points_set;
    if (faceIndexes_.size() == 4)   // 四面体直接，有计算公式
    {
        // 四面体只需遍历两个面即可将所有点集齐
        for (size_t i = 0; i < 2; ++i)
        {
            const Face& face = faces[faceIndexes_[i]];
            for (ULL pointIndex : face.getPointIndexes())
            {
                points_set.emplace(pointIndex);
            }
        }

        // 去除unordered_set里的元素放入vector
        std::vector<ULL> points_vec(
            std::make_move_iterator(points_set.begin()),
            std::make_move_iterator(points_set.end())
        );
        Point A = points[points_vec[0]];
        Point B = points[points_vec[1]];
        Point C = points[points_vec[2]];
        Point D = points[points_vec[3]];
        volume_ = std::abs((A - D) & ((B - D) ^ (C - D))) / 6.0;

        // 计算四面体形心
        Point bottomCenter = (B + C + D) / 3.0;     // 底面三角形形心
        // 四面体形心在底面形心与顶点的1/4处
        center_ = 0.75 * bottomCenter + 0.25 * A;
        return;
    }
    else if (faceIndexes_.size() > 4)   // 多面体
    {
        Scalar allVolume = 0.0;     // 总体积
        Point weightedCenter(0.0, 0.0, 0.0);    // 加权形心

        // 计算Cell中心（体心）
        Point volumeCenter(0.0, 0.0, 0.0);      // 体心（非形心）
        for (const Face& face : faces)      // 获取所有点
        {
            for (ULL pointIndex : face.getPointIndexes())
            {
                points_set.emplace(pointIndex);
            }
        }
        for (ULL pointIndex : points_set)
        {
            volumeCenter += points[pointIndex];
        }
        volumeCenter /= points_set.size();      // 计算体心


        for (const Face& face : faces)
        {
            Point faceCenter = face.getCenter();
            Point thisCenter = 0.75 * faceCenter + 0.25 * volumeCenter; // 本锥体形心

            Scalar thisVolume = std::abs(
                face.getArea() * (volumeCenter - faceCenter) & face.getNormal()
            ) / 3.0;
            allVolume += thisVolume;
            weightedCenter += thisCenter * thisVolume;
        }
        volume_ = allVolume;
        center_ = weightedCenter / allVolume;
        return;
    }
    std::cerr << "Error: cell has less than 4 faces" << std::endl;
    throw std::runtime_error("Error: cell has less than 4 faces");
}
