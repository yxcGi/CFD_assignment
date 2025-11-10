#ifndef FACE_H_
#define FACE_H_


#include <vector>
#include "Vector.hpp"
#include <iostream>

class Face
{
    using Scalar = double;
    using ULL = unsigned long long;
    using LL = long long;
public:
    Face() = delete;
    Face(const std::vector<ULL>& pointIndexs, ULL ownerIndex, ULL neighborIndex = -1);
    Face(const Face& face) = default;
    Face(Face&&) noexcept = default;
    Face& operator=(const Face& rhs);
    Face& operator=(Face&& rhs) noexcept = default;
    ~Face() = default;
    

public:
    // 获取面的点索引列表
    const std::vector<ULL>& getPointIndexes() const;
    // 获取面的法向量
    const Vector<Scalar>& getNormal() const;
    // 获取面积
    Scalar getArea() const;
    // 获取面心
    const Vector<Scalar>& getCenter() const;
    // 获取点的数量
    int getPointNum() const;
    // 获取主单元索引
    ULL getOwnerIndex() const;
    // 获取领单元索引
    LL getNeighborIndex() const;
    // 设置邻单元
    void setNeighbor(ULL neighborIndex);
    
    // 计算面信息，构造时调用
    void calculateFaceInfo(const std::vector<Vector<Scalar>>& points);

    // 输出流重载
    friend std::ostream& operator<<(std::ostream& out, const Face& face);


    // 测试用
    void printFaceInfo() const
    {
        std::cout << "Face Info:" << std::endl;
        std::cout << "  Point Num: " << pointNum_ << std::endl;
        std::cout << "  Point Indices: ";
        for (const auto& idx : pointIndexes_)
        {
            std::cout << idx << " ";
        }
        std::cout << std::endl;
        std::cout << "  Owner: " << owner_ << std::endl;
        std::cout << "  Neighbor: " << neighbor_ << std::endl;
        std::cout << "  Area: " << area_ << std::endl;
        std::cout << "  Center: " << center_ << std::endl;
        std::cout << "  Normal: " << normal_ << std::endl;
    }

private:
    // 计算面积辅助函数
    Scalar calculateArea(const std::vector<Vector<Scalar>>& points) const;
    // 计算坐标三角形中心辅助函数
    Point calculateCenter(const Point& p1, const Point& p2, const Point& p3) const;

private:
    // ULL id;
    int pointNum_;                      // 点的数量
    std::vector<ULL> pointIndexes_;     // 构成面的点的索引
    ULL owner_;                         // 主单元
    LL neighbor_;                       // 邻接单元
    Vector<Scalar> normal_;             // 法向量
    Scalar area_;                       // 面积
    Vector<Scalar> center_;             // 面心
};


#endif // FACE_H_

