#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <string>
#include <map>
#include <cmath>

// 点结构
struct Point {
    double x, y;
    Point(double x = 0.0, double y = 0.0) : x(x), y(y) {}
};

// 单元结构
struct Cell {
    int id;
    std::vector<int> nodes;  // 节点索引
    std::vector<int> neighbors;  // 相邻单元索引
    Point centroid;  // 单元中心
    double area;  // 单元面积
    int type;  // 单元类型：3=三角形，4=四边形
};

// 面结构
struct Face {
    int id;
    int node1, node2;  // 面的两个节点
    int cell1, cell2;  // 面相邻的两个单元（边界面的cell2为-1）
    Point centroid;  // 面中心
    Point normal;  // 面法向量
    double length;  // 面长度
};

class Geometry {
private:
    std::vector<Point> nodes;
    std::vector<Cell> cells;
    std::vector<Face> faces;
    std::map<int, std::string> boundary_conditions;
    
    // 内部函数
    void readMSHFile(const std::string& filename);
    void computeCellProperties();
    void computeFaceProperties();
    void buildFaceList();
    void computeGradients();
    
public:
    Geometry();
    ~Geometry();
    
    // 主要接口
    void loadMesh(const std::string& filename);
    void printMeshInfo() const;
    
    // 获取器
    const std::vector<Point>& getNodes() const { return nodes; }
    const std::vector<Cell>& getCells() const { return cells; }
    const std::vector<Face>& getFaces() const { return faces; }
    int getNumNodes() const { return static_cast<int>(nodes.size()); }
    int getNumCells() const { return static_cast<int>(cells.size()); }
    int getNumFaces() const { return static_cast<int>(faces.size()); }
    
    // 几何计算
    Point getCellCentroid(int cellId) const;
    double getCellArea(int cellId) const;
    Point getFaceCentroid(int faceId) const;
    Point getFaceNormal(int faceId) const;
    double getFaceLength(int faceId) const;
    
    // 边界条件
    void setBoundaryCondition(int boundaryId, const std::string& bcType);
    std::string getBoundaryCondition(int boundaryId) const;
};

#endif // GEOMETRY_H
