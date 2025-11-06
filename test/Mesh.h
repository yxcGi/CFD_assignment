#ifndef MESH_H_
#define MESH_H_


#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include "Vector.hpp"
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"

class Mesh
{
    using Scalar = double;
public:
    Mesh();
    Mesh(const std::string& path);
    Mesh(const Mesh&) = delete;
    Mesh(Mesh&&) noexcept = default;
    Mesh& operator=(const Mesh&) = delete;
    Mesh& operator=(Mesh&&) noexcept = default;
    ~Mesh() = default;

public:
    // 常用接口
    void readMesh(const std::string& path);
    void writeMesh(const std::string& path);
    void printMeshInfo() const;

    // 获取器
    const std::vector<Vector<Scalar>>& getPoints() const;
    const std::vector<Face>& getFaces() const;
    const std::vector<Cell>& getCells() const;

private:
    // 私有处理接口
    void readPoints(const std::string& pointsPath);
    void readBoundaryConditions(const std::string& boundaryPath);
    void readFaces(
        const std::string& facesPath,
        const std::string& ownerPath,
        const std::string& neighbourPath
    );


private:
    std::vector<Vector<Scalar>> points_;                                    // 点列表
    std::vector<Face> faces_;                                               // 面列表
    std::vector<Cell> cells_;                                               // 单元列表
    std::unordered_map<std::string, BoundaryPatch> boundaryConditions_;     // 边界条件映射
    bool isValid_;                                                          // 网格是否有效
    std::string meshPath_;                                                  // 网格路径
};



#endif // MESH_H_