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
#include <unordered_set>

class Mesh
{
    using ULL = unsigned long long;
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
    void readBoundaryPatch(const std::string& boundaryPath);
    void readFaces(
        const std::string& facesPath,
        const std::string& ownerPath,
        const std::string& neighbourPath
    );
    void readNeighbour(const std::string& neighbourPath, std::vector<ULL> internalFaceIndices);
    void buildCellsFromFaces();
    BoundaryPatch::BoundaryType stringToType(const std::string& name) const;
    void calculateMeshInfo();


private:
    std::vector<Vector<Scalar>> points_;                                    // 点列表
    std::vector<Face> faces_;                                               // 面列表
    std::vector<Cell> cells_;                                               // 单元列表
    std::unordered_set<ULL> internalCellIndexes_;                                 // 内部单元索引列表
    std::unordered_set<ULL> boundaryCellIndexes_;                                 // 边界面索引列表
    std::unordered_map<std::string, BoundaryPatch> boundaryPatch_;     // 边界条件映射
    bool isValid_;                                                          // 网格是否有效
    std::string meshPath_;                                                  // 网格路径

};



#endif // MESH_H_