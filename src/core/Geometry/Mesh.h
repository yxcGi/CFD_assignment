#ifndef MESH_H_
#define MESH_H_


// #include <iostream>
#include <vector>
#include <unordered_map>
#include "Vector.hpp"
#include "Face.h"
#include "Cell.h"
#include "BoundaryPatch.h"
#include <unordered_set>
// #include "BaseField.hpp"
#include "FieldType.hpp"

// using field::FieldType;




class Mesh
{
    using ULL = unsigned long long;
    using Scalar = double;
public:
    // 网格维度
    enum class Dimension
    {
        TWO_D,
        THREE_D
    };
    // 网格形状（tecplot里的非结构网格类型）
    enum class MeshShape
    {
        TRIANGLE,   // 三角形单元
        QUADRILATERAL,  // 四边形单元
        TETRAHEDRON,   // 四面体单元
        BRICK,      // 六面体单元
    };

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
    void writeMeshToFile(const std::string& path) const;

    // 获取器
    const std::vector<Vector<Scalar>>& getPoints() const;
    const std::vector<Face>& getFaces() const;
    const std::vector<ULL>& getInternalFaceIndexes() const;
    const std::vector<ULL>& getBoundaryFaceIndexes() const;
    const std::vector<Cell>& getCells() const;
    ULL getCellNumber() const;
    ULL getFaceNumber() const;
    ULL getPointNumber() const;
    Dimension getDimension() const;
    MeshShape getMeshShape() const;

    ULL getInternalCellNumber() const;
    ULL getBoundaryFaceNumber() const;
    const std::pair<ULL, ULL>& getEmptyFaceIndexesPair() const;

    const std::unordered_map<std::string, BoundaryPatch>&
        getBoundaryPatches() const;

    // 用于获取不同类型场的数量的统一接口
    ULL getNumber(field::FieldType type) const;





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
    
    void writePoints(const std::string& pointsPath) const;
    void writeFaces(const std::string& facesPath) const;
    void writeOwnerAndNeighbour(const std::string& ownerPath, const std::string& neighbourPath) const;
    void writeBoundaryPatch(const std::string& boundaryPath) const;

    void buildCellsFromFaces();
    BoundaryPatch::BoundaryType stringToType(const std::string& name) const;
    void setMeshShape();
    void calculateMeshInfo();




private:
    std::vector<Vector<Scalar>> points_;                                    // 点列表
    std::vector<Face> faces_;                                               // 面列表
    std::vector<ULL> internalFaceIndexes_;                                 // 内部面索引列表
    std::vector<ULL> boundaryFaceIndexes_;                                 // 边界面索引列表
    std::vector<Cell> cells_;                                               // 单元列表
    std::unordered_set<ULL> internalCellIndexes_;                                 // 内部单元索引列表
    std::unordered_set<ULL> boundaryCellIndexes_;                                 // 边界面索引列表
    std::unordered_map<std::string, BoundaryPatch> boundaryPatches_;     // 边界条件映射
    std::pair<ULL, ULL> emptyFaceIndexesPair_;                              // empty边界面的起止索引（二维专属，左闭右开）
    bool isValid_;                                                          // 网格是否有效
    std::string meshPath_;                                                  // 网格路径
    Dimension dimension_;                                                   // 网格维度，默认三维
    MeshShape meshShape_;
};



#endif // MESH_H_