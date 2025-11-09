#pragma once


#include "BaseField.hpp"
#include "FaceField.hpp"
#include "interpolation.hpp"


template<typename T>
class CellField : public BaseField<T>
{
    using ULL = unsigned long long;
public:
    CellField() = delete;
    // 构造但不初始化，场无效                     
    CellField(const std::string& name, Mesh* mesh);
    // 构造并初始化，场有效
    CellField(const std::string& name, Mesh* mesh, const T& initialValue);

public:
    // 从单元向面插值（默认线性插值）
    FaceField<T> cellToFace(interpolation::Scheme scheme = interpolation::Scheme::LINEAR) const;

private:
    Interpolation<T> interpolator_;
};

#pragma 函数实现

template<typename T>
inline CellField<T>::CellField(const std::string& name, Mesh* mesh)
    : BaseField<T>(name, mesh)
{
    this->type_ = field::FieldType::CELL_FIELD;
}

template<typename T>
inline CellField<T>::CellField(const std::string& name, Mesh* mesh, const T& initialValue)
    : BaseField<T>(name, mesh)
{
    this->type_ = field::FieldType::CELL_FIELD;
    // 初始化大小
    ULL cellNum = this->getDataNumer();
    this->data_.resize(cellNum, initialValue);
    this->isValid_ = true;
}

template<typename T>
inline FaceField<T> CellField<T>::cellToFace(interpolation::Scheme scheme) const
{
    using LL = long long;

    Mesh* mesh = this->mesh_;   // 网格
    const std::vector<Face>& faces = mesh->getFaces();  // 面列表
    const std::vector<Cell>& cells = mesh->getCells();
    // ULL faceNum = mesh->getFaceNumber();


    // 先遍历内部面
    std::vector<ULL> internalFaceIndexes = mesh->getInternalFaceIndexes();

    for (const ULL internalFaceIndex : internalFaceIndexes)
    {
        // 获取当前面，以及其相邻单元索引
        const Face& face = faces[internalFaceIndex];
        ULL ownerIndex = face.getOwnerIndex();
        LL neighborIndex = face.getNeighborIndex();

        // 获取当前面，以及其相邻单元中心坐标
        Vector<Scalar> faceCenter = face.getCenter();
        Vector<Scalar> ownerCenter = cells[ownerIndex].getCenter();
        Vector<Scalar> neighborCenter = cells[neighborIndex].getCenter();

        // 获取面的法向量
        Vector<Scalar> faceNormal = face.getNormal();

        // 计算面到两个面的距离
        Scalar ownerDistance = (faceCenter - ownerCenter) & faceNormal;
        Scalar neighborDistance = (faceCenter - neighborCenter) & faceNormal;

        // 计算插值权重
        Scalar alpha = ownerDistance / (ownerDistance + neighborDistance);

        // 插值
        T ownerValue = this->data_[ownerIndex];
        T neighborValue = this->data_[neighborIndex];

        T faceValue = interpolator_(ownerValue, neighborValue, scheme, alpha);
    }

    // 再遍历边界面，需考虑边界条件    挖坑
    std::vector<ULL> boundaryFaceIndexes = mesh->getBoundaryFaceIndexes();
}




