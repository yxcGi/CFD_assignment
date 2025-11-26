#ifndef CELLFIELD_H_
#define CELLFIELD_H_





#include "BaseField.hpp"
#include "FaceField.hpp"
#include "interpolation.hpp"


template<typename Tp>
class CellField : public BaseField<Tp>
{
    using ULL = unsigned long long;
public:
    CellField() = delete;
    CellField(CellField<Tp>&&) noexcept = default;
    // 构造但不初始化，场无效                     
    CellField(const std::string& name, Mesh* mesh);
    // 构造并初始化，场有效
    CellField(const std::string& name, Mesh* mesh, const Tp& initialValue);
    CellField<Tp>& operator=(const CellField<Tp>& other);
    CellField<Tp>& operator=(CellField<Tp>&&) noexcept = default;
public:
    // 从单元向面插值（默认线性插值）
    FaceField<Tp> cellToFace(interpolation::Scheme scheme = interpolation::Scheme::LINEAR) const;

private:
};

#pragma 函数实现
template<typename Tp>
inline CellField<Tp>::CellField(const std::string& name, Mesh* mesh)
    : BaseField<Tp>(name, mesh)
{
    this->type_ = field::FieldType::CELL_FIELD;
}

template<typename Tp>
inline CellField<Tp>::CellField(const std::string& name, Mesh* mesh, const Tp& initialValue)
    : BaseField<Tp>(name, mesh)
{
    this->type_ = field::FieldType::CELL_FIELD;
    // 初始化大小
    ULL cellNum = this->getDataNumer();
    this->data_.resize(cellNum, initialValue);
    this->isValid_ = true;
}

template<typename Tp>
inline CellField<Tp>& CellField<Tp>::operator=(const CellField<Tp>& other)
{
    if (this != &other)
    {
        BaseField<Tp>::operator=(other);  // 调用基类的赋值运算符
    }
    return *this;
}

template<typename Tp>
inline FaceField<Tp> CellField<Tp>::cellToFace(interpolation::Scheme scheme) const
{
    using LL = long long;

    Mesh* mesh = this->mesh_;   // 网格
    const std::vector<Face>& faces = mesh->getFaces();  // 面列表
    const std::vector<Cell>& cells = mesh->getCells();
    // ULL faceNum = mesh->getFaceNumber();

    // 创建面场
    FaceField<Tp> faceField(this->name_, mesh);


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
        Tp ownerValue = this->data_[ownerIndex];
        Tp neighborValue = this->data_[neighborIndex];

        Tp faceValue = interpolator_(ownerValue, neighborValue, scheme, alpha);

        // 设置面场
        faceField.setValue(internalFaceIndex, faceValue);
    }

    // 再遍历边界面，需考虑边界条件
    // using BoundaryConditionMap = std::unordered_map<std::string, BoundaryCondition<T>>;
    for (const auto& [name, bc] : this->boundaryConditions_)
    {
        // 如果是empty，则跳过
        if (bc.getType() == BoundaryPatch::BoundaryType::EMPTY)
        {
            continue;
        }

        ULL nFace = bc.getNFace();
        ULL startFace = bc.getStartFace();



        for (ULL boundaryFaceIndex = startFace;
            boundaryFaceIndex < startFace + nFace;
            ++boundaryFaceIndex)
        {
            const Face& face = faces[boundaryFaceIndex];
            ULL ownerIndex = face.getOwnerIndex();
            const Cell& ownerCell = cells[ownerIndex];
            const Vector<Tp>& faceCenter = face.getCenter();
            const Vector<Tp>& ownerCenter = ownerCell.getCenter();
            const Vector<Scalar>& normal = face.getNormal();   // 面法向量
            Vector<Scalar> V_CB = faceCenter - ownerCenter;

            // 计算中间量, normal = E + T
            Vector<Scalar> E = V_CB / (V_CB & normal);
            Vector<Scalar> T = normal - E;
            Scalar E_magnitude = E.magnitude();
            Scalar distance_CB = ownerCenter.getDistance(faceCenter);
            Scalar a = bc.get_a();
            Scalar b = bc.get_b();
            const Tp& c = bc.get_c();
            // const auto& gradientCell =       // 计算梯度，挖坑


            // 计算c1, c2
            auto c1 = (b * E_magnitude) / (a * distance_CB + b * E_magnitude);
            // auto c2 = (c - b * )    // 需要梯度计算，挖坑


            // faceField.setValue(boundaryFaceIndex, bc.)
        }
    }
}




#endif // CELLFIELD_H_
