#pragma once

#include "FaceField.hpp"
#include "CellField.hpp"
#include "FieldOperators/Gradient.hpp"


// 采用模板全特化，判断梯度类型
template<typename T>
struct GradientType;

template<>
struct GradientType<Scalar>     
{
    using Type = Vector<Scalar>;
};

template<>
struct GradientType<Vector<Scalar>>
{
    using Type = Tensor<Scalar>;
};


template <typename Tp>
class Field
{
    using ULL = unsigned long long;
public:
    Field() = delete;
    Field(const std::string& name, Mesh* mesh);
    Field(const Field<Tp>&) = delete;
    Field(Field<Tp>&&) noexcept = default;
    Field<Tp>& operator=(const Field<Tp>&) = delete;
    Field<Tp>& operator=(Field<Tp>&&) noexcept = default;
    ~Field() = default;

public:
    /* =========赋值======== */ // 与BaseField，接口一样
    // 统一赋值
    void setValue(const Tp& value);
    // 采用函数对象
    void setValue(const std::function<Tp(Scalar, Scalar, Scalar)>& func);


    // 获取器
    const FaceField<Tp>& getFaceField() const;
    const CellField<Tp>& getCellField() const;
    std::string name() const;
    Mesh* getMesh() const;

    // 场是否有效
    bool isValid() const;


    // cell场到face场的差值
    void cellToFace(interpolation::Scheme scheme = interpolation::Scheme::LINEAR);



    // 设置计算梯度方法
    void setGradientMethod(GradientMethod method);

    /* --------设置边界条件-------- */

    /* --------设置边界条件-------- */
    // a * φ + b * ∂φ/∂n = c
    void setBoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c);



private:
    FaceField<Tp> faceField_;        // 当前面场值
    CellField<Tp> cellField_;        // 当前单元场值
    CellField<Tp> cellField_0_;      // 上一步单元场的值
    // CellField<typename GradientType<Tp>::Type> cellGradientField_; // 单元场梯度
    CellField<decltype(Tp() * Vector<Scalar>())> cellGradientField_; // 单元场梯度
    std::string name_;
    GradientMethod gradientMethod_{GradientMethod::GAUSS_GREEN};
};

#pragma region 函数实现

template<typename Tp>
inline Field<Tp>::Field(const std::string& name, Mesh* mesh)
    : faceField_(name, mesh)
    , cellField_(name, mesh)
    , name_(name)
{}

template<typename Tp>
inline void Field<Tp>::setValue(const Tp& value)
{
    faceField_.setValue(value);
    cellField_.setValue(value);
}

template<typename Tp>
inline void Field<Tp>::setValue(const std::function<Tp(Scalar, Scalar, Scalar)>& func)
{
    faceField_.setValue(func);
    cellField_.setValue(func);
}


template<typename Tp>
inline const FaceField<Tp>& Field<Tp>::getFaceField() const
{
    return faceField_;
}

template<typename Tp>
inline const CellField<Tp>& Field<Tp>::getCellField() const
{
    return cellField_;
}

template<typename Tp>
inline std::string Field<Tp>::name() const
{
    return name_;
}

template<typename Tp>
inline Mesh* Field<Tp>::getMesh() const
{
    return faceField_.getMesh();
}

template<typename Tp>
inline bool Field<Tp>::isValid() const
{
    return faceField_.isValid() && cellField_.isValid();
}

template<typename Tp>
inline void Field<Tp>::cellToFace(interpolation::Scheme scheme)
{
    using LL = long long;
    Mesh* mesh = this->getMesh();      // 获取网格
    const std::vector<Face>& faces = mesh->getFaces();  // 面列表
    const std::vector<Cell>& cells = mesh->getCells();

    // 先遍历内部面
    std::vector<ULL> internalFaceIndexes = mesh->getInternalFaceIndexes();
    const CellField<Tp>& cellField = this->getCellField();

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
        Tp ownerValue = cellField[ownerIndex];
        Tp neighborValue = cellField[neighborIndex];

        Tp faceValue = interpolator_(ownerValue, neighborValue, scheme, alpha);

        // 设置面场
        faceField_.setValue(internalFaceIndex, faceValue);
    }

    // 再遍历边界面，需考虑边界条件
    for (const auto& [name, bc] : cellField_0_.getBoundaryConditions())
    {
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
            Scalar c1 = (b * E_magnitude) / (a * distance_CB + b * E_magnitude);
            // auto c2 = (c - b * )    // 需要梯度计算，挖坑
            auto cellGradient = grad(*this, gradientMethod_);
            // Tp = c2 = ()


            // faceField.setValue(boundaryFaceIndex, bc.)
        }
    }
}

template<typename Tp>
inline void Field<Tp>::setGradientMethod(GradientMethod method)
{
    gradientMethod_ = method;
}
