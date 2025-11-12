#pragma once

#include "FaceField.hpp"
#include "CellField.hpp"
#include "FieldOperators/Gradient.hpp"
#include <algorithm>


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
    std::string getName() const;
    Mesh* getMesh() const;

    // 场是否有效
    bool isValid() const;
    // 边界条件是否有效
    bool isBoundaryConditionValid() const;




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
    CellField<decltype(Tp()* Vector<Scalar>())> cellGradientField_0_; // 上一步单元场梯度
    // CellField<typename GradientType<Tp>::Type> cellGradientField_; // 单元场梯度
    // CellField<decltype(Tp()* Vector<Scalar>())> cellGradientField_; // 单元场梯度
    std::string name_;
    GradientMethod gradientMethod_{ GradientMethod::GAUSS_GREEN };
    std::unordered_map<std::string, BoundaryCondition<Tp>> boundaryConditions_;
    Interpolation<Tp> interpolator_;     // 插值函数对象
};

#pragma region 函数实现

template<typename Tp>
inline Field<Tp>::Field(const std::string& name, Mesh* mesh)
    : faceField_(name, mesh)
    , cellField_(name, mesh)
    , cellField_0_(name + "_0", mesh)
    , cellGradientField_0_(name + "_grad_0", mesh)
    , name_(name)
    , interpolator_(Interpolation<Tp>())
{
    const std::unordered_map<std::string, BoundaryPatch>& boundaryPatches = mesh->getBoundaryPatches();
    for (const auto& [name, patch] : boundaryPatches)
    {
        boundaryConditions_.emplace(name, BoundaryCondition<Tp>(patch));
    }
}

template<typename Tp>
inline void Field<Tp>::setValue(const Tp& value)
{
    faceField_.setValue(value);
    cellField_.setValue(value);

    // 利用面值计算单元中心梯度
    cellGradientField_0_ = grad(*this, gradientMethod_);
}

template<typename Tp>
inline void Field<Tp>::setValue(const std::function<Tp(Scalar, Scalar, Scalar)>& func)
{
    faceField_.setValue(func);
    cellField_.setValue(func);

    // 利用面值计算单元中心梯度
    cellGradientField_0_ = grad(*this, gradientMethod_);
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
inline std::string Field<Tp>::getName() const
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
inline bool Field<Tp>::isBoundaryConditionValid() const
{
    return std::all_of(
        boundaryConditions_.begin(),
        boundaryConditions_.end(),
        [](const std::pair<std::string, BoundaryCondition<Tp>>& bc) {
            return bc.second.isValid();
        }
    );
}

template<typename Tp>
inline void Field<Tp>::cellToFace(interpolation::Scheme scheme)
{
    if (!isValid())
    {
        std::cerr << "Field is not valid!" << std::endl;
        throw std::runtime_error("Field is not valid!");
    }
    if (!isBoundaryConditionValid())
    {
        std::cerr << "Boundary condition is not valid!" << std::endl;
        throw std::runtime_error("Boundary condition is not valid!");
    }

    using LL = long long;
    Mesh* mesh = this->getMesh();      // 获取网格
    const std::vector<Face>& faces = mesh->getFaces();  // 面列表
    const std::vector<Cell>& cells = mesh->getCells();

    // 先遍历内部面
    std::vector<ULL> internalFaceIndexes = mesh->getInternalFaceIndexes();
    const CellField<Tp>& cellField = this->getCellField();

    for (const ULL internalFaceIndex : internalFaceIndexes)     // 遍历内部面
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
            const decltype(Tp() * Vector<Scalar>())& cellGradient = cellField_0_[ownerIndex];
            Tp c2 = (c - b * cellGradient[ownerIndex] * T) * distance_CB / (a * distance_CB + b * E_magnitude);

            // 计算边界面值
            Tp boundaryFaceValue = c1 * cellField[ownerIndex] + c2;

            faceField_.setValue(boundaryFaceIndex, boundaryFaceValue);
        }
    }
    // 利用新的面值记录计算本时间步的梯度，用于下一时间步的边界面值计算
    cellField_0_ = grad(*this);
}

template<typename Tp>
inline void Field<Tp>::setGradientMethod(GradientMethod method)
{
    gradientMethod_ = method;
}

template<typename Tp>
inline void Field<Tp>::setBoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c)
{
    if (!isValid())
    {
        throw std::runtime_error("Field is not valid!");
    }

    // 判断是否存在该边界
    auto it = boundaryConditions_.find(name);
    if (it == boundaryConditions_.end())
    {
        std::cerr << "Boundary condition " << name << " does not exist!" << std::endl;
        throw std::runtime_error("Boundary condition does not exist!");
    }
    it->second.setBoundaryCondition(name, a, b, c);
    it->second.setValid();      // 设置该边界条件有效
}
