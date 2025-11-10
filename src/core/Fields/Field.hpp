#pragma once

#include "FaceField.hpp"
#include "CellField.hpp"


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
    const FaceField<Tp>& faceField() const { return faceField_; }
    const CellField<Tp>& cellField() const { return cellField_; }
    std::string name() const { return name_; }
    Mesh* mesh() const { return faceField_.getMesh(); }

    // 场是否有效
    bool isValid() const;

    /* --------设置边界条件-------- */
    // a * φ + b * ∂φ/∂n = c
    void setBoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c);



private:
    FaceField<Tp> faceField_;        // 当前面场值
    CellField<Tp> cellField_;        // 当前单元场值
    CellField<Tp> cellField_0_;      // 上一步单元场的值
    std::string name_;
};

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
inline bool Field<Tp>::isValid() const
{
    return faceField_.isValid() && cellField_.isValid();
}
