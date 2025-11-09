#pragma once

#include "FaceField.hpp"
#include "CellField.hpp"


template <typename T>
class Field
{
    using ULL = unsigned long long;
public:
    Field() = delete;
    Field(const std::string& name, Mesh* mesh);
    Field(const Field<T>&) = delete;
    Field(Field<T>&&) noexcept = default;
    Field<T>& operator=(const Field<T>&) = delete;
    Field<T>& operator=(Field<T>&&) noexcept = default;
    ~Field() = default;

public:
    /* =========赋值======== */ // 与BaseField，接口一样
    // 统一赋值
    void setValue(const T& value);
    // 采用函数对象
    void setValue(const std::function<T(Scalar, Scalar, Scalar)>& func);


    // 获取器
    const FaceField<T>& faceField() const { return faceField_; }
    const CellField<T>& cellField() const { return cellField_; }
    std::string name() const { return name_; }

    // 场是否有效
    bool isValid() const;

    /* --------设置边界条件-------- */
    // a * φ + b * ∂φ/∂n = c
    void setBoundaryCondition(const std::string& name, Scalar a, Scalar b, const T& c);



private:
    FaceField<T> faceField_;
    CellField<T> cellField_;
    std::string name_;
};

template<typename T>
inline Field<T>::Field(const std::string& name, Mesh* mesh)
    : faceField_(name, mesh)
    , cellField_(name, mesh)
    , name_(name)
{}

template<typename T>
inline void Field<T>::setValue(const T& value)
{
    faceField_.setValue(value);
    cellField_.setValue(value);
}

template<typename T>
inline void Field<T>::setValue(const std::function<T(Scalar, Scalar, Scalar)>& func)
{
    faceField_.setValue(func);
    cellField_.setValue(func);
}


template<typename T>
inline bool Field<T>::isValid() const
{
    return faceField_.isValid() && cellField_.isValid();
}
