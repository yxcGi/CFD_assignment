#ifndef BASEFIELD_H_
#define BASEFIELD_H_




#include "Mesh.h"
#include <string>
#include <vector>
#include <functional>
#include "FieldType.hpp"
#include "BoundaryCondition.hpp"





template <typename Tp>
class BaseField             // 基类
{
    using ULL = unsigned long long;
public:
    BaseField() = delete;
    BaseField(const std::string& name, Mesh* mesh);
    BaseField(const BaseField<Tp>& src) = delete;
    BaseField(BaseField<Tp>&& src) noexcept = default;
    BaseField<Tp>& operator=(const BaseField<Tp>& src) = delete;
    BaseField<Tp>& operator=(BaseField<Tp>&& src) noexcept = default;
    ~BaseField() = default;

public:

    // 获取器
    std::string getName() const;
    ULL getSize() const;
    field::FieldType getType() const;
    Mesh* getMesh() const;
    const std::unordered_map<std::string, BoundaryCondition<Tp>>& getBoundaryConditions() const;

    /* ---------赋值--------- */
    // 赋统一值
    void setValue(const Tp& value);
    // 对特定场点赋值
    void setValue(ULL id, const Tp& value);
    // 采用函数对象，lambda表达式赋值(传入的是坐标值)
    void setValue(const std::function<Tp(Scalar, Scalar, Scalar)>& func);

    // 场是否有效，必须要setValue初始化后才能有效
    bool isValid() const;


    /* --------设置边界条件-------- */
    // a * φ + b * ∂φ/∂n = c
    void setBoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c);



    // 访问
    Tp& operator[](ULL i);
    const Tp& operator[](ULL i) const;
protected:
    ULL getDataNumer() const;   // 直接获取当前特定场的数据数量


protected:
    std::string name_;          // 场名称
    Mesh* mesh_;                // 网格指针
    std::vector<Tp> data_;       // 场数据
    field::FieldType type_;            // 场类型

    // 存储边界条件
    std::unordered_map<std::string, BoundaryCondition<Tp>> boundaryConditions_;

    bool isValid_;              // 场是否有效
};

#pragma region 函数实现



template<typename Tp>
inline BaseField<Tp>::BaseField(const std::string& name, Mesh* mesh)
    : name_(name)
    , mesh_(mesh)
    , type_(field::FieldType::BASE)
    , isValid_(false)
{}

template<typename Tp>
inline std::string BaseField<Tp>::getName() const
{
    return name_;
}

template<typename Tp>
inline typename BaseField<Tp>::ULL BaseField<Tp>::getSize() const
{
    if (isValid_)
    {
        return data_.size();
    }
    std::cerr << "Error: Field is not valid!" << std::endl;
    return 0;
}

template<typename Tp>
inline field::FieldType BaseField<Tp>::getType() const
{
    return type_;
}

template<typename Tp>
inline Mesh* BaseField<Tp>::getMesh() const
{
    return mesh_;
}

template<typename Tp>
inline const std::unordered_map<std::string, BoundaryCondition<Tp>>& BaseField<Tp>::getBoundaryConditions() const
{
    return boundaryConditions_;
}

template<typename Tp>
inline void BaseField<Tp>::setValue(const Tp& value)
{
    if (isValid_)
    {
        data_.assign(this->data_.size(), value);
        return;
    }
    data_.resize(getDataNumer(), value);
    isValid_ = true;
}

template<typename Tp>
inline void BaseField<Tp>::setValue(ULL id, const Tp& value)
{
    if (isValid_)
    {
        data_[id] = value;
        return;
    }
    else
    {
        data_.resize(getDataNumer(), Tp());
        data_[id] = value;
        isValid_ = true;
    }
}

template<typename Tp>
inline void BaseField<Tp>::setValue(const std::function<Tp(Scalar, Scalar, Scalar)>& func)
{
    if (!isValid_)
    {
        std::cerr << "Error: Field is not valid!" << std::endl;
        throw std::runtime_error("Field is not valid!");
    }


    if (type_ == field::FieldType::CELL_FIELD)
    {
        for (ULL i = 0; i < data_.size(); ++i)
        {
            const Point& cellCenter = mesh_->getCells()[i].getCenter();
            data_[i] = func(cellCenter.x(), cellCenter.y(), cellCenter.z());
        }
    }
    else if (type_ == field::FieldType::FACE_FIELD)
    {
        for (ULL i = 0; i < data_.size(); ++i)
        {
            const Point& faceCenter = mesh_->getFaces()[i].getCenter();

            data_[i] = func(faceCenter.x(), faceCenter.y(), faceCenter.z());
        }
    }
    else if (type_ == field::FieldType::NODE_FIELD)
    {
        for (ULL i = 0; i < data_.size(); ++i)
        {
            const Point& node = mesh_->getPoints()[i];
            data_[i] = func(node.x(), node.y(), node.z());
        }
    }
    else if (type_ == field::FieldType::BASE)
    {
        std::cerr << "Error: Field type is not set!" << std::endl;
        throw std::runtime_error("Field type is not set!");
    }
}

template<typename Tp>
inline bool BaseField<Tp>::isValid() const
{
    return isValid_;
}

template<typename Tp>
inline void BaseField<Tp>::setBoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c)
{
    bool isExistBoundaryPatch = false;      // 是否存在边界
    const std::unordered_map<std::string, BoundaryPatch>& boundaryPatches = this->mesh_->getBoundaryPatches();
    for (const auto& [patchName, patch] : boundaryPatches)
    {
        if (name == patchName)
        {
            isExistBoundaryPatch = true;
            break;
        }
    }

    if (!isExistBoundaryPatch)      // 不存在边界，抛出异常
    {
        std::cerr << "Error: No such boundary patch named " << name << std::endl;
        throw std::runtime_error("No such boundary patch named " + name);
    }

    // 存在边界，设置边界条件
    boundaryConditions_.emplace(name, BoundaryCondition<Tp>(name, a, b, c));
}

template<typename Tp>
inline Tp& BaseField<Tp>::operator[](ULL i)
{
    if (isValid_)
    {
        return data_[i];
    }
    std::cerr << "Error: Field is not valid!" << std::endl;
    throw std::runtime_error("Field is not valid!");
}

template<typename Tp>
inline const Tp& BaseField<Tp>::operator[](ULL i) const
{
    if (isValid_)
    {
        return data_[i];
    }
    std::cerr << "Error: Field is not valid!" << std::endl;
    throw std::runtime_error("Field is not valid!");
}

template<typename Tp>
inline typename BaseField<Tp>::ULL BaseField<Tp>::getDataNumer() const
{
    if (type_ == field::FieldType::CELL_FIELD)
    {
        return mesh_->getCellNumber();
    }
    else if (type_ == field::FieldType::FACE_FIELD)
    {
        return mesh_->getFaceNumber();
    }
    else if (type_ == field::FieldType::NODE_FIELD)
    {
        return mesh_->getPointNumber();
    }
    else if (type_ == field::FieldType::BASE)
    {
        std::cerr << "Error: BaseField is not allowed to be used directly!" << std::endl;
        throw std::runtime_error("BaseField is not allowed to be used directly!");
    }
    std::cerr << "Unknown field type!" << std::endl;
    throw std::runtime_error("Unknown field type!");
}



#endif // BASEFIELD_H_
