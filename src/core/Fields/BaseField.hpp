#ifndef BASEFIELD_H_
#define BASEFIELD_H_




#include "Mesh.h"
#include <string>
#include <vector>
#include <functional>
#include "FieldType.hpp"





template <typename T>
class BaseField             // 基类
{
    using ULL = unsigned long long;
public:
    BaseField() = delete;
    BaseField(const std::string& name, Mesh* mesh);
    BaseField(const BaseField<T>& src) = delete;
    BaseField(BaseField<T>&& src) noexcept = default;
    BaseField<T>& operator=(const BaseField<T>& src) = delete;
    BaseField<T>& operator=(BaseField<T>&& src) = default;
    ~BaseField() = default;

public:

    // 获取器
    std::string getName() const;
    ULL getSize() const;
    field::FieldType getType() const;

    /* ---------赋值--------- */
    // 赋统一值
    void setValue(const T& value);
    // 对特定场点赋值
    void setValue(ULL id, const T& value);
    // 采用函数对象，lambda表达式赋值(传入的是坐标值)
    void setValue(std::function<T(Scalar, Scalar, Scalar)> func);



    // 访问
    T& operator[](ULL i);
    const T& operator[](ULL i) const;
protected:
    ULL getDataNumer() const;   // 直接获取当前特定场的数据数量


protected:
    std::string name_;          // 场名称
    Mesh* mesh_;                // 网格指针
    std::vector<T> data_;       // 场数据
    field::FieldType type_;            // 场类型
    bool isValid_;              // 场是否有效
};

template<typename T>
inline BaseField<T>::BaseField(const std::string& name, Mesh* mesh)
    : name_(name)
    , mesh_(mesh)
    , type_(field::FieldType::BASE)
    , isValid_(false)
{}

template<typename T>
inline std::string BaseField<T>::getName() const
{
    return name_;
}

template<typename T>
inline typename BaseField<T>::ULL BaseField<T>::getSize() const
{
    if (isValid_)
    {
        return data_.size();
    }
    std::cerr << "Error: Field is not valid!" << std::endl;
    return 0;
}

template<typename T>
inline field::FieldType BaseField<T>::getType() const
{
    return type_;
}

template<typename T>
inline void BaseField<T>::setValue(const T& value)
{
    if (isValid_)
    {
        data_.assign(this->data_.size(), value);
        return;
    }
    data_.resize(getDataNumer(), value);
    isValid_ = true;
}

template<typename T>
inline void BaseField<T>::setValue(ULL id, const T& value)
{
    if (isValid_)
    {
        data_[id] = value;
        return;
    }
    std::cerr << "Error: Field is not valid!" << std::endl;
    throw std::runtime_error("Field is not valid!");
}

template<typename T>
inline void BaseField<T>::setValue(std::function<T(Scalar, Scalar, Scalar)> func)
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

template<typename T>
inline T& BaseField<T>::operator[](ULL i)
{
    if (isValid_)
    {
        return data_[i];
    }
    std::cerr << "Error: Field is not valid!" << std::endl;
    throw std::runtime_error("Field is not valid!");
}

template<typename T>
inline const T& BaseField<T>::operator[](ULL i) const
{
    if (isValid_)
    {
        return data_[i];
    }
    std::cerr << "Error: Field is not valid!" << std::endl;
    throw std::runtime_error("Field is not valid!");
}

template<typename T>
inline typename BaseField<T>::ULL BaseField<T>::getDataNumer() const
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
