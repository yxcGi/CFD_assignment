#pragma once


#include "BaseField.hpp"


template <typename T>
class FaceField : public BaseField<T>
{
    using ULL = unsigned long long;
public:
    FaceField() = delete;
    // 构造但不初始化，场无效
    FaceField(const std::string& name, Mesh* mesh);
    // 构造初始化，有效
    FaceField(const std::string& name, Mesh* mesh, const T& initialValue);

public:
    // 给场赋值
    void setValue(const T& value);


};

template<typename T>
inline FaceField<T>::FaceField(const std::string& name, Mesh* mesh)
    : BaseField<T>(name, mesh)
{
    this->type_ = field::FieldType::FACE_FIELD;
}

template<typename T>
inline FaceField<T>::FaceField(const std::string& name, Mesh* mesh, const T& initialValue)
    : BaseField<T>(name, mesh)
{
    this->type_ = field::FieldType::FACE_FIELD;
    
    ULL faceNum = this->mesh_->getFaceNumber();
    this->data_.resize(faceNum, initialValue);
    this->isValid_= true;
}

template<typename T>
inline void FaceField<T>::setValue(const T& value)
{
    if (this->isValid_)
    {
        this->data_.assign(this->data_.size(), value);
        return;
    }
    this->data_.resize(this->mesh_->getFaceNumber(), value);
    this->isValid_ = true;
}
