#ifndef FACEFIELD_H_
#define FACEFIELD_H_






#include "BaseField.hpp"
#include "BoundaryCondition.hpp"


template <typename Tp>
class FaceField : public BaseField<Tp>
{
    using ULL = unsigned long long;
public:
    FaceField() = delete;
    // 构造但不初始化，场无效
    FaceField(const std::string& name, Mesh* mesh);
    // 构造初始化，有效
    FaceField(const std::string& name, Mesh* mesh, const Tp& initialValue);

public:
    // 给场赋值
    // void setValue(const T& value);




private:
    /* ========私有接口======== */
    // 获取BoundaryPatch
    const std::unordered_map<std::string, BoundaryPatch>&
    getBoundaryPatch() const;

private:


};

template<typename Tp>
inline FaceField<Tp>::FaceField(const std::string& name, Mesh* mesh)
    : BaseField<Tp>(name, mesh)
{
    this->type_ = field::FieldType::FACE_FIELD;
}

template<typename Tp>
inline FaceField<Tp>::FaceField(const std::string& name, Mesh* mesh, const Tp& initialValue)
    : BaseField<Tp>(name, mesh)
{
    this->type_ = field::FieldType::FACE_FIELD;
    
    ULL faceNum = this->mesh_->getFaceNumber();
    this->data_.resize(faceNum, initialValue);
    this->isValid_= true;
}

template<typename Tp>
inline const std::unordered_map<std::string, BoundaryPatch>& FaceField<Tp>::getBoundaryPatch() const
{
    return this->mesh_->getBoundaryPatches();
}

// template<typename T>
// inline void FaceField<T>::setValue(const T& value)
// {
//     if (this->isValid_)
//     {
//         this->data_.assign(this->data_.size(), value);
//         return;
//     }
//     this->data_.resize(this->mesh_->getFaceNumber(), value);
//     this->isValid_ = true;
// }



#endif // FACEFIELD_H_
