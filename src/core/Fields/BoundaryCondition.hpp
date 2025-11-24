#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_





#include <string>
#include "Mesh.h"
#include "BoundaryPatch.h"
// a * φ + b * ∂φ/∂n = c

template<typename Tp>
class BoundaryCondition
{
    using BType = BoundaryPatch::BoundaryType;
    using ULL = unsigned long long;
    using Scalar = double;
public:
    BoundaryCondition() = delete;
    // BoundaryCondition(Mesh* mesh); //
    BoundaryCondition(const BoundaryPatch& patch);
    BoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c);
    BoundaryCondition(const BoundaryCondition<Tp>&) = default;
    BoundaryCondition(BoundaryCondition&&) = default;
    BoundaryCondition& operator=(const BoundaryCondition<Tp>&) = delete;
    BoundaryCondition& operator=(BoundaryCondition<Tp>&&) = delete;

public:
    std::string getName() const;
    ULL getStartFace() const;
    ULL getNFace() const;
    BType getType() const;

    Scalar get_a() const;
    Scalar get_b() const;
    const Tp& get_c() const;
    void setValid(bool flag = true);

    // set
    void setBoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c);



    
    // 已经设置
    bool isValid() const;

private:
    std::string name_;      // 边界名
    Scalar a_;
    Scalar b_;
    Tp c_;                   // 标量/向量/张量

    ULL nFace_;         // 边界面数量   
    ULL startFace_;     // 边界起始面编号
    BType type_;
    bool isSet_ = false;
};

#pragma region 函数实现

// template<typename Tp>
// inline BoundaryCondition<Tp>::BoundaryCondition(Mesh* mesh)
// {
//     using BType = BoundaryPatch::BoundaryType;
//     const std::unordered_map<std::string, BoundaryPatch>& boundaryPatches = mesh->getBoundaryPatches();
//     for (const auto& [name, patch] : boundaryPatches)
//     {
        
//     }
// }

template<typename Tp>
inline BoundaryCondition<Tp>::BoundaryCondition(const BoundaryPatch& patch)
    : name_(patch.getName())
    , nFace_(patch.getNFace())
    , startFace_(patch.getStartFace())
    , type_(patch.getType())
{}
template<typename Tp>
inline BoundaryCondition<Tp>::BoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c)
    : name_(name)
    , a_(a)
    , b_(b)
    , c_(c)
{}

template<typename Tp>
inline std::string BoundaryCondition<Tp>::getName() const
{
    return name_;
}

template<typename Tp>
inline typename BoundaryCondition<Tp>::ULL BoundaryCondition<Tp>::getStartFace() const
{
    return startFace_;
}

template<typename Tp>
inline typename BoundaryCondition<Tp>::ULL BoundaryCondition<Tp>::getNFace() const
{
    return nFace_;
}

template<typename Tp>
inline typename BoundaryCondition<Tp>::BType BoundaryCondition<Tp>::getType() const
{
    return type_;
}


template<typename Tp>
inline typename BoundaryCondition<Tp>::Scalar BoundaryCondition<Tp>::get_a() const
{
    return a_;
}

template<typename Tp>
inline typename BoundaryCondition<Tp>::Scalar BoundaryCondition<Tp>::get_b() const
{
    return b_;
}

template<typename Tp>
inline const Tp& BoundaryCondition<Tp>::get_c() const
{
    return c_;
}

template<typename Tp>
inline void BoundaryCondition<Tp>::setValid(bool flag)
{
    isSet_ = flag;
}

template<typename Tp>
inline void BoundaryCondition<Tp>::setBoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c)
{
    name_ = name;
    a_ = a;
    b_ = b;
    c_ = c;
}

template<typename Tp>
inline bool BoundaryCondition<Tp>::isValid() const
{
    return isSet_;
}





#endif // BOUNDARYCONDITION_H_
