#pragma once

#include <string>


// a * φ + b * ∂φ/∂n = c

template<typename Tp>
class BoundaryCondition
{
    using ULL = unsigned long long;
    using Scalar = double;
public:
    BoundaryCondition() = delete;
    BoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c);
    BoundaryCondition(const BoundaryCondition<Tp>&) = delete;
    BoundaryCondition(BoundaryCondition&&) = delete;
    BoundaryCondition& operator=(const BoundaryCondition<Tp>&) = delete;
    BoundaryCondition& operator=(BoundaryCondition<Tp>&&) = delete;

public:
    std::string getName() const;
    ULL getStartFace() const;
    ULL getNFace() const;

    Scalar get_a() const;
    Scalar get_b() const;
    const Tp& get_c() const;

private:
    std::string name_;      // 边界名
    Scalar a_;
    Scalar b_;      
    Tp c_;                   // 标量/向量/张量

    ULL nFace_;         // 边界面数量   
    ULL startFace_;     // 边界起始面编号

};

#pragma region 函数实现

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

