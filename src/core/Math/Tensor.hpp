#ifndef TENSOR_H_
#define TENSOR_H_


#include <iostream>
#include "Vector.hpp"
#include <cmath>


template <typename Tp>
class Tensor
{
    using Scalar = double;
    static constexpr Scalar EPSILON = 1e-12;
public:
    Tensor();
    Tensor(     // 3D
        Tp xx, Tp xy, Tp xz,
        Tp yx, Tp yy, Tp yz,
        Tp zx, Tp zy, Tp zz
    );
    Tensor(     // 2D
        Tp xx, Tp xy,
        Tp yx, Tp yy
    );
    explicit Tensor(Tp diag);                    // 对角初始化
    Tensor(const Tensor<Tp>& src) = default;
    Tensor(Tensor<Tp>&&) noexcept = default;
    Tensor<Tp>& operator=(const Tensor<Tp>& src) = default;
    Tensor<Tp>& operator=(Tensor<Tp>&&) noexcept = default;

    ~Tensor() = default;
public:
    /* ----------基本运算----------- */
    // 张量与张量
    Tensor<Tp> operator+(const Tensor<Tp>& rhs) const;    // Tensor + Tensor
    Tensor<Tp> operator-(const Tensor<Tp>& rhs) const;    // Tensor - Tensor
    Tensor<Tp> operator*(const Tensor<Tp>& rhs) const;    // Tensor * Tensor 相乘
    Tensor<Tp> operator-() const;                        // 取反

    Tp operator&&(const Tensor<Tp>& rhs) const;           // 双点积 T :: T

    // 张量与矢量
    template<typename U>
    Vector<Scalar> operator*(const Vector<U>& rhs) const;   // Tensor * Vector

    // 张量与标量
    Tensor<Tp> operator*(const Scalar rhs) const;           // Tensor * Scalar
    template<typename U>
    friend Tensor<Scalar> operator*(const Scalar lhs, const Tensor<U>& rhs);
    Tensor<Tp> operator/(const Scalar rhs) const;           // Tensor / Scalar

    // += -= /= *=
    Tensor<Tp>& operator+=(const Tensor<Tp>& rhs);        // Tensor += Tensor
    Tensor<Tp>& operator-=(const Tensor<Tp>& rhs);        // Tensor -= Tensor
    Tensor<Tp>& operator*=(const Scalar rhs);            // Tensor *= Scalar
    Tensor<Tp>& operator/=(const Scalar rhs);            // Tensor /= Scalar


    // 输出流重载
    template<typename U>
    friend std::ostream& operator<<(std::ostream& out, const Tensor<U>& tensor);


    // 是否为0张量
    bool isZeroTensor(const Scalar epsilon = EPSILON) const;


private:
    // 判断标量（除数）是否为0
    static bool isZero(const Scalar value, const Scalar epsilon = EPSILON);



private:
    Tp xx_, xy_, xz_;
    Tp yx_, yy_, yz_;
    Tp zx_, zy_, zz_;
};





#pragma region 函数实现

template<typename Tp>
inline Tensor<Tp>::Tensor()
    : xx_(0), xy_(0), xz_(0)
    , yx_(0), yy_(0), yz_(0)
    , zx_(0), zy_(0), zz_(0)
{}

template<typename Tp>
inline Tensor<Tp>::Tensor(Tp xx, Tp xy, Tp xz, Tp yx, Tp yy, Tp yz, Tp zx, Tp zy, Tp zz)
    : xx_(xx), xy_(xy), xz_(xz)
    , yx_(yx), yy_(yy), yz_(yz)
    , zx_(zx), zy_(zy), zz_(zz)
{}

template<typename Tp>
inline Tensor<Tp>::Tensor(Tp xx, Tp xy, Tp yx, Tp yy)
    : xx_(xx), xy_(xy), xz_(0)
    , yx_(yx), yy_(yy), yz_(0)
    , zx_(0), zy_(0), zz_(0)
{}

template<typename Tp>
inline Tensor<Tp>::Tensor(Tp diag)
    : xx_(diag), xy_(0), xz_(0)
    , yx_(0), yy_(diag), yz_(0)
    , zx_(0), zy_(0), zz_(diag)
{}

template<typename Tp>
inline Tensor<Tp> Tensor<Tp>::operator+(const Tensor<Tp>& rhs) const
{
    return Tensor<Tp>(
        xx_ + rhs.xx_, xy_ + rhs.xy_, xz_ + rhs.xz_,
        yx_ + rhs.yx_, yy_ + rhs.yy_, yz_ + rhs.yz_,
        zx_ + rhs.zx_, zy_ + rhs.zy_, zz_ + rhs.zz_
    );
}

template<typename Tp>
inline Tensor<Tp> Tensor<Tp>::operator-(const Tensor<Tp>& rhs) const
{
    return Tensor<Tp>(
        xx_ - rhs.xx_, xy_ - rhs.xy_, xz_ - rhs.xz_,
        yx_ - rhs.yx_, yy_ - rhs.yy_, yz_ - rhs.yz_,
        zx_ - rhs.zx_, zy_ - rhs.zy_, zz_ - rhs.zz_
    );
}

template<typename Tp>
inline Tensor<Tp> Tensor<Tp>::operator*(const Tensor<Tp>& rhs) const
{
    return Tensor<Tp>(
        xx_ * rhs.xx_ + xy_ * rhs.yx_ + xz_ * rhs.zx_,
        xx_ * rhs.xy_ + xy_ * rhs.yy_ + xz_ * rhs.zy_,
        xx_ * rhs.xz_ + xy_ * rhs.yz_ + xz_ * rhs.zz_,

        yx_ * rhs.xx_ + yy_ * rhs.yx_ + yz_ * rhs.zx_,
        yx_ * rhs.xy_ + yy_ * rhs.yy_ + yz_ * rhs.zy_,
        yx_ * rhs.xz_ + yy_ * rhs.yz_ + yz_ * rhs.zz_,

        zx_ * rhs.xx_ + zy_ * rhs.yx_ + zz_ * rhs.zx_,
        zx_ * rhs.xy_ + zy_ * rhs.yy_ + zz_ * rhs.zy_,
        zx_ * rhs.xz_ + zy_ * rhs.yz_ + zz_ * rhs.zz_
    );
}

template<typename Tp>
inline Tensor<Tp> Tensor<Tp>::operator-() const
{
    return Tensor<Tp>(
        -xx_, -xy_, -xz_,
        -yx_, -yy_, -yz_,
        -zx_, -zy_, -zz_
    );
}

template<typename Tp>
inline Tp Tensor<Tp>::operator&&(const Tensor<Tp>& rhs) const
{
    return (
        xx_ * rhs.xx_ + xy_ * rhs.yx_ + xz_ * rhs.zx_ +
        yx_ * rhs.xy_ + yy_ * rhs.yy_ + yz_ * rhs.zy_ +
        zx_ * rhs.xz_ + zy_ * rhs.yz_ + zz_ * rhs.zz_
        );
}

template<typename Tp>
template<typename U>
inline Vector<typename Tensor<Tp>::Scalar> Tensor<Tp>::operator*(const Vector<U>& rhs) const
{
    return Vector<Scalar>(
        xx_ * rhs.x() + xy_ * rhs.y() + xz_ * rhs.z(),
        yx_ * rhs.x() + yy_ * rhs.y() + yz_ * rhs.z(),
        zx_ * rhs.x() + zy_ * rhs.y() + zz_ * rhs.z()
    );
}

template<typename Tp>
inline Tensor<Tp> Tensor<Tp>::operator*(const Scalar rhs) const
{
    return Tensor<Tp>(
        xx_ * rhs, xy_ * rhs, xz_ * rhs,
        yx_ * rhs, yy_ * rhs, yz_ * rhs,
        zx_ * rhs, zy_ * rhs, zz_ * rhs
    );
}

template<typename Tp>
inline Tensor<Tp> Tensor<Tp>::operator/(const Scalar rhs) const
{
    if (isZero(rhs))
    {
        std::cerr << "Error: Division by zero in Tensor::operator/." << std::endl;
        throw std::invalid_argument("Division by zero");
    }

    return Tensor<Tp>(
        xx_ / rhs, xy_ / rhs, xz_ / rhs,
        yx_ / rhs, yy_ / rhs, yz_ / rhs,
        zx_ / rhs, zy_ / rhs, zz_ / rhs
    );
}

template<typename Tp>
inline Tensor<Tp>& Tensor<Tp>::operator+=(const Tensor<Tp>& rhs)
{
    xx_ += rhs.xx_; xy_ += rhs.xy_; xz_ += rhs.xz_;
    yx_ += rhs.yx_; yy_ += rhs.yy_; yz_ += rhs.yz_;
    zx_ += rhs.zx_; zy_ += rhs.zy_; zz_ += rhs.zz_;
    return *this;
}

template<typename Tp>
inline Tensor<Tp>& Tensor<Tp>::operator-=(const Tensor<Tp>& rhs)
{
    xx_ -= rhs.xx_; xy_ -= rhs.xy_; xz_ -= rhs.xz_;
    yx_ -= rhs.yx_; yy_ -= rhs.yy_; yz_ -= rhs.yz_;
    zx_ -= rhs.zx_; zy_ -= rhs.zy_; zz_ -= rhs.zz_;
    return *this;
}

template<typename Tp>
inline Tensor<Tp>& Tensor<Tp>::operator*=(const Scalar rhs)
{
    xx_ *= rhs; xy_ *= rhs; xz_ *= rhs;
    yx_ *= rhs; yy_ *= rhs; yz_ *= rhs;
    zx_ *= rhs; zy_ *= rhs; zz_ *= rhs;
    return *this;
}

template<typename Tp>
inline Tensor<Tp>& Tensor<Tp>::operator/=(const Scalar rhs)
{
    if (isZero(rhs))
    {
        std::cerr << "Error: Division by zero in Tensor::operator/=." << std::endl;
        throw std::invalid_argument("Division by zero");
    }

    xx_ /= rhs; xy_ /= rhs; xz_ /= rhs;
    yx_ /= rhs; yy_ /= rhs; yz_ /= rhs;
    zx_ /= rhs; zy_ /= rhs; zz_ /= rhs;
    return *this;
}

template<typename Tp>
inline bool Tensor<Tp>::isZeroTensor(const Scalar epsilon) const
{
    return (
        isZero(xx_, epsilon) && isZero(xy_, epsilon) &&
        isZero(xz_, epsilon) && isZero(yx_, epsilon) &&
        isZero(yy_, epsilon) && isZero(yz_, epsilon) &&
        isZero(zx_, epsilon) && isZero(zy_, epsilon) &&
        isZero(zz_, epsilon)
        );
}

template<typename Tp>
inline bool Tensor<Tp>::isZero(const Scalar value, const Scalar epsilon)
{
    return (std::abs(value) < epsilon);
}



template<typename U1, typename U2, typename Scalar = double>
inline Tensor<Scalar> operator*(const Vector<U1>& lhs, const Vector<U2>& rhs)       // Vector * Vector  并矢
{
    // using Scalar = typename Vector<U1>::Scalar;
    return Tensor<Scalar>(
        lhs.x() * rhs.x(), lhs.x() * rhs.y(), lhs.x() * rhs.z(),
        lhs.y() * rhs.x(), lhs.y() * rhs.y(), lhs.y() * rhs.z(),
        lhs.z() * rhs.x(), lhs.z() * rhs.y(), lhs.z() * rhs.z()
    );
}

template<typename U>
inline Tensor<Scalar> operator*(const Scalar lhs, const Tensor<U>& rhs)
{
    return rhs * lhs;
}

template<typename U>
inline std::ostream& operator<<(std::ostream& out, const Tensor<U>& tensor)
{
    out << "\n"
        << "[ " << tensor.xx_ << " , " << tensor.xy_ << " , " << tensor.xz_ << " ]\n"
        << "[ " << tensor.yx_ << " , " << tensor.yy_ << " , " << tensor.yz_ << " ]\n"
        << "[ " << tensor.zx_ << " , " << tensor.zy_ << " , " << tensor.zz_ << " ]";
    return out;
}



#endif // TENSOR_H_
