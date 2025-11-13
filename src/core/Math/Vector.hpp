#ifndef VECTOR_H_
#define VECTOR_H_


#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>


using Scalar = double;

template <typename Tp>
class Vector;

using Point = Vector<double>;


// 前置声明
template<typename Tp>
class Tensor;


/**
 * @brief 向量/点的定义，模板类
 * @tparam T 数据类型
 */
template <typename Tp>
class Vector
{
    using Scalar = double;
    // using ScalarVector = Vector<Scalar>;
    static constexpr Scalar EPSILON = 1e-12;
public:
#pragma region 构造析构
    // 构造函数
    Vector();
    Vector(Tp x, Tp y, Tp z = 0);  // 3D/2D
    Vector(const Vector<Tp>& v);
    Vector(Vector<Tp>&&) noexcept = default;

    // 赋值重载
    Vector<Tp>& operator=(const Vector<Tp>& src);
    Vector<Tp>& operator=(Vector<Tp>&& src) noexcept = default;

    // 析构
    ~Vector() = default;
#pragma endregion

public:
#pragma region 运算符重载
    // 基本运算
    template<typename U>
    Vector<Scalar> operator+(const Vector<U>& rhs) const;        // Vector + Vector
    template<typename U>
    Vector<Scalar> operator-(const Vector<U>& rhs) const;        // Vector - Vector
    Vector<Tp> operator-() const;                            // -Vector  取反

    // 标量乘法
    Vector<Scalar> operator*(const Scalar rhs) const;                          // Vector * Scalar 右乘
    template<typename U>
    friend Vector<typename Vector<U>::Scalar> operator*(const typename Vector<U>::Scalar lhs, const Vector<U>& rhs);   // Scalar * Vector 左乘
    Vector<Scalar> operator/(const Scalar rhs) const;                          // Vector / Scalar

    // 复合复制运算
    Vector<Tp>& operator+=(const Vector<Tp>& rhs);        // Vector += Vector
    Vector<Tp>& operator-=(const Vector<Tp>& rhs);        // Vector -= Vector
    Vector<Tp>& operator*=(const Scalar rhs);        // Vector *= Scalar
    Vector<Tp>& operator/=(const Scalar rhs);        // Vector /= Scalar

    // 点乘
    Scalar operator&(const Vector<Tp>& rhs) const;          // Vector & Vector  点乘

    // 叉乘
    Vector<Tp> operator^(const Vector<Tp>& rhs)const;          // Vector ^ Vector  叉乘

    // 矢量与张量运算
    template<typename U>
    Vector<Scalar> operator*(const Tensor<U>& rhs) const;   // Vector * Tensor

    // 并矢：房子啊Tensor里为，自由函数
#pragma endregion

#pragma region 常用接口
    // 取模
    Scalar magnitude() const;
    // 模长平方
    Scalar magnitudeSquared() const;
    // 单位向量
    Vector<Scalar> unitVector() const;
    // 两点距离
    Scalar getDistance(const Vector<Tp>& other) const;
    // 取坐标
    Tp& x() { return x_; }
    Tp& y() { return y_; }
    Tp& z() { return z_; }
    const Tp& x() const { return x_; }
    const Tp& y() const { return y_; }
    const Tp& z() const { return z_; }

    // 判断是否为0向量
    bool isZeroVector(const Scalar epsilon = EPSILON) const;

    // 输出流
    template<typename U>
    friend std::ostream& operator<<(std::ostream& out, const Vector<U>& vec);
#pragma endregion

private:
    // 判断标量（除数）是否为0
    static bool isZero(const Scalar value, const Scalar epsilon = EPSILON); // 是否为0向量




private:
    Tp x_, y_, z_;   // 坐标
};







#pragma region 函数实现
template<typename Tp>
inline Vector<Tp>::Vector() : x_(0.0), y_(0.0), z_(0.0)
{}

template<typename Tp>
inline Vector<Tp>::Vector(Tp x, Tp y, Tp z) : x_(x), y_(y), z_(z)
{}

template <typename Tp>
inline Vector<Tp>::Vector(const Vector<Tp>& src)
{
    x_ = src.x_;
    y_ = src.y_;
    z_ = src.z_;
}


template<typename Tp>
inline Vector<Tp>& Vector<Tp>::operator=(const Vector<Tp>& src)
{
    if (&src == this)
    {
        return *this;
    }

    x_ = src.x_;
    y_ = src.y_;
    z_ = src.z_;
    return *this;
}


template<typename Tp>
template<typename U>
inline Vector<typename Vector<Tp>::Scalar> Vector<Tp>::operator+(const Vector<U>& rhs) const
{
    return Vector<Scalar>(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_);
}


template<typename Tp>
template<typename U>
inline Vector<typename Vector<Tp>::Scalar> Vector<Tp>::operator-(const Vector<U>& rhs) const
{
    return Vector<Scalar>(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_);
}

template<typename Tp>
template<typename U>
inline Vector<typename Vector<Tp>::Scalar> Vector<Tp>::operator*(const Tensor<U>& rhs) const
{
    return Vector<Scalar>(
        x_ * rhs.xx() + y_ * rhs.yx() + z_ * rhs.zx(),
        x_ * rhs.xy() + y_ * rhs.yy() + z_ * rhs.zy(),
        x_ * rhs.xz() + y_ * rhs.yz() + z_ * rhs.zz()
    );
}


template<typename U>
inline Vector<typename Vector<U>::Scalar> operator*(const typename Vector<U>::Scalar lhs, const Vector<U>& rhs)
{
    return rhs * lhs;
}

template<typename U>
inline std::ostream& operator<<(std::ostream& out, const Vector<U>& vec)
{
    // 保留两位小数
    out << std::fixed << std::setprecision(4);
    out << "(" << vec.x_ << ", " << vec.y_ << ", " << vec.z_ << ")";
    return out;
}


template<typename Tp>
inline Vector<Tp> Vector<Tp>::operator-() const
{
    return Vector<Tp>(-x_, -y_, -z_);
}

template<typename Tp>
inline Vector<typename Vector<Tp>::Scalar> Vector<Tp>::operator*(const Scalar rhs) const
{
    return Vector<Scalar>(x_ * rhs, y_ * rhs, z_ * rhs);
}

template<typename Tp>
inline Vector<typename Vector<Tp>::Scalar> Vector<Tp>::operator/(const Scalar rhs) const
{
    if (Vector<Tp>::isZero(rhs))
    {
        std::cerr << "Error: Division by zero in Vector::operator/." << std::endl;
        throw std::invalid_argument("Division by zero");
    }
    return Vector<Tp>(x_ / rhs, y_ / rhs, z_ / rhs);
}

template<typename Tp>
inline Vector<Tp>& Vector<Tp>::operator+=(const Vector<Tp>& rhs)
{
    x_ += rhs.x_;
    y_ += rhs.y_;
    z_ += rhs.z_;
    return *this;
}

template<typename Tp>
inline Vector<Tp>& Vector<Tp>::operator-=(const Vector<Tp>& rhs)
{
    x_ -= rhs.x_;
    y_ -= rhs.y_;
    z_ -= rhs.z_;
    return *this;
}

template<typename Tp>
inline Vector<Tp>& Vector<Tp>::operator*=(const Scalar rhs)
{
    x_ *= rhs;
    y_ *= rhs;
    z_ *= rhs;
    return *this;
}

template<typename Tp>
inline Vector<Tp>& Vector<Tp>::operator/=(const Scalar rhs)
{
    if (Vector<Tp>::isZero(rhs))
    {
        std::cerr << "Error: Division by zero in Vector::operator/=." << std::endl;
        throw std::invalid_argument("Division by zero");
    }
    x_ /= rhs;
    y_ /= rhs;
    z_ /= rhs;
    return *this;
}

template<typename Tp>
inline typename Vector<Tp>::Scalar Vector<Tp>::operator&(const Vector<Tp>& rhs) const
{
    return (x_ * rhs.x_ + y_ * rhs.y_ + z_ * rhs.z_);
}

template<typename Tp>
inline Vector<Tp> Vector<Tp>::operator^(const Vector<Tp>& rhs) const
{
    return Vector<Tp>(
        (y_ * rhs.z_ - z_ * rhs.y_),
        (z_ * rhs.x_ - x_ * rhs.z_),
        (x_ * rhs.y_ - y_ * rhs.x_)
    );
}

template<typename Tp>
inline typename Vector<Tp>::Scalar Vector<Tp>::magnitude() const
{
    return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
}

template<typename Tp>
inline typename Vector<Tp>::Scalar Vector<Tp>::magnitudeSquared() const
{
    return (x_ * x_ + y_ * y_ + z_ * z_);
}

template<typename Tp>
inline Vector<typename Vector<Tp>::Scalar> Vector<Tp>::unitVector() const
{
    if (this->isZeroVector())
    {
        std::cerr << "Error: Cannot compute unit vector of zero vector." << std::endl;
        throw std::invalid_argument("Unit vector of zero vector");
    }
    Scalar mag = this->magnitude();
    return Vector<Scalar>(x_ / mag, y_ / mag, z_ / mag);
}

template<typename Tp>
inline typename Vector<Tp>::Scalar Vector<Tp>::getDistance(const Vector<Tp>& other) const
{
    return std::sqrt(
        (x_ - other.x_) * (x_ - other.x_) +
        (y_ - other.y_) * (y_ - other.y_) +
        (z_ - other.z_) * (z_ - other.z_)
    );
}

template<typename Tp>
inline bool Vector<Tp>::isZero(Scalar value, Scalar epsilon)
{
    return std::abs(value) < epsilon;
}

template<typename Tp>
inline bool Vector<Tp>::isZeroVector(Scalar epsilon) const
{
    return isZero(x_, epsilon) && isZero(y_, epsilon) && isZero(z_, epsilon);
}



#endif // VECTOR_H_
