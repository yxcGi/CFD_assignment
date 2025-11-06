#ifndef VECTOR_H_
#define VECTOR_H_


#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>



template <typename T>
class Vector;

using Point = Vector<double>;


// 前置声明
template<typename T>
class Tensor;

/**
 * @brief 向量/点的定义，模板类
 * @tparam T 数据类型
 */
template <typename T>
class Vector
{
    using Scalar = double;
    // using ScalarVector = Vector<Scalar>;
    static constexpr Scalar EPSILON = 1e-12;
public:
#pragma region 构造析构
    // 构造函数
    Vector();
    Vector(T x, T y, T z = 0);  // 3D/2D
    Vector(const Vector<T>& v);
    Vector(Vector<T>&&) noexcept = default;

    // 赋值重载
    Vector<T>& operator=(const Vector<T>& src);
    Vector<T>& operator=(Vector<T>&& src) noexcept = default;

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
    Vector<T> operator-() const;                            // -Vector  取反

    // 标量乘法
    Vector<Scalar> operator*(const Scalar rhs) const;                          // Vector * Scalar 右乘
    template<typename U>
    friend Vector<Scalar> operator*(const Scalar lhs, const Vector<U>& rhs);   // Scalar * Vector 左乘
    Vector<Scalar> operator/(const Scalar rhs) const;                          // Vector / Scalar

    // 复合复制运算
    Vector<T>& operator+=(const Vector<T>& rhs);        // Vector += Vector
    Vector<T>& operator-=(const Vector<T>& rhs);        // Vector -= Vector
    Vector<T>& operator*=(const Scalar rhs);        // Vector *= Scalar
    Vector<T>& operator/=(const Scalar rhs);        // Vector /= Scalar

    // 点乘
    Scalar operator&(const Vector<T>& rhs);          // Vector & Vector  点乘

    // 叉乘
    Vector<T> operator^(const Vector<T>& rhs);          // Vector ^ Vector  叉乘

    // 矢量与张量运算
    template<typename U>
    Vector<Scalar> operator*(const Tensor<U>& rhs) const;   // Vector * Tensor

    // 并矢
    template<typename U>
    Tensor<Scalar> operator*(const Vector<U>& rhs) const;        // Vector * Vector  并矢
#pragma endregion

#pragma region 常用接口
    // 取模
    Scalar magnitude() const;
    // 模长平方
    Scalar magnitudeSquared() const;
    // 单位向量
    Vector<Scalar> unitVector() const;
    // 两点距离
    Scalar getDistance(const Vector<T>& other) const;
    // 取坐标
    T& x() { return x_; }
    T& y() { return y_; }
    T& z() { return z_; }
    const T& x() const { return x_; }
    const T& y() const { return y_; }
    const T& z() const { return z_; }

    // 判断标量（除数）是否为0
    static bool isZero(Scalar value, Scalar epsilon = EPSILON); // 是否为0向量
    bool isZeroVector(Scalar epsilon = EPSILON) const;

    // 输出流
    template<typename U>
    friend std::ostream& operator<<(std::ostream& out, const Vector<U>& vec);
#pragma endregion






private:
    T x_, y_, z_;   // 坐标
};



#pragma region 函数实现



template<typename T>
inline Vector<T>::Vector() : x_(0.0), y_(0.0), z_(0.0)
{}

template<typename T>
inline Vector<T>::Vector(T x, T y, T z) : x_(x), y_(y), z_(z)
{}

template <typename T>
inline Vector<T>::Vector(const Vector<T>& src)
{
    x_ = src.x_;
    y_ = src.y_;
    z_ = src.z_;
}


template<typename T>
inline Vector<T>& Vector<T>::operator=(const Vector<T>& src)
{
    if (&src == this)
    {
        return *this;
    }

    x_ = src.x_;
    y_ = src.y_;
    z_ = src.z_;
}


template<typename T>
template<typename U>
inline Vector<typename Vector<T>::Scalar> Vector<T>::operator+(const Vector<U>& rhs) const
{
    return Vector<Scalar>(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_);
}


template<typename T>
template<typename U>
inline Vector<typename Vector<T>::Scalar> Vector<T>::operator-(const Vector<U>& rhs) const
{
    return Vector<Scalar>(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_);
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


template<typename T>
inline Vector<T> Vector<T>::operator-() const
{
    return Vector<T>(-x_, -y_, -z_);
}

template<typename T>
inline Vector<typename Vector<T>::Scalar> Vector<T>::operator*(const Scalar rhs) const
{
    return Vector<Scalar>(x_ * rhs, y_ * rhs, z_ * rhs);
}

template<typename T>
inline Vector<typename Vector<T>::Scalar> Vector<T>::operator/(const Scalar rhs) const
{
    if (Vector<T>::isZero(rhs))
    {
        std::cerr << "Error: Division by zero in Vector::operator/." << std::endl;
        throw std::invalid_argument("Division by zero");
    }
    return Vector<T>(x_ / rhs, y_ / rhs, z_ / rhs);
}

template<typename T>
inline Vector<T>& Vector<T>::operator+=(const Vector<T>& rhs)
{
    x_ += rhs.x_;
    y_ += rhs.y_;
    z_ += rhs.z_;
    return *this;
}

template<typename T>
inline Vector<T>& Vector<T>::operator-=(const Vector<T>& rhs)
{
    x_ -= rhs.x_;
    y_ -= rhs.y_;
    z_ -= rhs.z_;
    return *this;
}

template<typename T>
inline Vector<T>& Vector<T>::operator*=(const Scalar rhs)
{
    x_ *= rhs;
    y_ *= rhs;
    z_ *= rhs;
    return *this;
}

template<typename T>
inline Vector<T>& Vector<T>::operator/=(const Scalar rhs)
{
    if (Vector<T>::isZero(rhs))
    {
        std::cerr << "Error: Division by zero in Vector::operator/=." << std::endl;
        throw std::invalid_argument("Division by zero");
    }
    x_ /= rhs;
    y_ /= rhs;
    z_ /= rhs;
    return *this;
}

template<typename T>
inline typename Vector<T>::Scalar Vector<T>::operator&(const Vector<T>& rhs)
{
    return (x_ * rhs.x_ + y_ * rhs.y_ + z_ * rhs.z_);
}

template<typename T>
inline Vector<T> Vector<T>::operator^(const Vector<T>& rhs)
{
    return Vector<T>(
        (y_ * rhs.z_ - z_ * rhs.y_),
        (z_ * rhs.x_ - x_ * rhs.z_),
        (x_ * rhs.y_ - y_ * rhs.x_)
    );
}

template<typename T>
inline typename Vector<T>::Scalar Vector<T>::magnitude() const
{
    return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
}

template<typename T>
inline typename Vector<T>::Scalar Vector<T>::magnitudeSquared() const
{
    return (x_ * x_ + y_ * y_ + z_ * z_);
}

template<typename T>
inline Vector<typename Vector<T>::Scalar> Vector<T>::unitVector() const
{
    if (this->isZeroVector())
    {
        std::cerr << "Error: Cannot compute unit vector of zero vector." << std::endl;
        throw std::invalid_argument("Unit vector of zero vector");
    }
    Scalar mag = this->magnitude();
    return Vector<Scalar>(x_ / mag, y_ / mag, z_ / mag);
}

template<typename T>
inline typename Vector<T>::Scalar Vector<T>::getDistance(const Vector<T>& other) const
{
    return std::sqrt(
        (x_ - other.x_) * (x_ - other.x_) +
        (y_ - other.y_) * (y_ - other.y_) +
        (z_ - other.z_) * (z_ - other.z_)
    );
}

template<typename T>
inline bool Vector<T>::isZero(Scalar value, Scalar epsilon)
{
    return std::fabs(value) < epsilon;
}

template<typename T>
inline bool Vector<T>::isZeroVector(Scalar epsilon) const
{
    return isZero(x_, epsilon) && isZero(y_, epsilon) && isZero(z_, epsilon);
}
#endif // VECTOR_H_
