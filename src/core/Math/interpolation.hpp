#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_





#include <iostream>


namespace interpolation
{
    enum class Scheme
    {
        // 线性差值
        LINEAR,
        // 算数平均
        ARITHMETIC_MEAN,
        // 调和平均
        HARMONIC_MEAN,
    };
}



/**
 * @brief 插值类函数对象，用于两点间的差值
 * @tparam T 数据类型
 */
template <typename Tp>
class Interpolation
{
    using Scalar = double;
public:
    // 传两个值就是算数平均
    // T operator()(const T& phi1, const T& phi2);
    Tp operator()(
        const Tp& phi1, const Tp& phi2,
        interpolation::Scheme scheme = interpolation::Scheme::LINEAR,
        Scalar alpha = 0.5) const;
};



template<typename Tp>
inline Tp Interpolation<Tp>::operator()(
    const Tp& phi1, const Tp& phi2,
    interpolation::Scheme scheme,
    Scalar alpha) const
{
    // phi1 * (1 - alpha) + phi2 * alpha
    if (scheme == interpolation::Scheme::LINEAR)
    {
        return phi1 + alpha * (phi2 - phi1);
    }
    else if (scheme == interpolation::Scheme::HARMONIC_MEAN)
    {
        return phi1 * phi2 / (phi2 * alpha + phi1 * (1 - alpha));
    }
    else if (scheme == interpolation::Scheme::ARITHMETIC_MEAN)
    {
        return (phi1 + phi2) * 0.5;
    }
    std::cerr << "interpolation::Scheme not supported!" << std::endl;
    throw std::runtime_error("interpolation::Scheme not supported!");
}



#endif // INTERPOLATION_H_
