#pragma once

#include <iostream>
#include "Vector.hpp"


template <typename T>
class Tensor
{
public:
    Tensor();
    Tensor(     // 3D
        T xx, T xy, T xz,
        T yx, T yy, T yz,
        T zx, T zy, T zz
    );
    Tensor(     // 2D
        T xx, T xy,
        T yx, T yy
    );
    Tensor(const Tensor<T>& src);
    Tensor(Tensor<T>&&) noexcept = default;
    Tensor<T>& operator=(const Tensor<T>& src);
    Tensor<T>& operator=(Tensor<T>&&) noexcept = default;

    ~Tensor() = default;
public:


private:
    T xx_, xy_, xz_;
    T yx_, yy_, yz_;
    T zx_, zy_, zz_;
};