#ifndef LAPLICIAN_H_
#define LAPLICIAN_H_

#include "SparseMatrix.hpp"
#include "Field.hpp"

namespace fvm
{
    using Scalar = double;

    template<typename Tp>
    void Laplician(
        SparseMatrix<Tp>& matrix,   // 稀疏矩阵
        FaceField<Tp>& gamma,        // 扩散系数
        Field<Tp>& phi,              // 场值
        CellField<decltype(Tp()* Vector<Scalar>())>& gradientField       // 单元梯度场
    );

    template<typename Tp>
    void Laplician(
        SparseMatrix<Tp>& matrix,
        FaceField<Tp>& gamma,
        Field<Tp>& phi,
        CellField<decltype(Tp()* Vector<Scalar>())>& gradientField
    )
    {
        
    }

}




#endif // LAPLICIAN_H_