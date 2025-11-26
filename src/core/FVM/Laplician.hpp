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
        using ULL = unsigned long long;
        const Mesh* mesh = matrix.getMesh();


        // 判断输入参数是否有效
        if (!matrix.isValid())
        {
            std::cerr << "Laplician() Error: matrix is not valid." << std::endl;
            throw std::invalid_argument("matrix is not valid.");
        }
        if (!gamma.isValid())
        {
            std::cerr << "Laplician() Error: gamma is not valid." << std::endl;
            throw std::invalid_argument("gamma is not valid.");
        }
        if (!phi.isValid())
        {
            std::cerr << "Laplician() Error: phi is not valid." << std::endl;
            throw std::invalid_argument("phi is not valid.");
        }
        if (!gradientField.isValid())
        {
            std::cerr << "Laplician() Error: gradientField is not valid." << std::endl;
            throw std::invalid_argument("gradientField is not valid.");
        }
        // 判断是否矩阵是否为空指针网格（稀疏矩阵只能由mesh构造）
        if (mesh == nullptr)
        {
            std::cerr << "Laplician() Error: mesh is nullptr." << std::endl;
            throw std::invalid_argument("mesh is nullptr.");
        }
        // 判断输出参数是否为同一网格
        if (mesh != gamma.getMesh() ||
            mesh != phi.getMesh() ||
            mesh != gradientField.getMesh())
        {
            std::cerr << "Laplician() Error: The mesh of the input fields parameters is not the same as the mesh of the matrix." << std::endl;
            throw std::invalid_argument("The mesh of the input fields parameters is not the same as the mesh of the matrix.");
        }


        // 开始对面进行遍历，内部面与边界面分开
        
    }
}




#endif // LAPLICIAN_H_