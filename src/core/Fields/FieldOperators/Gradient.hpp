#pragma once


// #include "Field.hpp"
#include "Vector.hpp"
#include "Tensor.hpp"
#include "CellField.hpp"
#include "FaceField.hpp"
// #include "Field.hpp"

template <typename Tp>
class Field;

// 定义计算梯度的方式，高斯格林公式、最小二乘
enum class GradientMethod
{
    GAUSS_GREEN,
    LEAST_SQUARES
};


// 计算场的梯度，默认为高斯格林公式
template <typename Tp>
inline auto grad(
    const Field<Tp>& field,
    GradientMethod method = GradientMethod::GAUSS_GREEN
) -> CellField<decltype(Tp() * Vector<Scalar>())>
{
    using ULL = unsigned long long;
    if (!field.getFaceField().isValid())        // 场必须有效
    {
        std::cerr << "Error: face field is not valid." << std::endl;
        throw std::runtime_error("face field is not valid.");
    }
    CellField<decltype(Tp()* Vector<Scalar>())> resultField(field.getName(), field.getMesh());
    // 取出当前面场，遍历每个cell, 计算梯度
    const FaceField<Tp>& currentFaceField = field.getFaceField();
    const std::vector<Cell>& cells = field.getMesh()->getCells();
    const std::vector<Face>& faces = field.getMesh()->getFaces();
    if (method == GradientMethod::GAUSS_GREEN)
    {

        for (int i = 0; i < cells.size(); ++i)  // 计算每个单元的梯度并赋值
        {
            const Cell& cell = cells[i];
            // 获取当前单元各个面的id
            const std::vector<ULL>& faceIds = cell.getFaceIndices();
            // 定义总的phi * S_f
            decltype(Tp() * Vector<Scalar>()) total_Phi_Sf{};
            for (ULL j : faceIds)       // 对每个面的Phi * S_f加和
            {
                const Face& face = faces[j];
                Vector<Scalar> Sf = face.getArea() * face.getNormal();
                Tp phi = currentFaceField[j];
                total_Phi_Sf += phi * Sf;
            }
            resultField[i] = total_Phi_Sf / cell.getVolume();
        }
        return resultField;
    }
    // else if (method == GradientMethod::LEAST_SQUARES)     // 挖坑
    // {
    // }
    return resultField;
}