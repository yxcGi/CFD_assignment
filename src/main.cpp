#include <iostream>
#include "Math/Vector.hpp"
#include "Geometry/Mesh.h"
#include "Math/Tensor.hpp"
#include <Field.hpp>
#include "SparseMatrix.hpp"
#include <cmath>


using Scalar = double;
int main()
{

#if 1
    try
    {
        Mesh mesh("/Users/yxc/Desktop/code/c++/CFD_assignment/tempFile/OpenFOAM_tutorials/pitzDailySteady/constant/polyMesh");

        SparseMatrix<Scalar> A_b(&mesh);

        

    }
    catch (std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
#endif





#if 0
    try
    {
        using Scalar = double;
        std::vector<std::vector<Scalar>> A{
            { 1, 0, 0, 0, 0 },
            { 0, 0, 3, 2, 0 },
            { 0, 2, 0, 0, 0 },
            { 0, 0, 3, 0, 1 },
            { 0, 9, 2, 0, 0 }
        };

        SparseMatrix<Scalar> A_sparse(A);
        A_sparse.printMatrix();
    }
    catch (std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
#endif

#if 0
    try
    {
        using Scalar = double;
        // 读取网格
        Mesh mesh("/Users/yxc/Desktop/code/c++/CFD_assignment/tempFile/OpenFOAM_tutorials/pitzDailySteady/constant/polyMesh");

        // 创建标量场
        Field<Scalar> phi("T", &mesh);

        phi.setValue(
            [](Scalar x, Scalar y, Scalar z) {
                return std::sin(200 * x * x) * 200;
            }
        );

        phi.setBoundaryCondition("inlet", 1, 0, 300);
        phi.setBoundaryCondition("outlet", 1, 0, 300);
        phi.setBoundaryCondition("upperWall", 1, 0, 300);
        phi.setBoundaryCondition("lowerWall", 1, 0, 300);
        // phi.setBoundaryCondition("frontAndBack", 0, 1, 100);     // empty边界不需要设置边界条件
        phi.cellToFace();       // 若是第一步，只是将边界面的场根据边界条件进行更新

        // phi.cellToFace();
        phi.writeToFile("phi.dat");
    }
    catch (const std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
#endif 





    // cellField.cellToFace();
    // phi.setBoundaryCondition("const std::string &name", Scalar a, Scalar b, const double &c)


    // // 测试

    // using namespace std;

    // // 测试 Vector
    // Vector<double> v1(1.0, 2.0, 3.0);
    // Vector<double> v2(4.0, 5.0, 6.0);

    // cout << "v1 = " << v1 << endl;
    // cout << "v2 = " << v2 << endl;

    // cout << "v1 + v2 = " << v1 + v2 << endl;
    // cout << "v1 - v2 = " << v1 - v2 << endl;
    // cout << "v1 * 2 = " << v1 * 2.0 << endl;
    // cout << "v1 · v2 = " << (v1 & v2) << endl;
    // cout << "v1 × v2 = " << (v1 ^ v2) << endl;
    // cout << "|v1| = " << v1.magnitude() << endl;
    // cout << "unit(v1) = " << v1.unitVector() << endl;

    // // 测试 Tensor
    // Tensor<double> t1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    // Tensor<double> t2(9, 8, 7, 6, 5, 4, 3, 2, 1);

    // cout << "t1 = " << t1 << endl;
    // cout << "t2 = " << t2 << endl;
    // cout << "t1 + t2 = " << t1 + t2 << endl;
    // cout << "t1 * v1 = " << t1 * v1 << endl;

    // // 测试并矢
    // Tensor<Scalar> dyad = v1 * v2;
    // cout << "v1 ⊗ v2 = " << dyad << endl;

    // return 0;
}