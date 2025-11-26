#include <iostream>
#include "Laplician.hpp"
#include "Vector.hpp"
#include "Geometry/Mesh.h"
#include "Tensor.hpp"
#include <Field.hpp>
#include "SparseMatrix.hpp"
#include <cmath>


using Scalar = double;
int main()
{


#if 0
    try
    {
        Mesh mesh("/Users/yxc/Desktop/code/c++/CFD_assignment/tempFile/OpenFOAM_tutorials/pitzDailySteady/constant/polyMesh");

        SparseMatrix<Scalar> A_b(&mesh);
        A_b.setValue(0, 0, 99);
        A_b.setValue(0, 1, 99);
        // A_b.setValue(0, 2, 99);
        // A_b.setValue(0, 3, 99);
        // A_b.setValue(0, 4, 99);

        // 打印第一行前5元素
        for (int i = 0; i < 5; i++)
        {
            std::cout << A_b.at(0, i) << " ";
        }
        std::cout << std::endl;

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
            { 0, 0, 0, 0, 0 },
            { 0, 1, 0, 0, 0 },
            { 0, 0, 0, 0, 0 },
            { 2, 0, 0, 0, 0 },
            { 0, 0, 0, 4, 0 }
        };

        std::vector<std::vector<Scalar>> B{
            { 0, 0 },
            { 1, 0 }
        };

        SparseMatrix<Scalar> B_sparse(B);
        B_sparse.printMatrix();
        B_sparse.compress();

        // B_sparse.setValue(0, 0, 99);
        // B_sparse.setValue(0, 1, 99);
        B_sparse.setValue(1, 0, 99);
        // B_sparse.setValue(1, 1, 99);
        
        B_sparse.printMatrix();

        
    }
    catch (std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
#endif

#if 1
    try
    {
        using Scalar = double;
        // 读取网格
        Mesh mesh("/Users/yxc/Desktop/code/c++/CFD_assignment/tempFile/OpenFOAM_tutorials/cavity/constant/polyMesh");

        // 创建标量场
        Field<Scalar> phi("T", &mesh);

        phi.setValue(
            [](Scalar x, Scalar y, Scalar z) {
                return 200 * (x * x + y * y);
            }
        );

        phi.setBoundaryCondition("fixedWalls", 0, 1, 0);
        phi.setBoundaryCondition("movingWall", 1, 0, 200);
        // phi.setBoundaryCondition("upperWall", 0, 1, 0);
        // phi.setBoundaryCondition("lowerWall", 0, 1, 0);
        // phi.setBoundaryCondition("frontAndBack", 0, 1, 100);     // empty边界不需要设置边界条件
        phi.cellToFace();       // 若是第一步，只是将边界面的场根据边界条件进行更新


        phi.writeToFile("phi.dat");
        
        Field<Vector<Scalar>> gradPhi(grad(phi));
        // gradPhi.writeToFile("phi.dat");

        SparseMatrix<Scalar> A_b(&mesh);
        FaceField<Scalar> gamma("gamma", &mesh);
        gamma.setValue(10);

        fvm::Laplician(A_b, gamma, phi);

        // getchar();
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