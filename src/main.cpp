#include <iostream>
#include "Math/Vector.hpp"
#include "Geometry/Mesh.h"
#include "Math/Tensor.hpp"
#include "Fields/CellField.hpp"

int main()
{

    
    using Scalar = double;
    // 读取网格
    // Mesh mesh("/Users/yxc/Desktop/code/c++/CFD_assignment/tempFile/OpenFOAM_tutorials/pitzDailySteady/constant/polyMesh");

    // mesh.writeMeshToFile("/Users/yxc/Desktop/code/c++/CFD_assignment/tempFile/outputPolyMesh");


    // 测试

    using namespace std;

    // 测试 Vector
    Vector<double> v1(1.0, 2.0, 3.0);
    Vector<double> v2(4.0, 5.0, 6.0);

    cout << "v1 = " << v1 << endl;
    cout << "v2 = " << v2 << endl;

    cout << "v1 + v2 = " << v1 + v2 << endl;
    cout << "v1 - v2 = " << v1 - v2 << endl;
    cout << "v1 * 2 = " << v1 * 2.0 << endl;
    cout << "v1 · v2 = " << (v1 & v2) << endl;
    cout << "v1 × v2 = " << (v1 ^ v2) << endl;
    cout << "|v1| = " << v1.magnitude() << endl;
    cout << "unit(v1) = " << v1.unitVector() << endl;

    // 测试 Tensor
    Tensor<double> t1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Tensor<double> t2(9, 8, 7, 6, 5, 4, 3, 2, 1);

    cout << "t1 = " << t1 << endl;
    cout << "t2 = " << t2 << endl;
    cout << "t1 + t2 = " << t1 + t2 << endl;
    cout << "t1 * v1 = " << t1 * v1 << endl;

    // 测试并矢
    Tensor<Scalar> dyad = v1 * v2;
    cout << "v1 ⊗ v2 = " << dyad << endl;

    return 0;
}