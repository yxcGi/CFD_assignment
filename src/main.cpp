#include <iostream>
#include "Math/Vector.hpp"
#include "Geometry/Mesh.h"

int main()
{
    // 读取网格
    Mesh mesh("/Users/yxc/Desktop/code/c++/CFD_assignment/tempFile/OpenFOAM_tutorials/pitzDailySteady/constant/polyMesh");

    mesh.writeMeshToFile("/Users/yxc/Desktop/code/c++/CFD_assignment/tempFile/outputPolyMesh");


    


    return 0;
}