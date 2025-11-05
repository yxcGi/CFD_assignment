#include <iostream>
#include "Vector.hpp"
#include "Mesh.h"

int main()
{
    // 读取网格
    Mesh mesh("complex_mesh.msh");
    mesh.printMeshInfo();

    


    return 0;
}