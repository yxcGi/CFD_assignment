#include "Mesh.h"

Mesh::Mesh() : isValid_(false)
{}

Mesh::Mesh(const std::string path)
{
    loadMesh(path);
}

void Mesh::loadMesh(const std::string& path)
{
    std::string points = path + "/points";
    std::string faces = path + "/faces";
    std::string boundary = path + "/boundary";
    std::string neighbour = path + "/neighbour";;
    std::string owner = path + "/owner";

    std::ifstream ifsP(points);
    std::ifstream ifsF(faces);
    std::ifstream ifsB(boundary);
    std::ifstream ifsN(neighbour);
    std::ifstream ifsO(owner);


    if (!ifsP.is_open() ||
        !ifsF.is_open() ||
        !ifsB.is_open() ||
        !ifsN.is_open() ||
        !ifsO.is_open())
    {
        std::cerr << "Error: Unable to open mesh files in path: " << path << std::endl;
        throw std::runtime_error("Mesh file open error");
    }

 

}
