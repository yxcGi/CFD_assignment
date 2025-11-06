#include "Mesh.h"
#include <sstream>


using ULL = unsigned long long;

Mesh::Mesh() : isValid_(false)
{}

Mesh::Mesh(const std::string& path) : isValid_(false), meshPath_(path)
{
    readMesh(path);
}

void Mesh::readMesh(const std::string& path)
{
    std::string Mesh_points = path + "/points";
    std::string Mesh_faces = path + "/faces";
    std::string Mesh_boundary = path + "/boundary";
    std::string Mesh_neighbour = path + "/neighbour";;
    std::string Mesh_owner = path + "/owner";

    readPoints(Mesh_points);
    readFaces(Mesh_faces, Mesh_boundary, Mesh_neighbour);

}

void Mesh::readPoints(const std::string& pointsPath)
{
    std::ifstream ifsP(pointsPath);
    if (!ifsP.is_open())
    {
        std::cerr << "Error: Unable to open points file: " << pointsPath << std::endl;
        throw std::runtime_error("Points file open error");
    }

    // 处理points文件
    std::string line;
    std::cout << "Reading points file..." << std::endl;

    // 跳过文件头, 注释和FoamFile
    while (std::getline(ifsP, line) && line.find("FoamFile") != std::string::npos) {}
    while (std::getline(ifsP, line) && line != "}") {}

    ULL pointNum = 0;
    // 记录点数，跳过(
    while (std::getline(ifsP, line) && line != "(")
    {
        // 若能转为数字
        if (std::isdigit(line[0]))
        {
            pointNum = std::stoull(line);
        }
    }

    points_.reserve(pointNum);   // 提前分配空间

    // std::getline(ifsP, line);

    // 读取坐标
    char ignore;    // 接收括号和逗号
    int px, py, pz;
    std::istringstream iss;
    while (std::getline(ifsP, line) && line != ")")
    {
        iss.str(line);
        iss >> ignore >> px >> ignore >> py >> ignore >> pz;
        points_.emplace_back(px, py, pz);
        iss.clear();
    }


    /* ==========测试========= */

    // for (const auto& ele : points_)
    // {
    //     std::cout << ele << std::endl;
    // }
}

void Mesh::readFaces(const std::string& facesPath, const std::string& ownerPath, const std::string& neighbourPath)
{}
