#include "Mesh.h"
#include <sstream>
#include "Face.h"


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
    readBoundaryConditions(Mesh_boundary);
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
    while (std::getline(ifsP, line) && line.find("FoamFile") == std::string::npos) {}
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

void Mesh::readBoundaryConditions(const std::string& boundaryPath)
{
    std::ifstream ifsB(boundaryPath);
    if (!ifsB.is_open())
    {
        std::cerr << "Error: Unable to open boundary file: " << boundaryPath << std::endl;
        throw std::runtime_error("Boundary file open error");
    }

    // 处理boundary文件
    std::string line;
    std::cout << "Reading boundary file..." << std::endl;

    // 跳过文件头, 注释和FoamFile
    while (std::getline(ifsB, line) && line.find("FoamFile") == std::string::npos) {}

    

}

void Mesh::readFaces(const std::string& facesPath, const std::string& ownerPath, const std::string& neighbourPath)
{
    std::ifstream ifsF(facesPath);
    std::ifstream ifsO(ownerPath);
    if (!ifsF.is_open())
    {
        std::cerr << "Error: Unable to open faces file: " << facesPath << std::endl;
        throw std::runtime_error("Faces file open error");
    }
    if (!ifsO.is_open())
    {
        std::cerr << "Error: Unable to open owner file: " << ownerPath << std::endl;
        throw std::runtime_error("Owner file open error");
    }

    // 处理faces文件
    std::string line;
    std::cout << "Reading faces file..." << std::endl;
    // 跳过faces文件头, 注释和FoamFile
    while (std::getline(ifsF, line) && line.find("FoamFile") == std::string::npos) {}
    while (std::getline(ifsF, line) && line != "}") {}

    // 处理owner文件
    // 跳过owner文件头, 注释和FoamFile
    while (std::getline(ifsO, line) && line.find("FoamFile") == std::string::npos) {}
    while (std::getline(ifsO, line) && line != "}") {}


    ULL faceNum = 0;
    while (std::getline(ifsF, line) && line != "(")
    {
        if (std::isdigit(line[0]))
        {
            faceNum = std::stoull(line);
        }
    }
    while (std::getline(ifsO, line) && line != "(") {}
    faces_.reserve(faceNum);

    // 开始读取面信息
    /*
    2(2 6)
    2(3 7)
    2(6 10)
    2(7 11)
    2(10 14)
    2(11 15)
    2(5 6)
    2(6 7)
     */
    char ignore;    // 接收括号
    int pointNum;
    ULL pIndex = 0;
    std::istringstream iss;
    while (std::getline(ifsF, line) && line != ")")
    {
        // 先处理faces文件
        std::vector<ULL> pointIndexs;   // 接收读取到的点索引
        iss.str(line);
        iss >> pointNum >> ignore;
        // std::cout << pointNum << ignore;
        // system("pause");
        // 创建Face对象
        pointIndexs.reserve(pointNum);
        for (int i = 0; i < pointNum; ++i)      // 读取点索引
        {
            iss >> pIndex;
            pointIndexs.push_back(pIndex);
        }

        // owner文件同步读取
        std::getline(ifsO, line);

        
        Face face(pointIndexs, std::stoull(line));
        faces_.emplace_back(std::move(face));
        

    }

    /* =========测试========= */
    // for (const auto& ele : faces_)
    // {
    //     std::cout << ele << std::endl;
    // }

}
