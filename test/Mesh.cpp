#include "Mesh.h"
#include <sstream>
#include "Face.h"
#include <sstream>
#include <map>
#include <cassert>
#include <thread>


using ULL = unsigned long long;
using LL = long long;

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
    readBoundaryPatch(Mesh_boundary);
    readFaces(Mesh_faces, Mesh_owner, Mesh_neighbour);
    buildCellsFromFaces();     // 建立网格拓扑关系
    calculateMeshInfo();       // 计算网格信息


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
        if (!line.empty() && std::isdigit(line[0]))
        {
            pointNum = std::stoull(line);
        }
    }

    points_.reserve(pointNum);   // 提前分配空间

    // std::getline(ifsP, line);

    // 读取坐标
    char ignore;    // 接收括号和逗号
    Scalar px, py, pz;
    std::istringstream iss;
    while (std::getline(ifsP, line) && line != ")")
    {
        iss.str(line);
        iss >> ignore >> px >> py >> pz;
        // std::cout << "Read point: (" << px << ", " << py << ", " << pz << ")" << std::endl;
        // std::this_thread::sleep_for(std::chrono::milliseconds(10));
        points_.emplace_back(px, py, pz);
        iss.clear();
    }


    /* ==========测试========= */

    // for (const auto& ele : points_)
    // {
    //     std::cout << ele << std::endl;
    // }
}

void Mesh::readBoundaryPatch(const std::string& boundaryPath)
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
    while (std::getline(ifsB, line) && line != "}") {}

    // 找 ( 括号，并记录边界类型数量
    int boundaryTypeNum = 0;
    std::string boundaryName;                   // 边界名称
    BoundaryPatch::BoundaryType type;           // 边界类型
    std::string typeStr;
    std::string ignoreStr;                     // 边界类型字符串(临时存储)
    ULL nFaces = 0;                             // 边界面数量
    ULL startFace = 0;                          // 起始面索引
    std::istringstream iss;
    while (std::getline(ifsB, line) && line != "(")
    {
        if (!line.empty() && std::isdigit(line[0]))
        {
            boundaryTypeNum = std::stoi(line);
        }
    }
    // boundaryPatch_.reserve(boundaryTypeNum);
    // 读取边界名称：inlet、outlet、upperWall

    for (int i = 0; i < boundaryTypeNum; ++i)
    {
        // 找{括号，读取上一行的名称name
        while (std::getline(ifsB, line) && line.find("{") == std::string::npos)
        {
            iss.str(line);
            iss >> boundaryName;
        }

        // 寻找type关键字
        while (std::getline(ifsB, line) && line.find("type") == std::string::npos) {}
        line.pop_back();
        iss.clear();
        iss.str(line);
        iss >> ignoreStr >> typeStr;        // 第二个line得到类型的字符串加分号，第一个得到type字符串
        // system("pause");
        type = stringToType(typeStr);

        // 寻找nFaces关键字
        while (std::getline(ifsB, line) && line.find("nFaces") == std::string::npos) {}
        line.pop_back();        // 去掉分号
        iss.clear();
        iss.str(line);
        iss >> ignoreStr >> nFaces;

        // 寻找startFace关键字
        while (std::getline(ifsB, line) && line.find("startFace") == std::string::npos) {}
        line.pop_back();
        iss.clear();
        iss.str(line);
        iss >> ignoreStr >> startFace;

        // 添加至boundaryPatch_ map中
        BoundaryPatch patch(boundaryName, nFaces, startFace, type);
        boundaryPatch_.emplace(boundaryName, std::move(patch));

        /* 测试用 */
        // std::cout << "Boundary Patch " << i + 1 << ": " << std::endl;
        // std::cout << "  Name: " << boundaryName << std::endl;
        // std::cout << "  Type: " << static_cast<int>(type) << std::endl;
        // std::cout << "  nFaces: " << nFaces << std::endl;
        // std::cout << "  startFace: " << startFace << std::endl;
    }
}

void Mesh::readFaces(const std::string& facesPath, const std::string& ownerPath, const std::string& neighbourPath)
{
    std::ifstream ifsF(facesPath, std::ios::in);
    std::ifstream ifsO(ownerPath, std::ios::in);
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
        if (!line.empty() && std::isdigit(line[0]))
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



    // 开始处理neighbour文件
    std::map<ULL, ULL> boundaryIndxRange;   // 边界面的起始和终止索引
    // 获取边界面的区间范围
    for (const auto& ele : boundaryPatch_)
    {
        boundaryIndxRange[ele.second.getStartFace()] = ele.second.getNFace() + ele.second.getStartFace();
    }
    std::vector<ULL> internalFaceIndices;    // 内部面索引列表
    internalFaceIndices.reserve(faces_.size()); // 分配大小

    // 首个内部区间
    int cycleIndex = 1;     // 记录循环次数
    ULL lastEnd = 0;        // 记录上一个边界面区间的结束索引

    for (const auto& [begin, end] : boundaryIndxRange)
    {
        if (cycleIndex == 1)    // 首个内部区间为0到第一个外部区间的左值
        {
            for (ULL i = 0; i < begin; ++i)
            {
                internalFaceIndices.push_back(i);
            }
            lastEnd = end;
        }
        else    // 中间的区间
        {
            for (ULL i = lastEnd; i < begin; ++i)
            {
                internalFaceIndices.push_back(i);
            }
            lastEnd = end;
        }
        cycleIndex++;
    }
    for (ULL i = lastEnd; i < faces_.size(); ++i)    // 最后一个内部面区间
    {
        internalFaceIndices.push_back(i);
    }

    // 将所有内部面的neighbor信息写入faces_
    readNeighbour(neighbourPath, internalFaceIndices);

    
    // 测试用
    // for (const auto& face : faces_)
    // {
    //     face.printFaceInfo();
    //     // 睡100毫秒
    //     std::this_thread::sleep_for(std::chrono::milliseconds(100));
    // }
}

void Mesh::readNeighbour(const std::string& neighbourPath, std::vector<ULL> internalFaceIndices)
{
    std::ifstream ifsN(neighbourPath);
    if (!ifsN.is_open())
    {
        std::cerr << "Error: Unable to open neighbour file: " << neighbourPath << std::endl;
        throw std::runtime_error("Neighbour file open error");
    }

    // 处理neighbour文件
    std::string line;
    std::cout << "Reading neighbour file..." << std::endl;
    // 跳过neighbour文件头, 注释和FoamFile
    while (std::getline(ifsN, line) && line.find("FoamFile") == std::string::npos) {}
    while (std::getline(ifsN, line) && line != "}") {}

    // 读取neighbour单元的个数
    ULL neighbourNum = 0;
    while (std::getline(ifsN, line) && line != "(")
    {
        if (!line.empty() && std::isdigit(line[0]))
        {
            neighbourNum = std::stoull(line);
        }
    }

    // assert(neighbourNum == faces_.size()- internalFaceIndices.size());     // 数量应相等
    int cycleIndex = 0;     // 循环次数
    while (std::getline(ifsN, line) && line != ")")
    {
        faces_[internalFaceIndices[cycleIndex]].setNeighbor(
            std::stoull(line));
        cycleIndex++;
    }

    // // 测试用
    // for (auto& ele : faces_)
    // {
    //     ele.printFaceInfo();
    // }
}

void Mesh::buildCellsFromFaces()
{
    // 确定单元总数（owner里的最大索引）
    ULL maxCellIndex = 0;
    for (const auto& face : faces_)
    {
        maxCellIndex = std::max(maxCellIndex, face.getOwnerIndex());

        auto faceNeighborIndex = face.getNeighborIndex();
        if (faceNeighborIndex != -1)
        {
            maxCellIndex = std::max(maxCellIndex, static_cast<ULL>(faceNeighborIndex));
        }
    }

    cells_.resize(maxCellIndex + 1);    // 索引从0开始，大小为最大索引+1

    // 遍历所有面
    for (ULL i = 0; i < faces_.size(); ++i)
    {
        ULL ownerIndex = faces_[i].getOwnerIndex();
        LL neighborIndex = faces_[i].getNeighborIndex();

        if (neighborIndex == -1)    // 说明是边界面
        {
            cells_[ownerIndex].addFaceIndex(i);
            boundaryCellIndexes_.emplace(ownerIndex);
        }
        else                        // 内部面
        {
            cells_[ownerIndex].addFaceIndex(i);
            cells_[neighborIndex].addFaceIndex(i);
            internalCellIndexes_.emplace(ownerIndex);
            internalCellIndexes_.emplace(neighborIndex);
        }
    }

    // 测试用
    // std::cout << "Internal Cell Indices: ";
    // for (const auto& ele : cells_)
    // {
    //     ele.printCellInfo();
    // }
}

/* =========测试========= */
// for (const auto& ele : faces_)
// {
//     std::cout << ele << std::endl;
// }


BoundaryPatch::BoundaryType Mesh::stringToType(const std::string& name) const
{
    if (name == "patch")
    {
        return BoundaryPatch::BoundaryType::PATCH;
    }
    else if (name == "wall")
    {
        return BoundaryPatch::BoundaryType::WALL;
    }
    else if (name == "symmetry")
    {
        return BoundaryPatch::BoundaryType::SYMMETRY;
    }
    else if (name == "cyclic")
    {
        return BoundaryPatch::BoundaryType::CYCLIC;
    }
    else if (name == "wedge")
    {
        return BoundaryPatch::BoundaryType::WEDGE;
    }
    else if (name == "empty")
    {
        return BoundaryPatch::BoundaryType::EMPTY;
    }
    else if (name == "processor")
    {
        return BoundaryPatch::BoundaryType::PROCESSOR;
    }
    throw std::invalid_argument("Unknown boundary type: " + name);
}

void Mesh::calculateMeshInfo()
{
    // 计算面信息
    for (auto& face : faces_)
    {
        face.calculateFaceInfo(points_);
    }
}
