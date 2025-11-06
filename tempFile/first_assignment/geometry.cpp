#include "geometry.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>

Geometry::Geometry() {}

Geometry::~Geometry() {}

void Geometry::loadMesh(const std::string& filename) {
    readMSHFile(filename);
    computeCellProperties();
    buildFaceList();
    computeFaceProperties();
    printMeshInfo();
}

void Geometry::readMSHFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open mesh file " << filename << std::endl;
        return;
    }
    
    std::string line;
    bool inNodes = false, inElements = false;
    int nodeCount = 0, elementCount = 0;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;
        
        if (token == "$Nodes") {
            inNodes = true;
            continue;
        } else if (token == "$EndNodes") {
            inNodes = false;
            continue;
        } else if (token == "$Elements") {
            inElements = true;
            continue;
        } else if (token == "$EndElements") {
            inElements = false;
            continue;
        }
        
        if (inNodes && !token.empty()) {
            nodeCount = std::stoi(token);
            nodes.resize(static_cast<size_t>(nodeCount));
            std::cout << "Reading " << nodeCount << " nodes..." << std::endl;
            // 循环读取所有节点数据
            for (int i = 0; i < nodeCount; i++) {
                std::getline(file, line);
                std::istringstream nodeIss(line);
                int nodeId;
                double x, y, z;
                nodeIss >> nodeId >> x >> y >> z;
                // GMSH节点ID从1开始，我们改为从0开始
                nodes[nodeId - 1] = Point(x, y);
            }
        }
        
        if (inElements && !token.empty()) {
            elementCount = std::stoi(token);
            cells.reserve(static_cast<size_t>(elementCount));
            std::cout << "Reading " << elementCount << " elements..." << std::endl;

            // 循环读取所有单元数据
            for (int i = 0; i < elementCount; i++) {
                std::getline(file, line);
                std::cout << "Reading element line " << i+1 << ": " << line << std::endl;
                std::istringstream elemIss(line);
                int elementId, elementType, numTags;
                elemIss >> elementId >> elementType >> numTags;
                std::cout << "Parsed: id=" << elementId << ", type=" << elementType << ", numTags=" << numTags << std::endl;
                
                // Skip tags
                for (int j = 0; j < numTags; j++) {
                    int tag;
                    elemIss >> tag;
                }
                
                // Process triangles and quadrilaterals
                // In GMSH format: type 1 = line, type 2 = triangle, type 3 = quadrilateral
                if (elementType == 2) {  // Triangle
                    Cell cell;
                    cell.id = static_cast<int>(cells.size());
                    cell.type = 3;  // Internal type: 3 = triangle
                    
                    for (int k = 0; k < 3; k++) {
                        int nodeId;
                        elemIss >> nodeId;
                        cell.nodes.push_back(nodeId - 1);  // GMSH starts from 1, we change to start from 0
                    }
                    
                    cells.push_back(cell);
                    std::cout << "Added triangle element " << elementId << " with nodes: " 
                              << cell.nodes[0] << " " << cell.nodes[1] << " " << cell.nodes[2] << std::endl;
                } else if (elementType == 3) {  // Quadrilateral
                    Cell cell;
                    cell.id = static_cast<int>(cells.size());
                    cell.type = 4;  // Internal type: 4 = quadrilateral
                    
                    for (int k = 0; k < 4; k++) {
                        int nodeId;
                        elemIss >> nodeId;
                        cell.nodes.push_back(nodeId - 1);  // GMSH starts from 1, we change to start from 0
                    }
                    
                    cells.push_back(cell);
                    std::cout << "Added quadrilateral element " << elementId << " with nodes: " 
                              << cell.nodes[0] << " " << cell.nodes[1] << " " << cell.nodes[2] << " " << cell.nodes[3] << std::endl;
                }
            }
        }
    }
    
    file.close();
    std::cout << "Successfully loaded mesh: " << nodes.size() << " nodes, " << cells.size() << " cells" << std::endl;
}

void Geometry::computeCellProperties() {
    for (auto& cell : cells) {
        // 计算单元中心
        double cx = 0.0, cy = 0.0;
        for (int nodeId : cell.nodes) {
            cx += nodes[nodeId].x;
            cy += nodes[nodeId].y;
        }
        cell.centroid.x = cx / cell.nodes.size();
        cell.centroid.y = cy / cell.nodes.size();
        
        // 计算单元面积
        if (cell.type == 3) {  // 三角形
            const Point& p1 = nodes[cell.nodes[0]];
            const Point& p2 = nodes[cell.nodes[1]];
            const Point& p3 = nodes[cell.nodes[2]];
            cell.area = 0.5 * std::abs((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y));
        } else if (cell.type == 4) {  // 四边形
            // 使用三角形分解计算面积
            const Point& p1 = nodes[cell.nodes[0]];
            const Point& p2 = nodes[cell.nodes[1]];
            const Point& p3 = nodes[cell.nodes[2]];
            const Point& p4 = nodes[cell.nodes[3]];
            
            double area1 = 0.5 * std::abs((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y));
            double area2 = 0.5 * std::abs((p3.x - p1.x) * (p4.y - p1.y) - (p4.x - p1.x) * (p3.y - p1.y));
            cell.area = area1 + area2;
        }
        
        // Initialize neighbor cell list
        cell.neighbors.clear();
    }
}

void Geometry::buildFaceList() {
    std::map<std::pair<int, int>, int> edgeToFace;
    
    for (const auto& cell : cells) {
        int numNodes = cell.nodes.size();
        for (int i = 0; i < numNodes; i++) {
            int node1 = cell.nodes[i];
            int node2 = cell.nodes[(i + 1) % numNodes];
            
            // Ensure consistent edge direction (smaller node ID first)
            if (node1 > node2) std::swap(node1, node2);
            
            std::pair<int, int> edge = {node1, node2};
            
            if (edgeToFace.find(edge) == edgeToFace.end()) {
                // Create new face
                Face face;
                face.id = static_cast<int>(faces.size());
                face.node1 = node1;
                face.node2 = node2;
                face.cell1 = cell.id;
                face.cell2 = -1;  // Default to boundary
                
                faces.push_back(face);
                edgeToFace[edge] = face.id;
            } else {
                // Face already exists, set adjacent cell
                int faceId = edgeToFace[edge];
                faces[faceId].cell2 = cell.id;
            }
        }
    }
    
    // Update cell neighbor information
    for (const auto& face : faces) {
        if (face.cell1 != -1) {
            cells[face.cell1].neighbors.push_back(face.cell2);
        }
        if (face.cell2 != -1) {
            cells[face.cell2].neighbors.push_back(face.cell1);
        }
    }
}

void Geometry::computeFaceProperties() {
    for (auto& face : faces) {
        // Calculate face center
        face.centroid.x = 0.5 * (nodes[face.node1].x + nodes[face.node2].x);
        face.centroid.y = 0.5 * (nodes[face.node1].y + nodes[face.node2].y);
        
        // Calculate face length
        double dx = nodes[face.node2].x - nodes[face.node1].x;
        double dy = nodes[face.node2].y - nodes[face.node1].y;
        face.length = std::sqrt(dx * dx + dy * dy);
        
        // Calculate face normal (pointing to cell1)
        face.normal.x = dy / face.length;
        face.normal.y = -dx / face.length;
        
        // Ensure normal points in correct direction
        if (face.cell1 != -1 && face.cell2 != -1) {
            Point cell1Center = cells[face.cell1].centroid;
            Point cell2Center = cells[face.cell2].centroid;
            
            Point faceToCell1 = {cell1Center.x - face.centroid.x, cell1Center.y - face.centroid.y};
            Point faceToCell2 = {cell2Center.x - face.centroid.x, cell2Center.y - face.centroid.y};
            
            double dot1 = face.normal.x * faceToCell1.x + face.normal.y * faceToCell1.y;
            double dot2 = face.normal.x * faceToCell2.x + face.normal.y * faceToCell2.y;
            
            if (dot1 < dot2) {
                face.normal.x = -face.normal.x;
                face.normal.y = -face.normal.y;
            }
        }
    }
}

void Geometry::printMeshInfo() const {
    std::cout << "\n=== Mesh Information ===" << std::endl;
    std::cout << "Number of nodes: " << nodes.size() << std::endl;
    std::cout << "Number of cells: " << cells.size() << std::endl;
    std::cout << "Number of faces: " << faces.size() << std::endl;
    
    int triangleCount = 0, quadCount = 0;
    for (const auto& cell : cells) {
        if (cell.type == 3) triangleCount++;
        else if (cell.type == 4) quadCount++;
    }
    
    std::cout << "Number of triangles: " << triangleCount << std::endl;
    std::cout << "Number of quadrilaterals: " << quadCount << std::endl;
    
    double totalArea = 0.0;
    for (const auto& cell : cells) {
        totalArea += cell.area;
    }
    std::cout << "Total domain area: " << totalArea << std::endl;
}

Point Geometry::getCellCentroid(int cellId) const {
    return cells[cellId].centroid;
}

double Geometry::getCellArea(int cellId) const {
    return cells[cellId].area;
}

Point Geometry::getFaceCentroid(int faceId) const {
    return faces[faceId].centroid;
}

Point Geometry::getFaceNormal(int faceId) const {
    return faces[faceId].normal;
}

double Geometry::getFaceLength(int faceId) const {
    return faces[faceId].length;
}

void Geometry::setBoundaryCondition(int boundaryId, const std::string& bcType) {
    boundary_conditions[boundaryId] = bcType;
}

std::string Geometry::getBoundaryCondition(int boundaryId) const {
    auto it = boundary_conditions.find(boundaryId);
    if (it != boundary_conditions.end()) {
        return it->second;
    }
    return "wall";  // 默认边界条件
}
