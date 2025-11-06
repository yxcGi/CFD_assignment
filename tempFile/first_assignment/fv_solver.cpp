#include "fv_solver.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

FVSolver::FVSolver(Geometry* geom) : geometry(geom) {
    // 设置默认参数
    params.convectionVelocityX = 1.0;
    params.convectionVelocityY = 0.0;
    params.diffusionCoefficient = 0.1;
    params.sourceTerm = 0.0;
    params.timeStep = 0.001;
    params.finalTime = 1.0;
    params.maxIterations = 1000;
    params.tolerance = 1e-6;
    params.reconstructionMethod = "gauss";
    params.fluxScheme = "upwind";
}

FVSolver::~FVSolver() {}

void FVSolver::loadParameters(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open parameter file " << filename << ", using default parameters" << std::endl;
        return;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string key, value;
        iss >> key >> value;
        
        if (key == "convection_velocity_x") {
            params.convectionVelocityX = std::stod(value);
        } else if (key == "convection_velocity_y") {
            params.convectionVelocityY = std::stod(value);
        } else if (key == "diffusion_coefficient") {
            params.diffusionCoefficient = std::stod(value);
        } else if (key == "source_term") {
            params.sourceTerm = std::stod(value);
        } else if (key == "time_step") {
            params.timeStep = std::stod(value);
        } else if (key == "final_time") {
            params.finalTime = std::stod(value);
        } else if (key == "max_iterations") {
            params.maxIterations = std::stoi(value);
        } else if (key == "tolerance") {
            params.tolerance = std::stod(value);
        } else if (key == "reconstruction_method") {
            params.reconstructionMethod = value;
        } else if (key == "flux_scheme") {
            params.fluxScheme = value;
        }
    }
    
    file.close();
    std::cout << "Successfully loaded solver parameters" << std::endl;
}

void FVSolver::setBoundaryCondition(int boundaryId, const std::string& type, double value) {
    BoundaryCondition bc;
    bc.type = type;
    bc.value = value;
    boundaryConditions[boundaryId] = bc;
}

void FVSolver::solve() {
    std::cout << "\n=== Starting Solution ===" << std::endl;
    std::cout << "Time step: " << params.timeStep << std::endl;
    std::cout << "Final time: " << params.finalTime << std::endl;
    std::cout << "Reconstruction method: " << params.reconstructionMethod << std::endl;
    std::cout << "Flux scheme: " << params.fluxScheme << std::endl;
    
    initializeSolution();
    
    double currentTime = 0.0;
    int timeStep = 0;
    
    while (currentTime < params.finalTime) {
        // 保存上一时间步的解
        solutionOld = solution;
        
        // 计算梯度
        computeGradients();
        
        // 时间推进
        updateSolution();
        
        currentTime += params.timeStep;
        timeStep++;
        
        // Output progress
        if (timeStep % 100 == 0) {
            std::cout << "Time step: " << timeStep << ", Time: " << currentTime 
                      << ", Max solution: " << *std::max_element(solution.begin(), solution.end()) << std::endl;
        }
        
        // Output Tecplot files
        if (timeStep % 500 == 0) {
            std::string filename = "solution_" + std::to_string(timeStep) + ".dat";
            writeTecplotOutput(filename, timeStep);
        }
    }
    
    // Output final result
    writeTecplotOutput("final_solution.dat", timeStep);
    std::cout << "\nSolution completed!" << std::endl;
}

void FVSolver::initializeSolution() {
    int numCells = geometry->getNumCells();
    solution.resize(numCells, 0.0);
    solutionOld.resize(numCells, 0.0);
    gradients.resize(numCells, std::vector<double>(2, 0.0));
    
    // 设置复杂的初始条件
    for (int i = 0; i < numCells; i++) {
        Point centroid = geometry->getCellCentroid(i);
        
        // 创建多个污染源点
        double concentration = 0.0;
        
        // 源点1：中心区域
        double dist1 = std::sqrt((centroid.x - 1.0) * (centroid.x - 1.0) + (centroid.y - 1.0) * (centroid.y - 1.0));
        if (dist1 < 0.3) {
            concentration += 1.0 * std::exp(-dist1 * dist1 / (0.1 * 0.1));
        }
        
        // 源点2：左下角
        double dist2 = std::sqrt((centroid.x - 0.5) * (centroid.x - 0.5) + (centroid.y - 0.5) * (centroid.y - 0.5));
        if (dist2 < 0.2) {
            concentration += 0.8 * std::exp(-dist2 * dist2 / (0.05 * 0.05));
        }
        
        // 源点3：右上角
        double dist3 = std::sqrt((centroid.x - 1.5) * (centroid.x - 1.5) + (centroid.y - 1.5) * (centroid.y - 1.5));
        if (dist3 < 0.25) {
            concentration += 0.6 * std::exp(-dist3 * dist3 / (0.08 * 0.08));
        }
        
        // 添加一些随机扰动
        concentration += 0.1 * std::sin(centroid.x * 3.14159) * std::cos(centroid.y * 3.14159);
        
        solution[i] = std::max(0.0, concentration);
    }
    
    std::cout << "Initial conditions set" << std::endl;
}

void FVSolver::computeGradients() {
    if (params.reconstructionMethod == "gauss") {
        computeGaussGradients();
    } else if (params.reconstructionMethod == "least_squares") {
        computeLeastSquaresGradients();
    }
}

void FVSolver::computeGaussGradients() {
    int numCells = geometry->getNumCells();
    
    for (int i = 0; i < numCells; i++) {
        gradients[i][0] = 0.0;  // ∂φ/∂x
        gradients[i][1] = 0.0;  // ∂φ/∂y
        
        const Cell& cell = geometry->getCells()[i];
        double cellArea = geometry->getCellArea(i);
        
        // 对每个面应用Gauss定理
        for (int faceIdx : cell.neighbors) {
            if (faceIdx == -1) continue;  // 边界面
            
            // 找到共享面
            int faceId = -1;
            for (int j = 0; j < geometry->getFaces().size(); j++) {
                const Face& face = geometry->getFaces()[j];
                if ((face.cell1 == i && face.cell2 == faceIdx) || 
                    (face.cell2 == i && face.cell1 == faceIdx)) {
                    faceId = j;
                    break;
                }
            }
            
            if (faceId != -1) {
                const Face& face = geometry->getFaces()[faceId];
                Point normal = geometry->getFaceNormal(faceId);
                double faceLength = geometry->getFaceLength(faceId);
                
                // 计算面中心处的解值（简单平均）
                double phiFace = 0.5 * (solution[i] + solution[faceIdx]);
                
                gradients[i][0] += phiFace * normal.x * faceLength;
                gradients[i][1] += phiFace * normal.y * faceLength;
            }
        }
        
        // 除以单元面积
        gradients[i][0] /= cellArea;
        gradients[i][1] /= cellArea;
    }
}

void FVSolver::computeLeastSquaresGradients() {
    int numCells = geometry->getNumCells();
    
    for (int i = 0; i < numCells; i++) {
        gradients[i][0] = 0.0;
        gradients[i][1] = 0.0;
        
        const Cell& cell = geometry->getCells()[i];
        Point center = geometry->getCellCentroid(i);
        
        // 构建最小二乘系统
        double A[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        double b[2] = {0.0, 0.0};
        
        for (int neighborId : cell.neighbors) {
            if (neighborId == -1) continue;
            
            Point neighborCenter = geometry->getCellCentroid(neighborId);
            double dx = neighborCenter.x - center.x;
            double dy = neighborCenter.y - center.y;
            double dphi = solution[neighborId] - solution[i];
            
            // 构建最小二乘矩阵
            A[0][0] += dx * dx;
            A[0][1] += dx * dy;
            A[1][0] += dx * dy;
            A[1][1] += dy * dy;
            
            b[0] += dx * dphi;
            b[1] += dy * dphi;
        }
        
        // 求解2x2线性系统
        double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        if (std::abs(det) > 1e-12) {
            gradients[i][0] = (A[1][1] * b[0] - A[0][1] * b[1]) / det;
            gradients[i][1] = (A[0][0] * b[1] - A[1][0] * b[0]) / det;
        }
    }
}

void FVSolver::updateSolution() {
    int numCells = geometry->getNumCells();
    
    // 隐式方法：形成线性方程组 A * phi^{n+1} = b
    // 其中 A = I - dt * J, b = phi^n + dt * S
    
    // 计算梯度（用于二阶迎风格式）
    computeLeastSquaresGradients();
    
    // 创建系数矩阵A和右端向量b
    std::vector<std::vector<double>> A(numCells, std::vector<double>(numCells, 0.0));
    std::vector<double> b(numCells, 0.0);
    
    // 构建线性方程组
    buildImplicitSystem(A, b);
    
    // 求解线性方程组
    solveLinearSystem(A, b, solution);
}

void FVSolver::buildImplicitSystem(std::vector<std::vector<double>>& A, std::vector<double>& b) {
    int numCells = geometry->getNumCells();
    int numFaces = geometry->getNumFaces();
    
    // 初始化矩阵A为单位矩阵，向量b为上一时间步的解
    for (int i = 0; i < numCells; i++) {
        A[i][i] = 1.0;  // 对角线元素
        b[i] = solutionOld[i];  // 右端向量
    }
    
    // 计算雅可比矩阵的贡献（对流和扩散项）
    for (int faceId = 0; faceId < numFaces; faceId++) {
        const Face& face = geometry->getFaces()[faceId];
        
        if (face.cell2 == -1) {
            // 边界面的处理
            buildBoundaryFaceJacobian(faceId, A, b);
        } else {
            // 内部面的处理
            buildInternalFaceJacobian(faceId, A, b);
        }
    }
    
    // 添加源项贡献
    for (int i = 0; i < numCells; i++) {
        double cellArea = geometry->getCellArea(i);
        b[i] += params.timeStep * params.sourceTerm * cellArea / cellArea;
    }
}

void FVSolver::buildInternalFaceJacobian(int faceId, std::vector<std::vector<double>>& A, std::vector<double>& b) {
    const Face& face = geometry->getFaces()[faceId];
    Point normal = geometry->getFaceNormal(faceId);
    double faceLength = geometry->getFaceLength(faceId);
    
    int leftCellId = face.cell1;
    int rightCellId = face.cell2;
    double leftArea = geometry->getCellArea(leftCellId);
    double rightArea = geometry->getCellArea(rightCellId);
    
    // 对流项的雅可比贡献
    double uDotN = params.convectionVelocityX * normal.x + params.convectionVelocityY * normal.y;
    
    if (params.fluxScheme == "upwind") {
        if (uDotN > 0) {
            // 通量从左单元流向右单元
            double coeff = params.timeStep * uDotN * faceLength / leftArea;
            A[leftCellId][leftCellId] -= coeff;
            A[rightCellId][leftCellId] += coeff;
        } else {
            // 通量从右单元流向左单元
            double coeff = params.timeStep * (-uDotN) * faceLength / rightArea;
            A[rightCellId][rightCellId] -= coeff;
            A[leftCellId][rightCellId] += coeff;
        }
    } else if (params.fluxScheme == "central") {
        // 对流项的中心差分
        double leftCoeff = params.timeStep * uDotN * faceLength * 0.5 / leftArea;
        double rightCoeff = params.timeStep * uDotN * faceLength * 0.5 / rightArea;
        
        A[leftCellId][leftCellId] -= leftCoeff;
        A[leftCellId][rightCellId] -= leftCoeff;
        A[rightCellId][leftCellId] += rightCoeff;
        A[rightCellId][rightCellId] += rightCoeff;
    } else if (params.fluxScheme == "second_order_upwind") {
        // 延迟修正方法：隐式部分采用一阶格式，显式部分计算二阶修正
        
        // 隐式部分：使用一阶迎风格式的雅可比矩阵
        if (uDotN > 0) {
            // 通量从左单元流向右单元
            double coeff = params.timeStep * uDotN * faceLength / leftArea;
            A[leftCellId][leftCellId] -= coeff;
            A[rightCellId][leftCellId] += coeff;
        } else {
            // 通量从右单元流向左单元
            double coeff = params.timeStep * (-uDotN) * faceLength / rightArea;
            A[rightCellId][rightCellId] -= coeff;
            A[leftCellId][rightCellId] += coeff;
        }
        
        // 显式部分：计算二阶修正项并添加到右端向量b
        // 二阶通量 - 一阶通量
        double firstOrderFlux, secondOrderFlux;
        
        // 计算一阶迎风通量
        if (uDotN > 0) {
            firstOrderFlux = uDotN * solutionOld[leftCellId];
        } else {
            firstOrderFlux = uDotN * solutionOld[rightCellId];
        }
        
        // 计算二阶迎风通量
        if (uDotN > 0) {
            // 从左单元重构
            Point leftCenter = geometry->getCellCentroid(leftCellId);
            Point faceCenter = geometry->getFaceCentroid(faceId);
            
            double dx = faceCenter.x - leftCenter.x;
            double dy = faceCenter.y - leftCenter.y;
            
            // 使用梯度进行线性重构
            double phiFace = solutionOld[leftCellId] + gradients[leftCellId][0] * dx + gradients[leftCellId][1] * dy;
            secondOrderFlux = uDotN * phiFace;
        } else {
            // 从右单元重构
            Point rightCenter = geometry->getCellCentroid(rightCellId);
            Point faceCenter = geometry->getFaceCentroid(faceId);
            
            double dx = faceCenter.x - rightCenter.x;
            double dy = faceCenter.y - rightCenter.y;
            
            // 使用梯度进行线性重构
            double phiFace = solutionOld[rightCellId] + gradients[rightCellId][0] * dx + gradients[rightCellId][1] * dy;
            secondOrderFlux = uDotN * phiFace;
        }
        
        // 延迟修正：将修正项添加到右端向量
        double correction = secondOrderFlux - firstOrderFlux;
        double correctionCoeff = params.timeStep * faceLength;
        
        // 根据流动方向分配修正项
        if (uDotN > 0) {
            // 从左单元流向右单元，修正项影响两个单元
            b[leftCellId] -= correctionCoeff * correction / leftArea;
            b[rightCellId] += correctionCoeff * correction / rightArea;
        } else {
            // 从右单元流向左单元，修正项影响两个单元
            b[rightCellId] -= correctionCoeff * correction / rightArea;
            b[leftCellId] += correctionCoeff * correction / leftArea;
        }
    } else {
        // 默认使用一阶迎风
        if (uDotN > 0) {
            double coeff = params.timeStep * uDotN * faceLength / leftArea;
            A[leftCellId][leftCellId] -= coeff;
            A[rightCellId][leftCellId] += coeff;
        } else {
            double coeff = params.timeStep * (-uDotN) * faceLength / rightArea;
            A[rightCellId][rightCellId] -= coeff;
            A[leftCellId][rightCellId] += coeff;
        }
    }
    
    // 扩散项的雅可比贡献（使用中心差分）
    double diffCoeff = params.timeStep * params.diffusionCoefficient * faceLength;
    
    // 计算面中心到单元中心的距离
    Point leftCenter = geometry->getCellCentroid(leftCellId);
    Point rightCenter = geometry->getCellCentroid(rightCellId);
    Point faceCenter = geometry->getFaceCentroid(faceId);
    
    double dxLeft = faceCenter.x - leftCenter.x;
    double dyLeft = faceCenter.y - leftCenter.y;
    double dxRight = faceCenter.x - rightCenter.x;
    double dyRight = faceCenter.y - rightCenter.y;
    
    // 扩散通量的雅可比矩阵
    double leftDiffCoeff = diffCoeff / (leftArea * (dxLeft * normal.x + dyLeft * normal.y));
    double rightDiffCoeff = diffCoeff / (rightArea * (dxRight * normal.x + dyRight * normal.y));
    
    A[leftCellId][leftCellId] += leftDiffCoeff;
    A[leftCellId][rightCellId] -= leftDiffCoeff;
    A[rightCellId][leftCellId] -= rightDiffCoeff;
    A[rightCellId][rightCellId] += rightDiffCoeff;
}

void FVSolver::buildBoundaryFaceJacobian(int faceId, std::vector<std::vector<double>>& A, std::vector<double>& b) {
    const Face& face = geometry->getFaces()[faceId];
    Point normal = geometry->getFaceNormal(faceId);
    double faceLength = geometry->getFaceLength(faceId);
    
    int cellId = face.cell1;
    double cellArea = geometry->getCellArea(cellId);
    
    // 获取边界条件类型
    std::string bcType = "dirichlet";  // 默认Dirichlet边界条件
    
    // 对流项的边界处理
    double uDotN = params.convectionVelocityX * normal.x + params.convectionVelocityY * normal.y;
    
    if (bcType == "dirichlet") {
        // Dirichlet边界条件
        if (uDotN > 0) {
            // 流入边界：使用边界值
            double boundaryValue = 0.0;  // 这里应该使用实际的边界值
            b[cellId] += params.timeStep * uDotN * boundaryValue * faceLength / cellArea;
        } else {
            // 流出边界：使用单元值
            // 对流项不贡献到雅可比矩阵
        }
    } else if (bcType == "neumann") {
        // Neumann边界条件：给定通量
        double flux = 0.0;  // 这里应该使用实际的Neumann通量值
        b[cellId] += params.timeStep * flux * faceLength / cellArea;
    } else if (bcType == "wall") {
        // 壁面边界条件：零通量
        // 不贡献到雅可比矩阵
    }
    
    // 扩散项的边界处理
    if (bcType == "dirichlet") {
        // Dirichlet边界条件
        double boundaryValue = 0.0;  // 这里应该使用实际的边界值
        double diffCoeff = params.timeStep * params.diffusionCoefficient * faceLength / cellArea;
        
        // 计算距离
        Point cellCenter = geometry->getCellCentroid(cellId);
        Point faceCenter = geometry->getFaceCentroid(faceId);
        double distance = std::sqrt((faceCenter.x - cellCenter.x) * (faceCenter.x - cellCenter.x) + 
                                   (faceCenter.y - cellCenter.y) * (faceCenter.y - cellCenter.y));
        
        if (distance > 1e-12) {
            double coeff = diffCoeff / distance;
            A[cellId][cellId] += coeff;
            b[cellId] += coeff * boundaryValue;
        }
    } else if (bcType == "neumann") {
        // Neumann边界条件：给定梯度
        // 已经在对流项中处理了通量贡献
    } else if (bcType == "wall") {
        // 壁面边界条件：零扩散通量
        // 不贡献到雅可比矩阵
    }
}

void FVSolver::solveLinearSystem(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x) {
    int n = A.size();
    
    // 使用简单的Gauss-Seidel迭代法求解线性方程组
    int maxIter = 1000;
    double tolerance = 1e-8;
    
    // 初始化解为上一时间步的值
    x = solutionOld;
    
    for (int iter = 0; iter < maxIter; iter++) {
        double maxError = 0.0;
        std::vector<double> xNew = x;
        
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += A[i][j] * xNew[j];
                }
            }
            
            if (std::abs(A[i][i]) > 1e-12) {
                xNew[i] = (b[i] - sum) / A[i][i];
                maxError = std::max(maxError, std::abs(xNew[i] - x[i]));
            }
        }
        
        x = xNew;
        
        if (maxError < tolerance) {
            //std::cout << "Linear system converged in " << iter + 1 << " iterations, error: " << maxError << std::endl;
            break;
        }
        
        if (iter == maxIter - 1) {
            std::cout << "Warning: Linear system did not converge, max error: " << maxError << std::endl;
        }
    }
}

double FVSolver::computeConvectiveFlux(int faceId, double phiLeft, double phiRight) {
    const Face& face = geometry->getFaces()[faceId];
    Point normal = geometry->getFaceNormal(faceId);
    
    double uDotN = params.convectionVelocityX * normal.x + params.convectionVelocityY * normal.y;
    
    if (params.fluxScheme == "upwind") {
        if (uDotN > 0) {
            return uDotN * phiLeft;
        } else {
            return uDotN * phiRight;
        }
    } else if (params.fluxScheme == "central") {
        return uDotN * 0.5 * (phiLeft + phiRight);
    } else if (params.fluxScheme == "second_order_upwind") {
        // 二阶迎风格式：使用梯度信息重构面中心值
        if (uDotN > 0) {
            // 从左单元重构
            Point leftCenter = geometry->getCellCentroid(face.cell1);
            Point faceCenter = geometry->getFaceCentroid(faceId);
            
            // 计算从单元中心到面中心的向量
            double dx = faceCenter.x - leftCenter.x;
            double dy = faceCenter.y - leftCenter.y;
            
            // 使用预先计算好的梯度进行线性重构
            double phiFace = phiLeft + gradients[face.cell1][0] * dx + gradients[face.cell1][1] * dy;
            
            return uDotN * phiFace;
        } else {
            // 从右单元重构
            Point rightCenter = geometry->getCellCentroid(face.cell2);
            Point faceCenter = geometry->getFaceCentroid(faceId);
            
            // 计算从单元中心到面中心的向量
            double dx = faceCenter.x - rightCenter.x;
            double dy = faceCenter.y - rightCenter.y;
            
            // 使用预先计算好的梯度进行线性重构
            double phiFace = phiRight + gradients[face.cell2][0] * dx + gradients[face.cell2][1] * dy;
            
            return uDotN * phiFace;
        }
    } else {
        // 默认使用一阶迎风
        if (uDotN > 0) {
            return uDotN * phiLeft;
        } else {
            return uDotN * phiRight;
        }
    }
}

double FVSolver::computeDiffusiveFlux(int faceId, double phiLeft, double phiRight, 
                                     const std::vector<double>& gradLeft, 
                                     const std::vector<double>& gradRight) {
    const Face& face = geometry->getFaces()[faceId];
    Point normal = geometry->getFaceNormal(faceId);
    
    // 计算梯度在法向量方向的分量
    double gradLeftDotN = gradLeft[0] * normal.x + gradLeft[1] * normal.y;
    double gradRightDotN = gradRight[0] * normal.x + gradRight[1] * normal.y;
    
    // 使用中心差分
    return -params.diffusionCoefficient * 0.5 * (gradLeftDotN + gradRightDotN);
}

void FVSolver::writeTecplotOutput(const std::string& filename, int timeStep) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot create output file " << filename << std::endl;
        return;
    }
    
    file << "TITLE = \"Finite Volume Solution\"" << std::endl;
    file << "VARIABLES = \"X\", \"Y\", \"PHI\"" << std::endl;
    file << "ZONE T=\"Time " << timeStep << "\", N=" << geometry->getNumNodes() 
         << ", E=" << geometry->getNumCells() << ", F=FEPOINT, ET=TRIANGLE" << std::endl;
    
    // Output nodes with solution values
    for (int i = 0; i < geometry->getNumNodes(); i++) {
        const Point& node = geometry->getNodes()[i];
        
        // Find the cell containing this node and get its solution value
        double phi = 0.0;
        int count = 0;
        for (int j = 0; j < geometry->getNumCells(); j++) {
            const Cell& cell = geometry->getCells()[j];
            for (int nodeId : cell.nodes) {
                if (nodeId == i) {
                    phi += solution[j];
                    count++;
                    break;
                }
            }
        }
        if (count > 0) phi /= count;
        
        file << std::fixed << std::setprecision(6) 
             << node.x << " " << node.y << " " << phi << std::endl;
    }
    
    // Output cell connectivity
    for (const auto& cell : geometry->getCells()) {
        if (cell.type == 3) {  // Triangle
            file << cell.nodes[0] + 1 << " " << cell.nodes[1] + 1 << " " << cell.nodes[2] + 1 << std::endl;
        } else if (cell.type == 4) {  // Quadrilateral, decompose into two triangles
            file << cell.nodes[0] + 1 << " " << cell.nodes[1] + 1 << " " << cell.nodes[2] + 1 << std::endl;
            file << cell.nodes[0] + 1 << " " << cell.nodes[2] + 1 << " " << cell.nodes[3] + 1 << std::endl;
        }
    }
    
    file.close();
    std::cout << "Output file saved: " << filename << std::endl;
}

void FVSolver::setInitialCondition(const std::vector<double>& initialValues) {
    if (initialValues.size() == solution.size()) {
        solution = initialValues;
    } else {
        std::cerr << "Error: Initial condition size mismatch" << std::endl;
    }
}

void FVSolver::setInitialCondition(double value) {
    std::fill(solution.begin(), solution.end(), value);
}
