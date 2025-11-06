#ifndef FV_SOLVER_H
#define FV_SOLVER_H

#include "geometry.h"
#include <vector>
#include <string>
#include <map>

// 求解参数结构
struct SolverParameters {
    double convectionVelocityX;
    double convectionVelocityY;
    double diffusionCoefficient;
    double sourceTerm;
    double timeStep;
    double finalTime;
    int maxIterations;
    double tolerance;
    std::string reconstructionMethod;  // "gauss" 或 "least_squares"
    std::string fluxScheme;  // "upwind" 或 "central"
};

// 边界条件参数
struct BoundaryCondition {
    std::string type;  // "dirichlet", "neumann", "wall"
    double value;
};

class FVSolver {
private:
    Geometry* geometry;
    SolverParameters params;
    std::map<int, BoundaryCondition> boundaryConditions;
    
    // 求解变量
    std::vector<double> solution;  // 当前解
    std::vector<double> solutionOld;  // 上一时间步解
    std::vector<std::vector<double>> gradients;  // 梯度 [cellId][component]
    
    // 内部函数
    void initializeSolution();
    void computeGradients();
    void computeGaussGradients();
    void computeLeastSquaresGradients();
    void reconstructSolution();
    double computeConvectiveFlux(int faceId, double phiLeft, double phiRight);
    double computeDiffusiveFlux(int faceId, double phiLeft, double phiRight, 
                               const std::vector<double>& gradLeft, 
                               const std::vector<double>& gradRight);
    void applyBoundaryConditions(std::vector<double>& residual);
    void updateSolution();
    double computeResidual();
    void writeTecplotOutput(const std::string& filename, int timeStep);
    
    // 隐式方法相关函数
    void buildImplicitSystem(std::vector<std::vector<double>>& A, std::vector<double>& b);
    void buildInternalFaceJacobian(int faceId, std::vector<std::vector<double>>& A, std::vector<double>& b);
    void buildBoundaryFaceJacobian(int faceId, std::vector<std::vector<double>>& A, std::vector<double>& b);
    void solveLinearSystem(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x);
    
public:
    FVSolver(Geometry* geom);
    ~FVSolver();
    
    // 主要接口
    void loadParameters(const std::string& filename);
    void setBoundaryCondition(int boundaryId, const std::string& type, double value);
    void solve();
    
    // 获取器
    const std::vector<double>& getSolution() const { return solution; }
    double getSolutionAtCell(int cellId) const { return solution[cellId]; }
    
    // 设置器
    void setInitialCondition(const std::vector<double>& initialValues);
    void setInitialCondition(double value);
};

#endif // FV_SOLVER_H
