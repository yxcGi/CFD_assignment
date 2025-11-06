#include "geometry.h"
#include "fv_solver.h"
#include <iostream>
#include <string>
#include <algorithm>

int main() {
    std::cout << "=== 2D Unstructured Grid Finite Volume Solver ===" << std::endl;
    std::cout << "Solving general convection-diffusion-source equation" << std::endl;
    std::cout << "Supporting GMSH triangular and quadrilateral mixed grids" << std::endl;
    std::cout << "Second-order accuracy with Gauss and least squares reconstruction" << std::endl;
    
    try {
        // Create geometry object
        Geometry geometry;
        
        // Load mesh
        std::string meshFile = "complex_mesh.msh";
        std::cout << "\nLoading mesh file: " << meshFile << std::endl;
        geometry.loadMesh(meshFile);
        
        // Create solver
        FVSolver solver(&geometry);
        
        // Load solver parameters
        std::string paramFile = "parameters.txt";
        std::cout << "\nLoading parameter file: " << paramFile << std::endl;
        solver.loadParameters(paramFile);
        
        // Set boundary conditions
        std::cout << "\nSetting boundary conditions..." << std::endl;
        // Inlet boundary (left boundary) - 复杂进口条件
        solver.setBoundaryCondition(1, "dirichlet", 2.0);
        // Outlet boundary (right boundary) - 自由出口
        solver.setBoundaryCondition(2, "neumann", 0.0);
        // Wall boundaries (top and bottom) - 无滑移壁面
        solver.setBoundaryCondition(3, "wall", 0.0);
        solver.setBoundaryCondition(4, "wall", 0.0);
        
        // Start solving
        solver.solve();
        
        // Output final result statistics
        const std::vector<double>& finalSolution = solver.getSolution();
        double maxValue = *std::max_element(finalSolution.begin(), finalSolution.end());
        double minValue = *std::min_element(finalSolution.begin(), finalSolution.end());
        
        std::cout << "\n=== Solution Statistics ===" << std::endl;
        std::cout << "Maximum concentration: " << maxValue << std::endl;
        std::cout << "Minimum concentration: " << minValue << std::endl;
        
        // Calculate total mass (integral)
        double totalMass = 0.0;
        for (int i = 0; i < geometry.getNumCells(); i++) {
            totalMass += solver.getSolutionAtCell(i) * geometry.getCellArea(i);
        }
        std::cout << "Total mass: " << totalMass << std::endl;
        
        std::cout << "\nSolution completed! Results saved in Tecplot format." << std::endl;
        std::cout << "You can use Tecplot, ParaView or other visualization software to view results." << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
