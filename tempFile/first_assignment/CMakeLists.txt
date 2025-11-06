cmake_minimum_required(VERSION 3.10)
project(FVSolver)

# 设置C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置编译选项
if(MSVC)
    # Windows MSVC编译器选项
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /O2")
else()
    # Linux/Mac GCC/Clang编译器选项
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O2")
endif()

# 添加可执行文件
add_executable(fv_solver
    main.cpp
    geometry.cpp
    fv_solver.cpp
)

# 设置输出目录
set_target_properties(fv_solver PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# 安装目标
install(TARGETS fv_solver DESTINATION bin)
