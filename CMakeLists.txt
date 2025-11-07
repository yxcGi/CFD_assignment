# SIMPLE 算法求解顶盖驱动流
cmake_minimum_required(VERSION 3.15)

get_filename_component(PRJ_NAME ${CMAKE_SOURCE_DIR} NAME)

project(${PRJ_NAME}
    VERSION 1.0.0
    DESCRIPTION "The SIMPLE arithmetic to solver a lid lid driven cavity flow"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 生成索引文件
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# 编译选项
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wextra -pedantic")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
endif()


# 设置可执行文件的输出路径
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
# 设置动态库输出路径
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
# 设置静态库输出路径
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# 添加子目录
add_subdirectory(src)
