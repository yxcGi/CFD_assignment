# 添加子文件库
add_subdirectory(core)

# 在CMakeLists.txt中设置
# set(CMAKE_BUILD_TYPE Debug)


# 生成可执行文件
add_executable(${PRJ_NAME} main.cpp)

# 链接动态库
target_link_libraries(${PRJ_NAME} PRIVATE core)

# 添加头文件目录
target_include_directories(
    ${PRJ_NAME}
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/core/
)
