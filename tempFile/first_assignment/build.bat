@echo off
echo === 二维非结构网格有限体积求解器构建脚本 ===

REM 创建构建目录
if not exist build mkdir build
echo 创建构建目录: build

cd build

REM 运行CMake
echo 运行CMake配置...
cmake ..

REM 编译项目
echo 编译项目...
cmake --build . --config Release

if %ERRORLEVEL% EQU 0 (
    echo.
    echo === 构建成功！ ===
    echo 可执行文件位置: build\bin\Release\fv_solver.exe
    echo.
    echo 运行求解器:
    echo cd build\bin\Release ^&^& fv_solver.exe
    echo.
    echo 或者从项目根目录运行:
    echo cd build\bin\Release ^&^& fv_solver.exe
) else (
    echo 构建失败！请检查错误信息。
    exit /b 1
)

pause
