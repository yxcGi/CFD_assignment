#ifndef BOUNDARYPATCH_H_
#define BOUNDARYPATCH_H_

#include <string>
#include <iostream>

class BoundaryPatch
{
    using Scalar = double;
    using ULL = unsigned long long;
public:
    enum class BoundaryType
    {
        PATCH,      // 净出口
        WALL,       // 无滑移壁面
        SYMMETRY,   // 对称面
        CYCLIC,     // 周期性边界
        WEDGE,      // 楔形边界
        EMPTY,      // 空边界
        PROCESSOR,  // 处理器边界，并行
        CUSTOM      // 自定义边界类型，未设置
    };
    
public:
    BoundaryPatch() = delete;
    BoundaryPatch(
        const std::string& name,
        ULL nFaces,
        ULL startFace,
        BoundaryType type = BoundaryType::CUSTOM
    );
    BoundaryPatch(const BoundaryPatch&) = delete;
    BoundaryPatch(BoundaryPatch&&) = default;
    BoundaryPatch& operator=(const BoundaryPatch&) = delete;
    BoundaryPatch& operator=(BoundaryPatch&&) = default;

    ~BoundaryPatch() = default;

public:
    // 常用接口
    std::string getName() const;
    BoundaryType getType() const;
    ULL getNFace() const;
    ULL getStartFace() const;
    void setBoundaryType(BoundaryType type);




private:
    std::string name_;      // 边界名称
    ULL nFaces_;            // 面的数量
    ULL startFace_;         // 起始面索引
    BoundaryType type_;     // 边界类型
    bool isValid_;          // 边界是否有效，默认未设置就是无效
};

#endif // BOUNDARYPATCH_H_