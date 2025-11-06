#ifndef BOUNDARYPATCH_H_
#define BOUNDARYPATCH_H_

#include <string>

class BoundaryPatch
{
    using ULL = unsigned long long;
public:
    enum class BoundaryType
    {
        DIRICHLET,  // 指定值边界
        NEUMANN,    // 指定通量边界
        WALL,       // 无滑移壁面
        SYMMETRY,   // 对称边界
        PERIODIC,   // 周期边界
        CUSTOM      // 自定义边界
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
    BoundaryPatch(BoundaryPatch&&) = delete;
    BoundaryPatch& operator=(const BoundaryPatch&) = delete;
    BoundaryPatch& operator=(BoundaryPatch&&) = delete;

    ~BoundaryPatch() = default;

public:
    // 常用接口
    std::string getName() const;
    BoundaryPatch getType() const;
    ULL getNFace() const;
    ULL getStartFace() const;


private:
    std::string name_;      // 边界名称
    ULL nFaces_;            // 面的数量
    ULL startFace_;         // 起始面索引
    BoundaryType type_;     // 边界类型
    bool isValid_;          // 边界是否有效，默认未设置就是无效
};

#endif // BOUNDARYPATCH_H_