#include "BoundaryPatch.h"
#include <iostream>

using Scalar = double;
using ULL = unsigned long long;

BoundaryPatch::BoundaryPatch(
    const std::string& name,
    ULL nFaces,
    ULL startFace,
    BoundaryType type
)
    : name_(name)
    , nFaces_(nFaces)
    , startFace_(startFace)
    , type_(type)
    , isValid_(false)
{
    if (type != BoundaryType::CUSTOM)
    {
        isValid_ = true;
    }
}

std::string BoundaryPatch::getName() const
{
    return name_;
}

BoundaryPatch::BoundaryType BoundaryPatch::getType() const
{
    return type_;
}

ULL BoundaryPatch::getNFace() const
{
    return nFaces_;
}

ULL BoundaryPatch::getStartFace() const
{
    return startFace_;
}

void BoundaryPatch::setBoundaryType(BoundaryType type)
{
    if (type == BoundaryType::CUSTOM)
    {
        isValid_ = false;
        std::cout << "[Warning] Setting boundary type to CUSTOM makes the boundary invalid." << std::endl;
        return;
    }
    type_ = type;
    isValid_ = true;
}
