#include "BoundaryPatch.h"

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