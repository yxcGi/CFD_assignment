#include "Face.h"

Face::Face(
    std::vector<ULL>& pointIndexs,
    ULL ownerIndex,
    ULL neighborIndex
)
    : pointNum_(pointIndexs.size())
    , owner_(ownerIndex)
    , neighbor_(neighborIndex)
{
    pointIndexes_ = std::move(pointIndexs);
}

std::ostream& operator<<(std::ostream& out, const Face& face)
{
    out << face.pointNum_ << "(";
    for (auto ele : face.pointIndexes_)
    {
        std::cout << ele << " ";
    }
    std::cout << ")";
    return out;
}


