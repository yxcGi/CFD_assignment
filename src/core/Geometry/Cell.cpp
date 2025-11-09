#include "Cell.h"

Cell::Cell()
{
    faceIndexes_.reserve(4);
}

const std::vector<Cell::ULL>& Cell::getFaceIndices() const
{
    return faceIndexes_;
}

Scalar Cell::getVolume() const
{
    return volume_;
}

const Vector<Scalar>& Cell::getCenter() const
{
    return center_;
}

int Cell::getFaceNum() const
{
    return faceIndexes_.size();
}

void Cell::addFaceIndex(ULL faceIndex)
{
    faceIndexes_.emplace_back(faceIndex);
}

void Cell::calculateCellInfo()
{

}
