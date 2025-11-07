#include "Cell.h"

Cell::Cell()
{
    faceIndexes_.reserve(4);
}

void Cell::addFaceIndex(ULL faceIndex)
{
    faceIndexes_.emplace_back(faceIndex);
}

void Cell::calculateCellInfo()
{

}
