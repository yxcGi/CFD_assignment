#ifndef BOUNDARYPATCH_H_
#define BOUNDARYPATCH_H_

#include <string>

class BoundaryPatch
{
    using ULL = unsigned long long;
public:

private:
    std::string name_;
    std::string type_;
    ULL nFaces_;
    ULL startFace_;
};

#endif // BOUNDARYPATCH_H_