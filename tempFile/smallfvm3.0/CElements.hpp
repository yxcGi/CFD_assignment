#pragma once

#include <vector>
#include "CFaces.hpp"

class CElements
{
public:
	//owner和neighbour文件中的面
	std::vector<int> faces;
	//单元中的节点
	std::vector<int> nodes;
};

