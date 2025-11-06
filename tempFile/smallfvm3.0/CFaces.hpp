#pragma once

#include <vector>
//#include "CPoints.hpp"

class CFaces
{
public:
	//这个面是由几个点组成
	int NumOf;
	//面的序号
	int index;
	std::vector<int> point_index;
public:
	CFaces() {}
	CFaces(int n,int i, std::vector<int> p) {
		NumOf = n;
		index = i;
		point_index = p;
	}
};

