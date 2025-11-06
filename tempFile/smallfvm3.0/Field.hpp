#pragma once

#include "CMesh.hpp"
#include "Point3D.hpp"


template<class Type>
class Field
{
public:
	//单元的中心值（质心），也是要求解的值
	std::vector<Type> phi;
	CMesh * msh;
public:
	Field() {}
	Field(CMesh * mesh, Type value);
};

template<class Type>
inline Field<Type>::Field(CMesh * mesh, Type value) {
	//定义phi的size和元素大小是一样
	phi.resize(mesh->element_number);
	for (size_t i = 0;i < phi.size();i++)
	{
		phi[i] = value;
	}
	msh = mesh;
}
