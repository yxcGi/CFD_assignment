/**
 * @file Point3D.h
 本文件定义了 点，向量的操作
 */
#ifndef _Point3D_
#define _Point3D_

#include <iostream>
#include <math.h>
#include <ostream>

/*标量定义为double类型*/
typedef double Scalar;

/**
 *@点的定义，模板类
 既是一点，还是一个向量
 */
template<class Type>
class Point3D
{
public:
	Type	px, py, pz;
	/**
	 * 2个构造函数
	 * 一个无参，一个3个参数
	 */
	Point3D() :px(0.0), py(0.0), pz(0.0) {}
	Point3D(Type X, Type Y, Type Z) :px(X), py(Y), pz(Z) {}
public:
	/*
	@取模函数
	*/
	Scalar GetMag()
	{
		Scalar ret;
		ret = px * px + py * py + pz * pz;
		ret = sqrt(ret);
		return ret;
	}

	/*===================================================================*/
	/*                                                                   */
	/*            求2个点的距离                                          */
	/*                                                                   */
	/*===================================================================*/
	Scalar GetDistance(Point3D<Type>& v){
		Scalar dx = this->px - v.px;
		Scalar dy = this->py - v.py;
		Scalar dz = this->pz - v.pz;
		
		return sqrt(dx*dx+dy*dy+dz*dz);
	}
	/*
	@输出当前类信息
	*/
	void ToString() {
		std::cout << "( " 
			<< px  << " , " 
			<< py  << " , " 
			<< pz  << " ) "
			<< std::endl;
	}
	/*
	点乘
	*/
	Scalar dotby(Point3D<Type>& b) {
		Scalar ret = 0;
		ret = px * b.px + py * b.py + pz * b.pz;
		return ret;
	}
	// 叉乘
	Point3D<Type>  crossby(Point3D<Type>& v2) {

		return Point3D<Type>
			((this->py * v2.pz - this->pz * v2.py),
			 (this->pz * v2.px - this->px * v2.pz),
			 (this->px * v2.py - this->py * v2.px)
			);
	}
	// <<运算符号重载
	inline friend std::ostream& operator << (std::ostream& out, const Point3D<Type> & rhs)
	{
		out << " " << rhs.px << "  " << rhs.py << "  " << rhs.pz << "  ";
		return out;
	}
	inline friend Point3D<Type> operator * (const Scalar d, const Point3D<Type>& p)
	{
		return Point3D(d * p.px, d * p.py, d * p.pz);
	}
public:
	/*
	加减乘除重载
	*/
	 Point3D<Type> operator + (const Point3D<Type>& v)
	{
		return Point3D(px+ v.px, py+v.py, pz+ v.pz);
	}

	Point3D<Type> operator - (const Point3D<Type>& v)
	{
		return Point3D(px - v.px, py - v.py, pz - v.pz);
	}

	Point3D<Type> operator * (const Scalar d)
	{
		return Point3D(d * px, d * py, d * pz);
	}
	Scalar operator * (const Point3D<Type>& v)
	{
		return (px * v.px + py * v.py + pz * v.pz);
	}
	Point3D<Type> operator / (const Scalar d)
	{
		return Point3D(px /d, py /d, pz/d);
	}
};

typedef Point3D<Scalar> Vector;
typedef Point3D<Scalar> Vertex;

#endif

