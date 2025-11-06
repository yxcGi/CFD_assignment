#pragma once

#include "Point3D.hpp"
#include "CElements.hpp"
#include "CFaces.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <memory>
#include <cstring>

using namespace std;

class CMesh
{
public:
	//网格的单元总数
	size_t element_number;
	//单元的类型，0-mixed，1-三角形，2-四面体，3-四边形，4-六面体，5-金字塔，6-楔形，7-多面体。	
	int type;
	//网格的面总数
	int face_number;
	std::vector<CElements> CellList; //单元list
	std::vector<CElements> NodeList; //单元list
	std::vector<int> owner;
	std::vector<int> neighbour;
	std::vector<CFaces> FaceList; //面list
	std::vector<Vertex> PointList; //点List
	std::vector<Vector> vlist;      //两个点组成的向量，即面向量
	std::vector<Scalar> AreaList; //面积list
	std::vector<Vector> NorList; //法向量List

	std::vector<Scalar> CellVolumnList;//体积List

	std::vector<Vertex> FaceGeometryCenterList; //面几何中心
	std::vector<Vertex> FaceCentroidList; //面质心

	std::vector<Vertex> CellGeometryCenterList; //单元几何中心
	std::vector<Vertex> CellCentroidList; //单元质心
	vector<vector<int>> PointInFaceList;

public:
	void ReadOwner(const char * filename);
	void ReadNeighbour(const char * filename);
	void ReadFaces(const char * filename);
	void ReadPoints(const char * filename);
	void ComputeArea();
	void ComputeGeometryCenter();
	void ComputeCentriod();
	void RemoveChar(std::string & ss);
	void RemoveSpace(std::string & ss);
	void WritePltFile(const char * filename);
	vector<int> FindFaces(int index);
	void FindPointFromFace();
	void ComputeVolumn();
	void ReadFoamMesh(const char * filename);
};

/*
读取owner文件
*/
void CMesh::ReadOwner(const char * filename) {

	char str[256] = { 0 };
	char temp[256] = { 0 };
	ifstream file;
	bool rc = false;
	int ncount = 0;
	file.open(filename, ios::in);
	int Max = 0;

	if (!file)
		return;
	while (true)
	{
		int num = 0;

		if (file.eof())
			break;
		file.getline(str, 256);
		if (strcmp(str, ")") == 0) {
			break;
		}
		if (rc) {
			num = atoi(str);
			owner.push_back(num);
			Max = (Max > num) ? Max : num;
			CellList[num].faces.push_back(ncount);
			ncount++;
		}

		if (strcmp(str, "(") == 0) {
			num = atoi(temp);
			CellList.resize(num);
			cout << temp << endl;
			rc = true;
		}

		strcpy(temp, str);
	};

	this->element_number = Max + 1;
	CellList.resize(element_number);
	file.close();

	return;
}
/*
读取neighbour文件
*/
void CMesh::ReadNeighbour(const char * filename)
{
	char str[256] = { 0 };
	char temp[256] = { 0 };
	ifstream file;
	bool rc = false;
	int ncount = 0;

	file.open(filename, ios::in);

	if (!file)
		return;
	while (true)
	{
		int num = 0;

		if (file.eof())
			break;
		file.getline(str, 256);

		if (strcmp(str, ")") == 0) {
			break;
		}
		if (rc) {
			num = atoi(str);

			neighbour.push_back(num);
			CellList[num].faces.push_back(ncount);
			ncount++;
		}

		if (strcmp(str, "(") == 0) {
			cout << temp << endl;
			rc = true;
		}

		strcpy(temp, str);
	};


	file.close();

	return;

}
/*
读取faces文件
*/
void CMesh::ReadFaces(const char * filename) {

	char str[256] = { 0 };
	char temp[256] = { 0 };
	ifstream file;
	bool rc = false;
	int ncount = 0;
	file.open(filename, ios::in);

	while (true)
	{
		int num = 0;

		if (file.eof())
			break;
		file.getline(str, 256);

		if (strcmp(str, ")") == 0) {
			break;
		}
		if (rc) {
			char *s;
			char * s2 = str;
			char * s3;
			CFaces f;

			strtok_r(str, "(", &s);
			num = atoi(str);
			f.NumOf = num;
			f.index = ncount;

			s3 = strtok_r(s2 + 2, " ", &s);
			num = atoi(s3);
			f.point_index.push_back(num);

			while (*s != '\0') {
				s3 = strtok_r(NULL, " ", &s);
				num = atoi(s3);
				f.point_index.push_back(num);
			}

			FaceList.push_back(f);
			ncount++;
		}
		//cout << str << endl;
		if (strcmp(str, "(") == 0) {
			cout << temp << endl;
			rc = true;
		}

		strcpy(temp, str);
	};

	file.close();
}
/*
读取points文件
*/
void CMesh::ReadPoints(const char * filename) {
	char str[256] = { 0 };
	char temp[256] = { 0 };
	
	ifstream file;
	bool rc = false;
	int ncount = 0;

	file.open(filename, ios::in);

	while (true)
	{
		int num = 0;
		Point3D<Scalar> pnt;
		vector<double> p2;
		if (file.eof())
			break;
		file.getline(str, 256);

		if (strcmp(str, ")") == 0) {
			break;
		}
		if (rc) {
			char *s;
			char * s3;
			Vertex f;

			s3 = strtok_r(str + 1, " ", &s);
			double px = atof(s3);
			p2.push_back(px);

			while (*s != '\0') {
				s3 = strtok_r(NULL, " ", &s);
				px = atof(s3);
				p2.push_back(px);
			}
			pnt.px = p2[0];
			pnt.py = p2[1];
			pnt.pz = p2[2];
			f = pnt;
		

			PointList.push_back(f);
			ncount++;
		}

		if (strcmp(str, "(") == 0) {
			cout << temp << endl;
			rc = true;
		}

		strcpy(temp, str);
	};

	file.close();
}


/*
计算所有面的面积
*/
void CMesh::ComputeArea()
{
	size_t len = FaceList.size();
	std::vector<CFaces> List = FaceList;
	std::vector<Vertex> PList = PointList;
	size_t i, j;
	
	for (i = 0; i < len; i++) {
		size_t len2 = List[i].point_index.size();
		//
		std::vector<int> point_index = List[i].point_index;
		for (j = 0; j < 2; j++) {
			int m = point_index[j];
			Vertex pnt = PList[m];
			int n = point_index[j + 1];
			Vertex pnt2 = PList[n];
			vlist.push_back(pnt- pnt2);
		}
		//求面的单位法向量，先叉乘，再除以模
		Vector v1 = vlist[i * 2].crossby(vlist[i * 2 + 1]);
		//叉乘的模就是此面的面积，（由于是四边形，如果是三角形则mag/2）
		Scalar mag = v1.GetMag();
		AreaList.push_back(mag);
		Vector normal = v1 / mag;
		NorList.push_back(v1);
	}


}
/*
这个函数计算所有面的几何中心
和单元的几何中心
*/
void CMesh::ComputeGeometryCenter()
{
	size_t len = FaceList.size();
	size_t i, j, k;
	CellCentroidList.resize(element_number);
	//所有面
	for (i = 0; i < len; i++) {
		//面包含点的个数
		size_t p_len = FaceList[i].point_index.size();
		//面包含的点
		std::vector<int> index = FaceList[i].point_index;

		Scalar px = 0;
		Scalar py = 0;
		Scalar pz = 0;
		//计算面的几何中心，即xi = (Σx1-xn)/n
		for (j = 0; j < p_len; j++) {
			int m = index[j];
			Vertex pnt = PointList[m];
			px += pnt.px;
			py += pnt.py;
			pz += pnt.pz;
		}
		//仿真List中
		Vertex point(px/ p_len, py / p_len, pz / p_len);
		FaceGeometryCenterList.push_back(point);
	}
	//所有单元
	for (i = 0; i < element_number; i++) {
		//遍历每个单元中包含的面,正六面体，有6个面
		std::vector<int> faceindex = CellList[i].faces;
		Scalar px = 0;
		Scalar py = 0;
		Scalar pz = 0;
		int count = 0;
		for (j = 0; j < faceindex.size(); j++) {
			//在遍历6个面
			int index = faceindex[j];
			//FaceList是面集合
			std::vector<int> pindex = FaceList[index].point_index;
			size_t p_len = pindex.size();
			//计算单元几何中心（所有点/count）
			for (k = 0; k < p_len; k++) {
				int m = pindex[k];
				//PointList是点集合
				Vertex pnt = PointList[m];
				px += pnt.px;
				py += pnt.py;
				pz += pnt.pz;
				count++;
			}
		}
		Vertex point(px / count, py / count, pz / count);
		CellGeometryCenterList.push_back(point);
		CellCentroidList[i] = point;
	}
}

void CMesh::RemoveChar(std::string & ss)
{
	ss.erase(std::remove_if(ss.begin(), ss.end(), 
		[](char ch) {
		if (ch >= 0x30 && ch <= 0x39)
			return false;
		return true;
	}), ss.end());
}

void CMesh::RemoveSpace(std::string & ss)
{
	ss.erase(std::remove_if(ss.begin(), ss.end(),
		[](char ch) {
		if (ch == ' ' || ch  == '\t' || ch == '\n')
			return true;
		return false;
	}), ss.end());
}

void CMesh::ComputeCentriod() {

	//得到面和点的List
	size_t len = FaceList.size();
	std::vector<CFaces> List = FaceList;
	std::vector<Vertex> PList = PointList;
	size_t i, j;
	//质心的数组大小
	FaceCentroidList.resize(len);

	for (i = 0; i < len; i++) {
		//每个面是由几个点组成的
		size_t len2 = List[i].point_index.size();
		//面的几何中心点
		Vertex cpnt = FaceGeometryCenterList[i];
		//点的指针
		std::vector<int> point_index = List[i].point_index;
		//面积等于所有小三角形面积之和
		Scalar AreaSum = 0;
		//面的法向量
		Vector normal;
		//如果是三角形，那么只有3个点，直接计算
		//add 计算质心
		if(len2 == 3){
			Vector a = PList[point_index[0]];
			Vector b = PList[point_index[1]];
			Vector c = PList[point_index[2]];
			a = b - a;
			c = b - c;
			Vector varea = a.crossby(c)/2;
			AreaSum = varea.GetMag();
			normal = varea / varea.GetMag();
			//对于三角形，几何中心就是质心
			FaceCentroidList[i] = FaceGeometryCenterList[i];
		}else{
			//下面是求多边形质心的通用程序
			Scalar sumPx =0;
			Scalar sumPy =0;
			Scalar sumPz =0;
			for (j = 0; j < len2; j++) {
				int m = point_index[j];
				int n = 0;
				if (j == len2 - 1) {
					n = point_index[0];
				}
				else {
					n = point_index[j + 1];
				}
				Vertex pnt = PList[m];
				Vertex pnt2 = PList[n];
	
				Vector v1 = pnt - cpnt;
				Vector v2 = pnt2 - cpnt;
				//求小三角形的质心
				Scalar px = (pnt.px + pnt2.px + cpnt.px)/3;
				Scalar py = (pnt.py + pnt2.py + cpnt.py)/3;
				Scalar pz = (pnt.pz + pnt2.pz + cpnt.pz)/3;
				//求小三角形面积
				Vector varea = v1.crossby(v2)/2;
				//求模
				Scalar Mag = varea.GetMag();
				sumPx += px * Mag;
				sumPy += py * Mag;
				sumPz += pz * Mag;
				//面积求和
				AreaSum += Mag;
				//得到法向量
				normal = varea / Mag;
			}
			Scalar ctpx = sumPx / AreaSum;
			Scalar ctpy = sumPy / AreaSum;
			Scalar ctpz = sumPz / AreaSum;
			Vector ctPnt(ctpx,ctpy,ctpz);
			FaceCentroidList[i] = ctPnt;
		}
		//将面积和面积矢量压入List中
		AreaList.push_back(AreaSum);
		normal = normal * AreaSum;
		NorList.push_back(normal);
	}
	//求多面体的体心（质心）和体积
	//所有单元(有bug，暂不用)
	
	for (i = 0; i < element_number; i++) {
		//遍历每个单元中包含的面,比如正六面体，有6个面
		std::vector<int> faceindex = CellList[i].faces;
		int count = 0;
		Vector Xcepy(0,0,0);
		Vector XSum(0,0,0);
		Scalar LocalVolumnSum = 0;
		for (j = 0; j < faceindex.size(); j++) {
			//在遍历6个面，求每个金字塔的质心
			int index = faceindex[j];
			Xcepy = 0.75 * FaceCentroidList[index] + 0.25 * CellGeometryCenterList[i];
			Vector Cf = FaceCentroidList[index] - CellGeometryCenterList[i];
			Vector Sf = NorList[index];
			Scalar LocalVolumn = (Cf.dotby(Sf))/3;
			XSum = XSum + Xcepy * LocalVolumn;
			LocalVolumnSum += LocalVolumn;
			
		}
		//计算的体积和质心保存
		CellVolumnList.push_back(LocalVolumnSum);
		XSum = XSum / LocalVolumnSum;
		CellCentroidList.push_back(XSum);
	}

}


void str_replace(char *str, char old_char1, char old_char2, char new_char) {
	int i = 0;
	while (str[i] != '\0') {
		if (str[i] == old_char1 || str[i] == old_char2) {
			str[i] = new_char;
		}
		i++;
	}
}

vector<string> split(const string& s, char delimiter) {
	vector<string> tokens;
	string token;
	for (char c : s) {
		if (c == delimiter) {
			if (token.size() != 0)
				tokens.push_back(token);
			token.clear();
		}
		else {
			token += c;
		}
	}
	tokens.push_back(token); // 添加最后一个token
	return tokens;
}

inline void CMesh::WritePltFile(const char * filename)
{
	char str[256] = { 0 };
	fstream file;
	file.open(filename, ios::out);
	//写入变量
	const char * str1 = "Variables = x, y, T\n";
	file.write(str1, strlen(str1));
	stringstream fmt1;
	fmt1 << "Zone N=" << PointList.size() << ", E =" << element_number <<
		",datapacking=point,varlocation=([" << element_number << "]=cellcentered),zonetype=fetriangle\n";
	//写入节点，元素等信息
	string targetString = fmt1.str();
	file.write(targetString.c_str(), targetString.size());//写入文件
	//写入点
	for (size_t i = 0; i < PointList.size(); i++) {
		stringstream fmt;
		fmt << PointList[i].px << " "
			<< PointList[i].py << " "
			<< 0 << std::endl;
		string targetString = fmt.str();
		file.write(targetString.c_str(), targetString.size());//写入文件
	}
	//写入单元信息，每个单元由哪几个点组成
	for (size_t i = 0; i < CellList.size(); i++) {
		stringstream fmt;
			fmt << CellList[i].faces[0] + 1 << " "
				<< CellList[i].faces[1] + 1 << " "
				<< CellList[i].faces[2] + 1 << std::endl;
		string targetString = fmt.str();
		file.write(targetString.c_str(), targetString.size());//写入文件
	}
	file.close();
}

inline vector<int> CMesh::FindFaces(int index)
{
	vector<int> face;
	for (unsigned int i = 0; i < FaceList.size(); i++) {
		//查找点
		for (unsigned int j = 0; j < FaceList[i].point_index.size(); j++) {
			if (index == FaceList[i].point_index[j])
			{
				face.push_back(i);
				break;
			}
		}
	}
	return face;
}

inline void CMesh::FindPointFromFace()
{
	//PointInFaceList.resize(PointList.size());
	//遍历所有节点
	for (size_t i = 0; i < PointList.size(); i++) {
		//查找围绕这个节点的所有单元和面
		PointInFaceList.push_back(FindFaces(i));
	}
}

inline void CMesh::ComputeVolumn()
{
	for (size_t i = 0; i < element_number; i++)
	{
		//三个点直接求面积
		Vertex p1 = PointList[NodeList[i].nodes[0]];
		Vertex p2 = PointList[NodeList[i].nodes[1]];
		Vertex p3 = PointList[NodeList[i].nodes[2]];

		Vector v1 = p2 - p1;
		Vector v2 = p3 - p2;

		v1 = v1.crossby(v2);
		CellVolumnList.push_back(v1.GetMag());
	}
}

inline void CMesh::ReadFoamMesh(const char * filename)
{
	wstring filepath;
	filepath = filename;
	filepath.append(L"/owner");
	ReadOwner(filepath.c_str());
	filepath = filename;
	filepath.append(L"/neighbour");
	ReadNeighbour(filepath.c_str());
	filepath = filename;
	filepath.append(L"/faces");
	ReadFaces(filepath.c_str());
	filepath = filename;
	filepath.append(L"/points");
	ReadPoints(filepath.c_str());

	ComputeGeometryCenter();
	ComputeCentriod();
}


