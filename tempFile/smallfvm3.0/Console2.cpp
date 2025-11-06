// Console2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <windows.h>
#include "G:/CFD-PRO/smallfvm/CMesh.hpp"
#include "G:/CFD-PRO/smallfvm/Field.hpp"

vector<double> Gauss_Seidel(vector<vector<double>> &A, vector<double> & B, size_t it);

int main()
{
	CMesh * msh = new CMesh;
	msh->ReadFoamMesh(L"D:/1/mesh");

	Field<Scalar> phiT(msh, 0);
	double gamma = 0.1;

	//Ax = b
	vector<vector<Scalar>> A(msh->element_number);
	vector<Scalar> b(msh->element_number, 0);

	int i = 0;
	for (i = 0; i < msh->element_number; i++)
		A[i].resize(msh->element_number, 0);
	int len = msh->neighbour.size();
	for (i = 0; i < len; i++) {
		//内部面的主单元和邻单元号，面的序号是i
		int cellC = msh->owner[i];
		int cellF = msh->neighbour[i];
		//然后得到单元的质心（结构化网格，几何中心和质心重合）
		Vertex cPoint = msh->CellGeometryCenterList[cellC];
		Vertex fPoint = msh->CellGeometryCenterList[cellF];
		Vertex faceCtPoint = msh->FaceCentroidList[i];
		//
		Scalar FaceArea = msh->AreaList[i];
		Scalar Distance = cPoint.GetDistance(fPoint);

		Scalar dSum = -gamma * FaceArea / Distance;

		A[cellC][cellC] += dSum;
		A[cellC][cellF] -= dSum;
		A[cellF][cellF] += dSum;
		A[cellF][cellC] -= dSum;
	}
	for (i = 180; i < 190; i++) {
		int cellC = msh->owner[i];
		//得到点
		Vertex cPoint = msh->CellGeometryCenterList[cellC];
		Vertex facePoint = msh->FaceGeometryCenterList[i];

		//面心到单元中心的距离
		Scalar dcf = cPoint.GetDistance(facePoint);
		Scalar Area = msh->AreaList[i];

		Scalar ab = -gamma * Area / dcf;

		A[cellC][cellC] += ab;
		b[cellC] = ab * 0;
	}

	for (i = 200; i < 210; i++) {
		int cellC = msh->owner[i];
		//得到点
		Vertex cPoint = msh->CellGeometryCenterList[cellC];
		Vertex facePoint = msh->FaceGeometryCenterList[i];

		//面心到单元中心的距离
		Scalar dcf = cPoint.GetDistance(facePoint);
		Scalar Area = msh->AreaList[i];

		Scalar ab = -gamma * Area / dcf;

		A[cellC][cellC] += ab;

		b[cellC] = ab * 300;
	}

	phiT.phi = Gauss_Seidel(A, b, 100);

	for (i = 0; i < msh->element_number; i++)
		cout << phiT.phi[i] << '\t';
	if (i % 9 == 0)
		cout << endl;
	cout << endl;
}

//高斯-赛戴尔迭代法函数
vector<double> Gauss_Seidel(vector<vector<double>> &A, vector<double> & B, size_t it) //高斯-赛戴尔迭代法函数
{
	size_t i, j, k;
	size_t n = B.size();
	double temp = 0;
	//要求的解
	vector<double> X(n);
	//误差
	double err = 0;
	//方法主要部分
	for (k = 0; k < it; k++) { //最大迭代次数为1000
		err = 0;
		for (i = 0; i < n; i++) { //双重for循环遍历数组
			double sum = 0;
			for (j = 0; j < n; j++) {
				if (A[i][j] == 0)
					continue;
				sum += A[i][j] * X[j]; // 对除对角线的其余元素求和
			}
			temp = B[i] - sum;
			X[i] = X[i] + temp / A[i][i]; //计算x
			err += abs(temp); //计算误差
		}
		std::cout << "迭代次数：" << k << ";\t残差为: " << err << endl;

		if (err < 1e-7)
		{
			printf("求解已收敛\n");
			break;
		}

	}
	return X;
}