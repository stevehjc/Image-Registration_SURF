#define _CRT_SECURE_NO_WARNINGS

#include<stdio.h>
#include<iostream>
#include"proc.h"
#include"mSURF.h"

using namespace std;

int main()
{
	int rows, cols;

	unsigned char* srcImg1_ = ReadBMPFile(&rows, &cols, "3.bmp");
	unsigned char* srcImg2_ = ReadBMPFile(&rows, &cols, "4.bmp");

	double * srcImg1 = new double[rows*cols]();
	double * srcImg2 = new double[rows*cols]();
	for (int i = 0; i < rows*cols; i++)
	{
		srcImg1[i] = srcImg1_[i];
		srcImg2[i] = srcImg2_[i];
	}
	const int Rows = rows;
	const int Cols = cols;

	//1.特征点检测和匹配
	vector<Point2D> matchpoints1;
	vector<Point2D> matchpoints2;
	int matches_num = 150;
	bool isFeatureDection =
		mSURF_Detection(srcImg1, srcImg2, Rows, Cols, matchpoints1, matchpoints2, matches_num);
	//2.变换模型选择与重映射
	vector<double> transParameter;
	transParameter = getRANSAC2(matchpoints1, matchpoints2, 9, 0);//【*透视矩阵】
	if (!transParameter.size())
	{
		printf("透视变换矩阵求解错误！");
		return 0;
	}
	printf("变换矩阵参数：\n");
	for (int i = 0; i < transParameter.size(); i++)
	{
		printf("%lf ", transParameter[i]);		
		if ((i + 1) % 3 == 0)
		{
			printf("\n");
		}
	}
	printf("1\n");
	double * dstImg1_adj = new double[rows*cols]();
	mWarpPerspective(srcImg1, dstImg1_adj, transParameter, Rows, Cols);
	WriteBMPFile(dstImg1_adj, Rows, Cols, "registrated.bmp");

	return 0;

}