#ifndef _MATH_M_H_0123456_
#define _MATH_M_H_0123456_


#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

#include<time.h>

using namespace std;
const int MAX_ITER = 100000;
//const double eps = 0.0000001;
const double eps = 2.2204e-15;
const double PI = 3.1415926535897932384626433832795;

double get_norm(double *x, int n);
double normalize(double *x, int n);
inline double product(double*a, double *b, int n);
void orth(double *a, double *b, int n);
void print(vector<vector<double> > &A);

bool svd(vector<vector<double> > A, int K,
	std::vector<std::vector<double> > &U,
	std::vector<double> &S, std::vector<std::vector<double> > &V);

bool Inv_svd(vector<vector<double> > A, vector<vector<double> > &invA, int n);
bool Inv_svdMN(vector<vector<double> > A, vector<vector<double> > &invA, int n);
bool MatrixMultiplyV(vector<vector<double> > X1, vector<vector<double> > X2, vector<vector<double> > &Y);

/*矩阵运算函数*/
bool InvertGaussJordan(double *m_pData, int m_nNumColumns);
void swap(double *a, double *b);

bool Inv_matrix(double *p, int n);

double *MatrixMultiply(double *L_Data, int Lrow, int Lcol, double *R_Data, int Rcol);

//方差计算
double Deviation(double *in_Data, int row, int col);
double Deviation(double *in_Data, int col, int row_begin, int row_count, int col_begin, int col_count);
//均值计算
double MeanImg(double *in_Data, int Rows, int Cols);

/*计算相关性系数*/
double correlation_coefficient(double *Bx, double *By, int N);
/*计算绝对差*/
double SAD(double *Bx, double *By, int N);

/*最大类间差方法*/
double getOstu(double *in_Data, int row, int col);
double get1DMaxEntropyThreshold(double *in_Data, int Rows, int Cols);
double getGuDi(double *in_Data, int Rows, int Cols);
bool IsDimodal(double* HistGram);      // 检测直方图是否为双峰的


/*得到四个随机数*/
void getRandom(int a[4], int max);

double* mgetGaussianKernel(int n, double sigma);





#endif