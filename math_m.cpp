

#include "math_m.h"

double* mgetGaussianKernel(int n, double sigma)
{
	const int SMALL_GAUSSIAN_SIZE = 7;
	static const float small_gaussian_tab[][SMALL_GAUSSIAN_SIZE] =
	{
		{ 1.f },
		{ 0.25f, 0.5f, 0.25f },
		{ 0.0625f, 0.25f, 0.375f, 0.25f, 0.0625f },
		{ 0.03125f, 0.109375f, 0.21875f, 0.28125f, 0.21875f, 0.109375f, 0.03125f }
	};

	const float* fixed_kernel = n % 2 == 1 && n <= SMALL_GAUSSIAN_SIZE && sigma <= 0 ?
		small_gaussian_tab[n >> 1] : 0;


	double* kernel = (double*)calloc(n, sizeof(double));
	//Mat kernel(n, 1, ktype);

	//float* cf = (float*)kernel.data;
	double* cd = kernel;

	double sigmaX = sigma > 0 ? sigma : ((n - 1)*0.5 - 1)*0.3 + 0.8;
	double scale2X = -0.5 / (sigmaX*sigmaX);
	double sum = 0;

	int i;
	for (i = 0; i < n; i++)
	{
		double x = i - (n - 1)*0.5;
		double t = fixed_kernel ? (double)fixed_kernel[i] : std::exp(scale2X*x*x);

		cd[i] = t;
		sum += cd[i];
		//}
	}

	sum = 1. / sum;
	for (i = 0; i < n; i++)
	{

		cd[i] *= sum;
	}

	return kernel;
};



/*计算矩阵相乘*/
//[8][8]x[8][1]
double *MatrixMultiply(double *L_Data, int Lrow, int Lcol, double *R_Data, int Rcol)
{
	double *result = (double *)calloc(Lrow*Rcol, sizeof(double));

	double temp = 0;
	for (int i = 0; i < Lrow; i++)
	{
		temp = 0;
		for (int j = 0; j < Lcol; j++)
		{
			temp += L_Data[i*Lcol + j] * R_Data[j];
		}
		result[i] = temp;
	}
	return result;
}

/*交换指针地址（矩阵求逆中用到）*/
void swap(double *a, double *b)
{
	double c;
	c = *a;
	*a = *b;
	*b = c;
}
//////////////////////////////////////////////////////////////////////
// 实矩阵求逆的全选主元高斯－约当法
//
// 参数：无
//
// 返回值：BOOL型，求逆是否成功
//////////////////////////////////////////////////////////////////////
bool InvertGaussJordan(double *m_pData, int m_nNumColumns)
{
	int *pnRow, *pnCol, i, j, k, l, u, v;
	double d = 0, p = 0;

	// 分配内存
	//pnRow = new int[m_nNumColumns];
	//pnCol = new int[m_nNumColumns];
	pnRow = (int*)calloc(m_nNumColumns, sizeof(int));
	pnCol = (int*)calloc(m_nNumColumns, sizeof(int));
	if (pnRow == NULL || pnCol == NULL)
		return false;

	// 消元
	for (k = 0; k <= m_nNumColumns - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= m_nNumColumns - 1; i++)
		{
			for (j = k; j <= m_nNumColumns - 1; j++)
			{
				l = i*m_nNumColumns + j; p = fabs(m_pData[l]);
				if (p > d)
				{
					d = p;
					pnRow[k] = i;
					pnCol[k] = j;
				}
			}
		}

		// 失败
		if (d == 0.0)
		{
			/*delete[] pnRow;
			delete[] pnCol;*/
			free(pnRow);
			free(pnCol);
			return false;
		}

		if (pnRow[k] != k)
		{
			for (j = 0; j <= m_nNumColumns - 1; j++)
			{
				u = k*m_nNumColumns + j;
				v = pnRow[k] * m_nNumColumns + j;
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}

		if (pnCol[k] != k)
		{
			for (i = 0; i <= m_nNumColumns - 1; i++)
			{
				u = i*m_nNumColumns + k;
				v = i*m_nNumColumns + pnCol[k];
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}

		l = k*m_nNumColumns + k;
		m_pData[l] = 1.0 / m_pData[l];
		for (j = 0; j <= m_nNumColumns - 1; j++)
		{
			if (j != k)
			{
				u = k*m_nNumColumns + j;
				m_pData[u] = m_pData[u] * m_pData[l];
			}
		}

		for (i = 0; i <= m_nNumColumns - 1; i++)
		{
			if (i != k)
			{
				for (j = 0; j <= m_nNumColumns - 1; j++)
				{
					if (j != k)
					{
						u = i*m_nNumColumns + j;
						m_pData[u] = m_pData[u] - m_pData[i*m_nNumColumns + k] * m_pData[k*m_nNumColumns + j];
					}
				}
			}
		}

		for (i = 0; i <= m_nNumColumns - 1; i++)
		{
			if (i != k)
			{
				u = i*m_nNumColumns + k;
				m_pData[u] = -m_pData[u] * m_pData[l];
			}
		}
	}

	// 调整恢复行列次序
	for (k = m_nNumColumns - 1; k >= 0; k--)
	{
		if (pnCol[k] != k)
		{
			for (j = 0; j <= m_nNumColumns - 1; j++)
			{
				u = k*m_nNumColumns + j;
				v = pnCol[k] * m_nNumColumns + j;
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}

		if (pnRow[k] != k)
		{
			for (i = 0; i <= m_nNumColumns - 1; i++)
			{
				u = i*m_nNumColumns + k;
				v = i*m_nNumColumns + pnRow[k];
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}
	}

	// 清理内存
	/*delete[] pnRow;
	delete[] pnCol;*/
	free(pnRow);
	free(pnCol);

	// 成功返回
	return true;
}

bool Inv_matrix(double *p, int n)
{
	void swap(double *a, double *b);
	int *is, *js, i, j, k, l;
	/*for (i = 0; i<n; i++)
	{
	putchar('\n');
	for (j = 0; j<n; j++)
	printf("%f  ", *(p + i*n + j));
	}
	puts("\n\n\n\n");*/
	double temp, fmax;
	is = (int *)malloc(n*sizeof(int));
	js = (int *)malloc(n*sizeof(int));
	for (k = 0; k < n; k++)
	{
		fmax = 0.0;
		for (i = k; i < n; i++)
		for (j = k; j<n; j++)
		{
			temp = fabs(*(p + i*n + j));//找最大值
			if (temp>fmax)
			{
				fmax = temp;
				is[k] = i; js[k] = j;
			}
		}
		if ((fmax + 1.0) == 1.0)
		{
			free(is); free(js);
			printf("no inv");
			return false;
		}
		if ((i = is[k]) != k)
		for (j = 0; j < n; j++)
			swap((p + k*n + j), (p + i*n + j));//交换指针
		if ((j = js[k]) != k)
		for (i = 0; i < n; i++)
			swap((p + i*n + k), (p + i*n + j));  //交换指针
		p[k*n + k] = 1.0 / p[k*n + k];
		for (j = 0; j < n; j++)
		if (j != k)
			p[k*n + j] *= p[k*n + k];
		for (i = 0; i < n; i++)
		if (i != k)
		for (j = 0; j < n; j++)
		if (j != k)
			p[i*n + j] = p[i*n + j] - p[i*n + k] * p[k*n + j];
		for (i = 0; i < n; i++)
		if (i != k)
			p[i*n + k] *= -p[k*n + k];
	}
	for (k = n - 1; k >= 0; k--)
	{
		if ((j = js[k]) != k)
		for (i = 0; i < n; i++)
			swap((p + j*n + i), (p + k*n + i));
		if ((i = is[k]) != k)
		for (j = 0; j < n; j++)
			swap((p + j*n + i), (p + j*n + k));
	}
	free(is);
	free(js);
	//打印出逆矩阵
	/*for (i = 0; i<n; i++)
	{
	putchar('\n');
	for (j = 0; j<n; j++)
	printf("%f  ", *(p + i*n + j));
	}
	puts("\n\n\n\n");*/
	return true;
}

/*计算单个图像的灰度（能量）方差*/
double Deviation(double *in_Data, int row, int col)
{
	double SumX = 0;
	double SumXX = 0;
	double dev;

	for (int i = 0; i < row*col; i++)
	{
		SumX += in_Data[i];
		SumXX += in_Data[i] * in_Data[i];
	}

	dev = sqrt(((SumXX - SumX * SumX / (row * col))) / (row * col - 1));

	return dev;
}
/*计算图像中某一位置矩形块内灰度方差*/
double Deviation(double *in_Data, int col, int row_begin, int row_count, int col_begin, int col_count)
{
	double SumX = 0;
	double SumXX = 0;
	double temp = 0;
	double dev;
	for (int r = 0; r < row_count; r++)
	{

		for (int i = 0; i < col_count; i++)
		{
			temp = in_Data[(r + row_begin)*col + col_begin + i];
			SumX += temp;
			SumXX += temp *temp;
		}
	}


	dev = sqrt(((SumXX - SumX * SumX / (row_count * col_count))) / (row_count * col_count - 1));

	return dev;
}


/*计算相关性系数*/
double correlation_coefficient(double *Bx, double *By, int N)
{
	//计算相关系数
	double Sum_xy = 0;
	double Sum_x = 0;
	double Sum_y = 0;
	double Sum_x2 = 0;
	double Sum_y2 = 0;
	for (int i = 0; i < N; i++)
	{
		Sum_x += Bx[i];
		Sum_y += By[i];
		Sum_xy += Bx[i] * By[i];
		Sum_x2 += Bx[i] * Bx[i];
		Sum_y2 += By[i] * By[i];
	}
	return (Sum_xy*N - Sum_x*Sum_y) / sqrt((N*Sum_x2 - pow(Sum_x, 2))*(N*Sum_y2 - pow(Sum_y, 2)));
}
/*计算绝对差*/
double SAD(double *Bx, double *By, int N)
{
	//计算绝对差
	double sum = 0;
	for (int i = 0; i < N; i++)
	{
		sum += fabs(Bx[i] - By[i]);
	}
	return sum;
}
//double SAD(double *Bx, double *By, int N,double threshold)
//{
//	//计算绝对差
//	double sum = 0;
//	for (int i = 0; i < N; i++)
//	{
//		sum += abs(Bx[i] - By[i]);
//		if (sum>threshold)
//			return 
//	}
//	return sum;
//}
/*最大类间差方法*/
double getOstu(double *in_Data, int row, int col)
{
	double max_gray, min_gray;
	max_gray = min_gray = in_Data[0];
	for (int i = 0; i < row*col; i++)
	{
		if (in_Data[i] < min_gray)
			min_gray = in_Data[i];
		if (in_Data[i]> max_gray)
			max_gray = in_Data[i];
	}

	double nHistogram[500] = { 0 };

	//
	for (int i = 0; i < row*col; i++)
	{
		nHistogram[int(in_Data[i] - min_gray)]++;
	}
	int threshold;//阈值
	long sum0 = 0, sum1 = 0; //存储前景的灰度总和及背景灰度总和    
	long cnt0 = 0, cnt1 = 0; //前景的总个数及背景的总个数    
	double w0 = 0, w1 = 0; //前景及背景所占整幅图像的比例    
	double u0 = 0, u1 = 0;  //前景及背景的平均灰度    
	double variance = 0; //最大类间方差 

	double maxVariance = 0;
	for (int i = 0; i < 500; i++)
	{
		sum0 = 0;
		sum1 = 0;
		cnt0 = 0;
		cnt1 = 0;
		w0 = 0;
		w1 = 0;
		for (int j = 0; j < i; j++)
		{
			cnt0 += nHistogram[j];
			sum0 += j * nHistogram[j];
		}

		u0 = (double)sum0 / cnt0;
		w0 = (double)cnt0 / (row*col);

		for (int j = i; j <= 255; j++)
		{
			cnt1 += nHistogram[j];
			sum1 += j * nHistogram[j];
		}
		u1 = (double)sum1 / cnt1;
		w1 = 1 - w0; // (double)cnt1 / size;   
		variance = w0 * w1 *  (u0 - u1) * (u0 - u1);
		if (variance > maxVariance)
		{
			maxVariance = variance;
			threshold = i;
		}
	}
	return threshold + min_gray;
}

double get1DMaxEntropyThreshold(double *in_Data, int Rows, int Cols)
{
	double max, min, max_min;
	max = min = in_Data[0];
	for (int i = 0; i < Rows*Cols; i++)
	{
		if (max<in_Data[i]) max = in_Data[i];
		if (min>in_Data[i]) min = in_Data[i];
	}
	max_min = max - min;

	//计算直方图
	int HistGram[256] = { 0 };
	for (int i = 0; i < Rows*Cols; i++)
	{
		int temp = (int)((in_Data[i] - min) * 255 / (max_min + 0.00001));
		if (temp <= 256 && temp >= 0)
			HistGram[temp]++;
	}
	int X, Y, Amount = 0;
	double HistGramD[256] = { 0 };
	double SumIntegral, EntropyBack, EntropyFore, MaxEntropy;
	int MinValue = 255, MaxValue = 0;
	double Threshold = 0;

	for (MinValue = 0; MinValue < 256 && HistGram[MinValue] == 0; MinValue++);
	for (MaxValue = 255; MaxValue > MinValue && HistGram[MinValue] == 0; MaxValue--);
	if (MaxValue == MinValue) return MaxValue;          // 图像中只有一个颜色             
	if (MinValue + 1 == MaxValue) return MinValue;      // 图像中只有二个颜色

	for (Y = MinValue; Y <= MaxValue; Y++) Amount += HistGram[Y];        //  像素总数

	for (Y = MinValue; Y <= MaxValue; Y++) HistGramD[Y] = (double)HistGram[Y] / Amount + 1e-17;

	MaxEntropy = DBL_MIN;
	for (Y = MinValue + 1; Y < MaxValue; Y++)
	{
		SumIntegral = 0;
		for (X = MinValue; X <= Y; X++) SumIntegral += HistGramD[X];
		EntropyBack = 0;
		for (X = MinValue; X <= Y; X++) EntropyBack += (-HistGramD[X] / SumIntegral * log(HistGramD[X] / SumIntegral));
		EntropyFore = 0;
		for (X = Y + 1; X <= MaxValue; X++) EntropyFore += (-HistGramD[X] / (1 - SumIntegral) * log(HistGramD[X] / (1 - SumIntegral)));
		if (MaxEntropy < EntropyBack + EntropyFore)
		{
			Threshold = Y;
			MaxEntropy = EntropyBack + EntropyFore;
		}
	}
	return (double)(Threshold / 255 * max_min + min);
}


double getGuDi(double *in_Data, int Rows, int Cols)
{
	double max, min, max_min;
	max = min = in_Data[0];
	for (int i = 0; i < Rows*Cols; i++)
	{
		if (max<in_Data[i]) max = in_Data[i];
		if (min>in_Data[i]) min = in_Data[i];
	}
	max_min = max - min;

	//计算直方图
	int hist[256] = { 0 };
	for (int i = 0; i < Rows*Cols; i++)
	{
		int temp = (int)((in_Data[i] - min) * 255 / (max_min + 0.00001) + 0.5);
		if (temp <= 256 && temp >= 0)
			hist[temp]++;
	}
	//
	double histC[256] = { 0 };
	double histCC[256] = { 0 };
	for (int i = 0; i < 256; i++)
	{
		histC[i] = hist[i];
		histCC[i] = hist[i];
	}
	int Iter = 0;
	while (IsDimodal(histC) == false)
	{
		histCC[0] = (histC[0] + histC[0] + histC[1]) / 3;
		for (int i = 1; i < 255; i++)
			histCC[i] = (histC[i - 1] + histC[i] + histC[i + 1]) / 3;
		histCC[255] = (histC[254] + histC[255] + histC[255]) / 3;
		for (int i = 0; i < 256; i++)
		{
			histC[i] = histCC[i];
		}
		Iter++;
		if (Iter >= 1000) return -1;
	}
	double histS[256] = { 0 };
	for (int i = 0; i < 256; i++)
	{
		histS[i] = (int)histCC[i];
	}
	//阈值为两峰之间的最小值
	bool Peakfound = false;
	for (int i = 0; i < 255; i++)
	{
		if (histCC[i - 1] < histCC[i] && histCC[i + 1] < histCC[i]) Peakfound = true;
		if (Peakfound == true && histCC[i - 1] >= histCC[i] && histCC[i + 1] >= histCC[i])
			return i - 1;
	}
	return -1;

}
bool IsDimodal(double * HistGram)       // 检测直方图是否为双峰的
{
	// 对直方图的峰进行计数，只有峰数位2才为双峰 
	int Count = 0;
	for (int Y = 1; Y < 255; Y++)
	{
		if (HistGram[Y - 1] < HistGram[Y] && HistGram[Y + 1] < HistGram[Y])
		{
			Count++;
			if (Count > 2) return false;
		}
	}
	if (Count == 2)
		return true;
	else
		return false;
}


void getRandom(int a[4], int max)
{
	srand(time(NULL)); //用时间做种，每次产生随机数不一样
	int i = 0;
	int temp = 0;
	a[0] = a[1] = a[2] = a[3] = max + 1;
	while (i < 4)
	{
		temp = rand() % max;
		if (temp != a[0] && temp != a[1] && temp != a[2] && temp != a[3])
		{
			a[i] = temp;
			i++;
		}
	}
}



/**
函数原型:
bool svd(vector<vector<double> > A, int K, std::vector<std::vector<double> > &U, std::vector<double> &S, std::vector<std::vector<double> > &V);
输入矩阵A,分解矩阵的秩K
输出U,S,V
本函数将A分解为U diag(S) V'
S[i],U[i],V[i]是A的第i大奇异值，及其对应的左歧义向量和右奇异向量
S,U,V的size由K指定
K是需要分解的rank，0<K<=min(m,n)

本程序采用的是最基本幂迭代算法，在linux g++下编译通过
**/
double get_norm(double *x, int n){
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * x[i];
	return sqrt(r);
}
double normalize(double *x, int n){
	double r = get_norm(x, n);
	if (r < eps)
		return 0;
	for (int i = 0; i < n; i++)
		x[i] /= r;
	return r;
}

inline double product(double*a, double *b, int n){
	double r = 0;
	for (int i = 0; i < n; i++)
		r += a[i] * b[i];
	return r;
}
void orth(double *a, double *b, int n){//|a|=1
	double r = product(a, b, n);
	for (int i = 0; i < n; i++)
		b[i] -= r*a[i];

}
void print(vector<vector<double> > &A){
	for (int i = 0; i < A.size(); i++){
		for (int j = 0; j < A[i].size(); j++){
			cout << setprecision(6) << A[i][j] << ' ';
		}
		cout << endl;
	}
}


bool svd(vector<vector<double> > A, int K,
	std::vector<std::vector<double> > &U,
	std::vector<double> &S, std::vector<std::vector<double> > &V)
{
	int M = A.size();
	int N = A[0].size();
	U.clear();
	V.clear();
	S.clear();
	S.resize(K, 0);
	U.resize(K);
	for (int i = 0; i < K; i++)
		U[i].resize(M, 0);
	V.resize(K);
	for (int i = 0; i < K; i++)
		V[i].resize(N, 0);


	srand(time(0));
	double *left_vector = new double[M];
	double *next_left_vector = new double[M];
	double *right_vector = new double[N];
	double *next_right_vector = new double[N];
	while (1){
		for (int i = 0; i<M; i++)
			left_vector[i] = (float)rand() / RAND_MAX;
		if (normalize(left_vector, M)>eps)
			break;
	}
	int col = 0;
	for (int col = 0; col < K; col++){
		double diff = 1;
		double r = -1;
		for (int iter = 0; diff >= eps && iter < MAX_ITER; iter++){
			memset(next_left_vector, 0, sizeof(double)*M);
			memset(next_right_vector, 0, sizeof(double)*N);
			for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
				next_right_vector[j] += left_vector[i] * A[i][j];

			r = normalize(next_right_vector, N);
			if (r < eps) break;
			for (int i = 0; i < col; i++)
				orth(&V[i][0], next_right_vector, N);
			normalize(next_right_vector, N);

			for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
				next_left_vector[i] += next_right_vector[j] * A[i][j];
			r = normalize(next_left_vector, M);
			if (r < eps) break;
			for (int i = 0; i < col; i++)
				orth(&U[i][0], next_left_vector, M);
			normalize(next_left_vector, M);
			diff = 0;
			for (int i = 0; i < M; i++){
				double d = next_left_vector[i] - left_vector[i];
				diff += d*d;
			}

			memcpy(left_vector, next_left_vector, sizeof(double)*M);
			memcpy(right_vector, next_right_vector, sizeof(double)*N);
		}
		if (r >= eps){
			S[col] = r;
			memcpy((char *)&U[col][0], left_vector, sizeof(double)*M);
			memcpy((char *)&V[col][0], right_vector, sizeof(double)*N);
		}
		else
			break;
	}
	delete[] next_left_vector;
	delete[] next_right_vector;
	delete[] left_vector;
	delete[] right_vector;

	return true;
}


//int main(){
//	int m = 10;
//	int n = 10;
//	srand(time(0));
//	vector<vector<double> > A;
//	A.resize(m);
//
//	for (int i = 0; i < m; i++){
//		A[i].resize(n);
//		for (int j = 0; j < n; j++)
//			A[i][j] = (float)rand() / RAND_MAX;
//	}
//	print(A);
//	cout << endl;
//
//	vector<vector<double> > U;
//	vector<double> S;
//	vector<vector<double> > V;
//	svd(A, 10, U, S, V);
//	cout << "U=" << endl;
//	print(U);
//	cout << endl;
//	cout << "S=" << endl;
//	for (int i = 0; i < S.size(); i++){
//		cout << S[i] << ' ';
//	}
//	cout << endl;
//	cout << "V=" << endl;
//	print(V);
//	return 0;
//}

bool Inv_svd(vector<vector<double> > A, vector<vector<double> > &invA, int n)
{
	vector<vector<double> > U;
	vector<double> S;
	vector<vector<double> > V;
	invA.resize(n);
	for (int i = 0; i < n; i++)
	{
		invA[i].resize(n, 0);
	}
	//print(A);
	//cout << endl;
	//SVD分解
	svd(A, n, U, S, V);//U[i][j]表示：第i列第j行
	//cout << "S=" << endl;
	//for (int i = 0; i < S.size(); i++){
	//	cout << S[i] << ' ';
	//}
	//cout << endl;
	/*将svd中转换成U[行][列]表示*/
	vector<vector<double> > u;
	vector<vector<double> > s_reciprocal;
	vector<vector<double> > v;
	u.resize(n); v.resize(n); s_reciprocal.resize(n);//初始化n行
	for (int i = 0; i < n; i++)
	{
		u[i].resize(n); v[i].resize(n);//每行初始化n列
		s_reciprocal[i].resize(n);
		if (S[i] != 0)
			s_reciprocal[i][i] = 1. / S[i];
		for (int j = 0; j < n; j++)
		{
			u[i][j] = U[j][i];
			v[i][j] = V[j][i];
		}
	}
	//cout << "u=" << endl;
	//print(u);
	//cout << endl;
	//cout << "v=" << endl;
	//print(v);
	//cout << "s_=" << endl;
	//print(s_reciprocal);

	vector<vector<double> > Temp;
	MatrixMultiplyV(v, s_reciprocal, Temp);
	MatrixMultiplyV(Temp, U, invA);
	//cout << "invA=" << endl;
	//print(invA);

	return true;

}
bool Inv_svdMN(vector<vector<double> > A, vector<vector<double> > &invA, int n)
{
	int M = A.size();//行数
	int N = A[0].size();//列数
	invA.resize(N);//mxn的矩阵A的伪逆为nxm
	for (int i = 0; i < N; i++)
	{
		invA[i].resize(M, 0);
	}
	vector<vector<double> > U;
	vector<double> S;
	vector<vector<double> > V;
	//print(A);
	//cout << endl;
	//SVD分解
	svd(A, n, U, S, V);//U[i][j]表示：第i列第j行
	//cout << "S=" << endl;
	//for (int i = 0; i < S.size(); i++){
	//	cout << S[i] << ' ';
	//}
	//cout << endl;
	/*将svd中转换成U[行][列]表示*/
	vector<vector<double> > u;
	vector<vector<double> > ut;
	vector<vector<double> > v;
	vector<vector<double> > dt_s;//S的伪逆转置
	u.resize(M); ut.resize(M); v.resize(N);
	dt_s.resize(N);//初始化n行
	//s_reciprocal[i].resize(n);
	for (int i = 0; i < N; i++)
	{
		dt_s[i].resize(M, 0);
		if (S[i] != 0)
			dt_s[i][i] = 1. / S[i];
	}
	for (int i = 0; i < M; i++)//真正的u
	{
		u[i].resize(M, 0);
		for (int j = 0; j < U.size(); j++)
		{
			u[i][j] = U[j][i];
		}
	}
	for (int i = 0; i < M; i++)//u的转置
	{
		ut[i].resize(M, 0);
		for (int j = 0; j < M; j++)
		{
			ut[i][j] = u[j][i];
		}
	}

	for (int i = 0; i < N; i++)
	{
		v[i].resize(N, 0);
		for (int j = 0; j < N; j++)
		{
			v[i][j] = V[j][i];
		}
	}
	//cout << "u=" << endl;
	//print(u);
	//cout << endl;
	//cout << "v=" << endl;
	//print(v);
	//cout << endl;
	//cout << "dt_s=" << endl;
	//print(dt_s);
	//cout << endl;

	vector<vector<double> > Temp;
	MatrixMultiplyV(v, dt_s, Temp);
	MatrixMultiplyV(Temp, ut, invA);
	//cout << "invA=" << endl;
	//print(invA);

	return true;

}
bool MatrixMultiplyV(vector<vector<double> > X1, vector<vector<double> > X2, vector<vector<double> > &Y)
{

	int row1 = X1.size();
	int col1 = X1[0].size();
	int row2 = X2.size();
	int col2 = X2[0].size();
	if (col1 != row2)
		return false;
	Y.resize(row1);
	for (int i = 0; i < row1; i++)
		Y[i].resize(col2);

	for (int i = 0; i < Y.size(); i++)
	{
		for (int j = 0; j < Y[0].size(); j++)
		{
			double temp = 0;
			for (int ki = 0; ki < col1; ki++)
			{
				temp += X1[i][ki] * X2[ki][j];
			}
			Y[i][j] = temp;
		}
	}
	return true;

}
//均值计算
double MeanImg(double *in_Data, int Rows, int Cols)
{
	double sum = 0, mean = 0;
	for (int i = 0; i < Rows*Cols; i++)
	{
		sum += in_Data[i];
	}
	mean = sum / (Rows*Cols);
	return mean;
}