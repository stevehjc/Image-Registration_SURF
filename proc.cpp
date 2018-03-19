#define _CRT_SECURE_NO_WARNINGS


#include"proc.h"
#include"math_m.h"
#include"mSURF.h"

void * fspace_1d(int col, int length)
{
	void *b;

	b = (void *)calloc(col, length);
	if (!b) return NULL;
	return(b);
}

int WriteBMPFile(unsigned char *image, int Row, int Col, char *FileName)
{
	FILE *fp;
	_BITMAPFILEHEADER FileHeader;
	_BITMAPINFOHEADER InfoHeader;
	_RGBQUAD bmiColors[256];

	short i;
	short linebytes, skipstep, skip[3];
	long off;

	for (i = 0; i<3; i++)
		skip[i] = 0;

	FileHeader.bfType = (short)0x4d42;
	FileHeader.bfReserved1 = 0;
	FileHeader.bfReserved2 = 0;

	InfoHeader.biSize = 40;
	InfoHeader.biPlanes = 1;
	InfoHeader.biCompression = 0;
	InfoHeader.biXPelsPerMeter = 0;
	InfoHeader.biYPelsPerMeter = 0;
	InfoHeader.biClrUsed = 0;
	InfoHeader.biClrImportant = 0;

	for (i = 0; i<256; i++)
	{
		bmiColors[i].rgbRed = (unsigned char)i;
		bmiColors[i].rgbGreen = (unsigned char)i;
		bmiColors[i].rgbBlue = (unsigned char)i;
		bmiColors[i].rgbReserved = 0;
	}

	/*attempt to open the file*/
	if ((fp = fopen(FileName, "wb")) == NULL) return 0;

	linebytes = ((Col + 3) / 4) * 4;
	skipstep = linebytes - Col;

	FileHeader.bfSize = Row*Col + 54 + 256 * 4;
	FileHeader.bfOffBits = 54 + 256 * 4;

	InfoHeader.biWidth = Col;
	InfoHeader.biHeight = Row;
	InfoHeader.biBitCount = 8;
	InfoHeader.biSizeImage = Row*Col;


	//fwrite(&FileHeader,sizeof(BITMAPFILEHEADER),1,fp);
	fwrite(&FileHeader.bfType, 1, sizeof(unsigned short), fp);
	fwrite(&FileHeader.bfSize, 1, sizeof(unsigned long), fp);
	fwrite(&FileHeader.bfReserved1, 1, sizeof(unsigned short), fp);
	fwrite(&FileHeader.bfReserved2, 1, sizeof(unsigned short), fp);
	fwrite(&FileHeader.bfOffBits, 1, sizeof(unsigned long), fp);

	//	fwrite(&FileHeader,14,1,fp);
	fwrite(&InfoHeader, sizeof(_BITMAPINFOHEADER), 1, fp);
	fwrite(bmiColors, sizeof(_RGBQUAD), 256, fp);

	for (i = Row - 1; i >= 0; i--)
	{
		off = (long)(Col)*(long)(i + 1);
		fwrite(&image[i*Col], 1, Col, fp);
		fwrite(skip, 1, skipstep, fp);
	}
	fclose(fp);

	return 1;

}
unsigned char * ReadBMPFile(int *Row, int *Col, char *FileName)
{
	FILE *fp;
	_BITMAPFILEHEADER BmpFileHeader;
	_BITMAPINFOHEADER BmpInfoHeader;
	int iRow, iCol;
	unsigned char *Image;


	short i, j;
	short bytes, linebytes, skip;
	unsigned char d1, d2, d3;
	long off;

	/*attempt to open the file*/
	if ((fp = fopen(FileName, "rb")) == NULL) return NULL;
	//    if (fread(&BmpFileHeader,1,sizeof(BITMAPFILEHEADER),fp)!=sizeof(BITMAPFILEHEADER))
	if (fread(&BmpFileHeader.bfType, 1, sizeof(unsigned short), fp) != sizeof(unsigned short)) return NULL;
	if (fread(&BmpFileHeader.bfSize, 1, sizeof(unsigned long), fp) != sizeof(unsigned long)) return NULL;
	if (fread(&BmpFileHeader.bfReserved1, 1, sizeof(unsigned short), fp) != sizeof(unsigned short)) return NULL;
	if (fread(&BmpFileHeader.bfReserved2, 1, sizeof(unsigned short), fp) != sizeof(unsigned short)) return NULL;
	if (fread(&BmpFileHeader.bfOffBits, 1, sizeof(unsigned long), fp) != sizeof(unsigned long)) return NULL;

	if (fread(&BmpInfoHeader, 1, sizeof(_BITMAPINFOHEADER), fp) != sizeof(_BITMAPINFOHEADER))
		return NULL;

	iCol = BmpInfoHeader.biWidth;
	iRow = BmpInfoHeader.biHeight;

	Image = (unsigned char *)fspace_1d(iRow*iCol, sizeof(unsigned char));

	if (BmpInfoHeader.biBitCount == 8) {
		bytes = iCol;
		linebytes = ((bytes + 3) / 4) * 4;
		skip = linebytes - bytes;
		for (i = 0; i<iRow; i++)
		{
			off = (long)(linebytes)*(long)(i + 1);
			fseek(fp, -off, SEEK_END);
			fread(&Image[i*iCol], iCol, 1, fp);
		}
	}
	else if (BmpInfoHeader.biBitCount == 24) {
		bytes = iCol * 3;
		linebytes = ((bytes + 3) / 4) * 4;
		for (i = 0; i<iRow; i++)
		{
			j = 0;
			off = (long)(linebytes)*(long)(i + 1);
			fseek(fp, -off, SEEK_END);
			while (j<iCol)
			{
				fread(&d1, 1, 1, fp);
				fread(&d2, 1, 1, fp);
				fread(&d3, 1, 1, fp);
				Image[i*iCol + j] = (unsigned char)(d1*0.299 + d2*0.587 + d3*0.114);
				j++;
			}
		}
	}

	fclose(fp);
	*Row = iRow;
	*Col = iCol;
	return Image;

}

bool mSURF_Detection(double *srcImg1, double *srcImg2, int Rows, int Cols,
	vector<Point2D> &matchpoints1, vector<Point2D> &matchpoints2, int matches_num)
{

	//（一）转换格式为IMAGE,计算积分图像
	IMAGE* srcIMAGE1 = toIMAGE(srcImg1, Rows, Cols);//IMAGE 格式
	IMAGE* srcIMAGE2 = toIMAGE(srcImg2, Rows, Cols);//IMAGE 格式
	IMAGE* intergalSum1 = CreatImage(Rows + 1, Cols + 1);//积分和图像
	IMAGE* intergalSum2 = CreatImage(Rows + 1, Cols + 1);//积分和图像
	intergal(srcIMAGE1, intergalSum1);//计算积分图像
	intergal(srcIMAGE2, intergalSum2);//计算积分图像
	//Mat sumMat1 = toMat((double*)intergalSum1->gray, intergalSum1->row, intergalSum1->col);
	//Mat sumMat2 = toMat((double*)intergalSum2->gray, intergalSum2->row, intergalSum2->col);
	//（二）检测SURF特征点
	vector<mKeyPoint> keypoints1, keypoints2;
	float hessianThreshold = 800;//hessian阈值1300
	//SURF检测特征点
	SURF_detect(srcIMAGE1, intergalSum1, keypoints1, hessianThreshold);
	SURF_detect(srcIMAGE2, intergalSum2, keypoints2, hessianThreshold);
	printf("----------\n检测到的特征点：%d \t %d \n", keypoints1.size(), keypoints2.size());
	//SURF计算角度和特征点描述矩阵
	if (!keypoints1.size() || !keypoints2.size())//如果检测到的特征点个数为0，返回
	{
		return false;
	}
	/**********通过限制特征点个数，缩短匹配时间，提高效率************/

	int Nk = 2;//保留匹配点数系数

	if (1)
	{
		if (keypoints1.size() > Nk * matches_num)
		{
			keypoints1.erase(keypoints1.begin() + Nk*matches_num, keypoints1.end());
			printf("keypoints1 erase：%d ", Nk*matches_num);
		}
		if (keypoints2.size() > Nk * matches_num)
		{
			keypoints2.erase(keypoints2.begin() + Nk * matches_num, keypoints2.end());
			printf("keypoints2 erase：%d\n", Nk * matches_num);
		}
	}
	if (0)
	{
		if (keypoints1.size() > Nk * matches_num)
		{
			double temp = keypoints1[Nk*matches_num].response;
			int k2index = Nk*matches_num;
			for (int i = 0; i < keypoints2.size(); i++)
			{
				if (keypoints2[i].response < temp)
				{
					k2index = i;
					break;
				}
			}
			keypoints1.erase(keypoints1.begin() + Nk*matches_num, keypoints1.end());
			printf("keypoints1 erase：%d ", Nk*matches_num);
			keypoints2.erase(keypoints2.begin() + k2index, keypoints2.end());
			printf("keypoints2 erase：%d\n", Nk * matches_num);
		}
	}



	//（三）特征点描述

	IMAGE* descriptor1 = CreatImage(keypoints1.size(), 64);
	IMAGE* descriptor2 = CreatImage(keypoints2.size(), 64);
	SURF_descriptor(srcIMAGE1, intergalSum1, keypoints1, descriptor1);
	SURF_descriptor(srcIMAGE2, intergalSum2, keypoints2, descriptor2);

	//（四）匹配特征点
	vector<mDMatch> matchePoints;
	for (int i = 0; i < keypoints1.size(); i++)
	{
		//descriptor1的i行首
		double *p1 = (double*)(descriptor1->gray + i * 64 * sizeof(double));
		double distance_min = 100;
		//double var_max = -1;//
		double var_t = 0;
		int j_min = 0, i_min = 0;
		int j;
		double dx = 0, dy = 0;
		for (j = 0; j < keypoints2.size(); j++)
		{
			dx = fabs(keypoints1[i].pt.x - keypoints2[j].pt.x);
			dy = fabs(keypoints1[i].pt.y - keypoints2[j].pt.y);

			double *p2 = (double*)(descriptor2->gray + j * 64 * sizeof(double));
			//Mat tempH2(1, 64, CV_64FC1, p2);
			//计算距离
			double distance_t = calcDistance(p1, p2, 64);
			//double distance_t = SAD(p1, p2, 64);
			//double var_t = correlation_coefficient(p1, p2, 64);
			//if (var_t > var_max)
			//{
			//	var_max = var_t;
			//	//var_max = var_t;
			//	i_min = i;
			//	j_min = j;
			//}
			if (distance_t < distance_min)
			{
				distance_min = distance_t;
				i_min = i;
				j_min = j;
			}
		}
		matchePoints.push_back(mDMatch(i_min, j_min, distance_min));
	}

	sort(matchePoints.begin(), matchePoints.end(), mathcePointsBetter()); //特征点排序

	//转换为Point2D格式
	//vector<Point2D> p1, p2;
	int goodPoints;
	//= (int)(matchePoints.size() * 2 / 3);
	goodPoints = matches_num<matchePoints.size() ? matches_num:matchePoints.size();
	for (int i = 0; i < goodPoints; i++)
	{
		matchpoints1.push_back(Point2D(keypoints1[matchePoints[i].queryIdx].pt.x, keypoints1[matchePoints[i].queryIdx].pt.y));
		matchpoints2.push_back(Point2D(keypoints2[matchePoints[i].trainIdx].pt.x, keypoints2[matchePoints[i].trainIdx].pt.y));
	}

	//for (int i = 0; i < matchpoints1.size(); i++)
	//{
	//	printf("%f %f \t %f %f\n", matchpoints1[i].x, matchpoints1[i].y,
	//		matchpoints2[i].x, matchpoints2[i].y);
	//}


	//3.绘制特征点，分析分布特征
	/*******************************/

	//重要：释放内存
	DeleteImage(intergalSum1);
	DeleteImage(intergalSum2);
	DeleteImage(descriptor1);
	DeleteImage(descriptor2);


	/********配准性能评价-附加内容**********/

	if (matchpoints1.size() && matchpoints2.size())
		return true;
	else
		return false;

}
void mWarpPerspective(double* srcImg1, double* imageTransform1, vector<double>para, int Rows, int Cols)
{

	double homo[3][3];
	if (para.size() == 8)
	{
		homo[0][0] = para[0]; homo[0][1] = para[1]; homo[0][2] = para[2];
		homo[0][3] = para[3]; homo[0][4] = para[4]; homo[0][5] = para[5];
		homo[0][6] = para[6]; homo[0][7] = para[7]; homo[0][8] = 1;
	}
	Inv_matrix(homo[0], 3);

	double u, v, w;
	double xx, yy;//前一帧图像中坐标，（带小数的）


	for (int y = 0; y < Rows; y++)
	{
		for (int x = 0; x < Cols; x++)
		{
			u = homo[0][0] * x + homo[0][1] * y + homo[0][2];
			v = homo[0][3] * x + homo[0][4] * y + homo[0][5];
			w = homo[0][6] * x + homo[0][7] * y + homo[0][8];
			xx = u / w;
			yy = v / w;

			int _i, _j;
			double _idec, _jdec;
			if (xx < 0)
			{
				_i = (int)xx - 1;
			}
			else
			{
				_i = (int)xx;
			}
			if (yy < 0)
			{
				_j = (int)yy - 1;
			}
			else
			{
				_j = (int)yy;
			}
			_idec = xx - _i;
			_jdec = yy - _j;

			double P_ji, P_ji1, P_j1i, P_j1i1;//双线性插值法

			//P_ji;
			if (_j >= 0 && _j < Rows && _i >= 0 && _i < Cols)
			{
				P_ji = srcImg1[_j*Cols + _i];
			}
			else
			{
				P_ji = 0;
			}
			//P_ji1;
			if (_j >= 0 && _j < Rows && _i >= 0 && _i < (Cols - 1))
			{
				P_ji1 = srcImg1[_j*Cols + _i + 1];
			}
			else
			{
				P_ji1 = 0;
			}
			//P_j1i;
			if (_j >= 0 && _j < (Rows - 1) && _i >= 0 && _i < Cols)
			{
				P_j1i = srcImg1[(_j + 1)*Cols + _i];
			}
			else
			{
				P_j1i = 0;
			}
			//P_j1i1;
			if (_j >= 0 && _j < (Rows - 1) && _i >= 0 && _i < (Cols - 1))
			{
				P_j1i1 = srcImg1[(_j + 1)*Cols + _i + 1];
			}
			else
			{
				P_j1i1 = 0;
			}

			imageTransform1[y*Cols + x] = (1 - _idec)*(1 - _jdec)*P_ji
				+ _idec*(1 - _jdec)*P_ji1
				+ _jdec*(1 - _idec)*P_j1i
				+ _idec*_jdec*P_j1i1;

			//imageTransform1[y*Cols + x] = (1 - _idec)*(1 - _jdec)*srcImg1[_j*Cols + _i]
			//	+ _idec*(1 - _jdec)*srcImg1[_j*Cols + _i + 1]
			//	+ _jdec*(1 - _idec)*srcImg1[(_j + 1)*Cols + _i]
			//	+ _idec*_jdec*srcImg1[(_j + 1)*Cols + _i + 1];

		}
	}

}

int WriteBMPFile(double *image, int Row, int Col, char *FileName)
{
	unsigned char *img;
	double max, min, max_min;
	int i;
	int res;

	img = (unsigned char *)fspace_1d(Row*Col, sizeof(unsigned char));
	max = min = image[0];
	for (i = 0; i<Row*Col; i++)
	{
		if (max<image[i]) max = image[i];
		if (min>image[i]) min = image[i];
	}
	max_min = max - min;

	for (i = 0; i<Row*Col; i++)
	{
		img[i] = (unsigned char)((image[i] - min) * 255 / (max_min + 0.0001));
	}

	res = WriteBMPFile(img, Row, Col, FileName);

	free((void*)img);

	return res;
}
vector<double> getRANSAC2(vector<Point2D>p1, vector<Point2D>p2, double threshold, int count)
{
	vector<double> homoV;
	if (p1.size()==0)
		return homoV;
	//计算最大内点集合
	int max_iters = 2000;
	int iters = max_iters;
	int innerP, max_innerP = 0;
	vector<int> innerPvInd;//内点集合索引-临时
	vector<int> innerPvInd_i;
	vector<Point2D> selectP1, selectP2;
	double *homo = NULL;//透视变换参数
	double *trans;
	int k = 0;
	//生成随机表
	int selectIndex[2000][4];
	srand(time(NULL)); //用时间做种，每次产生随机数不一样
	int pCount = p1.size();

	for (int i = 0; i < 2000; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			int ii = 0;
			int temp = 0;
			selectIndex[i][0] = selectIndex[i][1] = selectIndex[i][2] = selectIndex[i][3] = pCount + 1;
			while (ii < 4)
			{
				temp = rand() % pCount;
				if (temp != selectIndex[i][0] &&
					temp != selectIndex[i][1] &&
					temp != selectIndex[i][2] &&
					temp != selectIndex[i][3])
				{
					selectIndex[i][ii] = temp;
					ii++;
				}
			}
		}
	}
	for (; k < iters; k++)
	{
		//int selectIndex[4] = { 0, 1, 2, 3 };
		//getRandom(selectIndex, p1.size());//获得不重复的四个坐标点索引
		//测试随机性
		//printf("Random:%d %d %d %d\n", selectIndex[0], selectIndex[1], selectIndex[2], selectIndex[3]);
		selectP1.push_back(p1[selectIndex[k][0]]);
		selectP1.push_back(p1[selectIndex[k][1]]);
		selectP1.push_back(p1[selectIndex[k][2]]);
		selectP1.push_back(p1[selectIndex[k][3]]);
		selectP2.push_back(p2[selectIndex[k][0]]);
		selectP2.push_back(p2[selectIndex[k][1]]);
		selectP2.push_back(p2[selectIndex[k][2]]);
		selectP2.push_back(p2[selectIndex[k][3]]);

		//trans = getPerspectiveTransformI(selectP1, selectP2);//计算透视变换矩阵（目前对奇异矩阵逆无法求解，后续将改进奇异分解求逆方法）
		trans = getPerspectiveTransformIata(selectP1, selectP2);
		if (trans == NULL)
		{
			printf("透视变换矩阵变换错误");
			selectP1.clear();
			selectP2.clear();
			continue;
		}
		//计算模型参数误差，如果误差大于阈值，放弃这组模型参数
		innerP = 0;
		double u, v, w;
		double errX, errY;
		for (int i = 0; i < p1.size(); i++)
		{
			errX = errY = 0;
			u = p1[i].x*trans[0] + p1[i].y*trans[1] + trans[2];
			v = p1[i].x*trans[3] + p1[i].y*trans[4] + trans[5];
			w = p1[i].x*trans[6] + p1[i].y*trans[7] + 1;
			errX = fabs(u / w - p2[i].x);
			errY = fabs(v / w - p2[i].y);
			if (threshold>(errX*errX + errY*errY))
			{
				innerP++;
				innerPvInd.push_back(i);
			}
		}
		if (innerP>max_innerP)
		{
			max_innerP = innerP;
			homo = trans;
			innerPvInd_i = innerPvInd;
			//更新迭代次数
			double p = 0.995;
			double ep = (double)(p1.size() - innerP) / p1.size();
			// avoid inf's & nan's
			double num_ = std::max(1. - p, DBL_MIN);
			double denom_ = 1. - pow(1. - ep, 4);

			if (denom_ < DBL_MIN)
				iters = 0;
			else
			{
				double num = log(num_);
				double denom = log(denom_);
				iters = (denom >= 0 || -num >= max_iters*(-denom) ? max_iters : (int)(num / denom));
			}

		}
		//if (innerP / p1.size()>0.95)//如果内点比例大于给定值，则跳出
		//{
		//	break;
		//}
		selectP1.clear();
		selectP2.clear();
		innerPvInd.clear();
	}
	printf("RANSAC内点比例-循环次数：%d %d %d \t\n", max_innerP, p1.size(), k);

	//根据内点计算投影矩阵参数

	vector<Point2D > _p1, _p2;
	for (int i = 0; i < max_innerP; i++)
	{
		_p1.push_back(Point2D(p1[innerPvInd_i[i]].x, p1[innerPvInd_i[i]].y));
		_p2.push_back(Point2D(p2[innerPvInd_i[i]].x, p2[innerPvInd_i[i]].y));
	}
	
	homoV = getPerspectiveTransformLSM2(_p1, _p2);
	return homoV;
}
double* getPerspectiveTransformIata(vector<Point2D> src, vector<Point2D> dst)
{

	const int NN = 8;//
	int count = src.size();
	vector<vector<double> > a;
	vector<vector<double> > b;
	a.resize(2 * count); b.resize(2 * count);
	for (int i = 0; i < 2 * count; i++)
	{
		a[i].resize(NN, 0); b[i].resize(1, 0);
	}
	for (int i = 0; i < count; i++)
	{
		a[i][0] = a[i + count][3] = src[i].x;
		a[i][1] = a[i + count][4] = src[i].y;
		a[i][2] = a[i + count][5] = 1;
		a[i][3] = a[i][4] = a[i][5] = a[i + count][0] = a[i + count][1] = a[i + count][2] = 0;
		a[i][6] = -src[i].x*dst[i].x;
		a[i][7] = -src[i].y*dst[i].x;
		a[i + count][6] = -src[i].x*dst[i].y;
		a[i + count][7] = -src[i].y*dst[i].y;
		b[i][0] = dst[i].x;
		b[i + count][0] = dst[i].y;
	}

	/*cout << "a=" << endl;
	print(a);
	cout << endl;
	cout << "b=" << endl;
	print(b);
	cout << endl;*/
	/*正规矩阵求最小二乘*/
	double aa[8][8], bb[8];
	vector<vector<double> > at(NN, vector<double>(a.size()));
	for (int i = 0; i < NN; i++)//求a的转置
	{
		for (int j = 0; j < 2 * count; j++)
		{
			at[i][j] = a[j][i];
		}
	}
	vector<vector<double>> ata(NN, vector<double>(NN, 0));

	vector<vector<double>> atb(NN, vector<double>(1, 0));
	if (MatrixMultiplyV(at, a, ata) && MatrixMultiplyV(at, b, atb))
	{
		for (int i = 0; i < NN; i++)
		{
			for (int j = 0; j < NN; j++)
			{
				aa[i][j] = ata[i][j];
			}
		}
		for (int i = 0; i < NN; i++)
		{
			bb[i] = atb[i][0];
		}
	}

	vector<double> transV;
	if (Inv_matrix(aa[0], NN))
	{
		double *trans = NULL;
		trans = MatrixMultiply(aa[0], NN, NN, bb, 1);
		return trans;
	}
	return NULL;
}
vector<double> getPerspectiveTransformLSM(vector<Point2D> src, vector<Point2D> dst)
{
	int count = src.size();
	vector<vector<double> > a;
	vector<vector<double> > b;
	a.resize(2 * count); b.resize(2 * count);
	for (int i = 0; i < 2 * count; i++)
	{
		a[i].resize(8, 0); b[i].resize(1, 0);
	}
	for (int i = 0; i < count; i++)
	{
		a[i][0] = a[i + count][3] = src[i].x;
		a[i][1] = a[i + count][4] = src[i].y;
		a[i][2] = a[i + count][5] = 1;
		a[i][3] = a[i][4] = a[i][5] = a[i + count][0] = a[i + count][1] = a[i + count][2] = 0;
		a[i][6] = -src[i].x*dst[i].x;
		a[i][7] = -src[i].y*dst[i].x;
		a[i + count][6] = -src[i].x*dst[i].y;
		a[i + count][7] = -src[i].y*dst[i].y;
		b[i][0] = dst[i].x;
		b[i + count][0] = dst[i].y;
	}
	//cout << "a=" << endl;
	//print(a);
	//cout << endl;
	//cout << "b=" << endl;
	//print(b);
	//cout << endl;
	vector<vector<double> > invA;
	vector<vector<double> > _transV;
	vector<double> transV;
	if (Inv_svdMN(a, invA, 8))
	{
		if (MatrixMultiplyV(invA, b, _transV))
		{
			for (int i = 0; i < 8; i++)
				transV.push_back(_transV[i][0]);
		}
	}

	return transV;

}
vector<double> getPerspectiveTransformLSM2(vector<Point2D> src, vector<Point2D> dst)
{
	const int NN = 8;//
	int count = src.size();
	vector<vector<double> > a;
	vector<vector<double> > b;
	a.resize(2 * count); b.resize(2 * count);
	for (int i = 0; i < 2 * count; i++)
	{
		a[i].resize(NN, 0); b[i].resize(1, 0);
	}
	for (int i = 0; i < count; i++)
	{
		a[i][0] = a[i + count][3] = src[i].x;
		a[i][1] = a[i + count][4] = src[i].y;
		a[i][2] = a[i + count][5] = 1;
		a[i][3] = a[i][4] = a[i][5] = a[i + count][0] = a[i + count][1] = a[i + count][2] = 0;
		a[i][6] = -src[i].x*dst[i].x;
		a[i][7] = -src[i].y*dst[i].x;
		a[i + count][6] = -src[i].x*dst[i].y;
		a[i + count][7] = -src[i].y*dst[i].y;
		b[i][0] = dst[i].x;
		b[i + count][0] = dst[i].y;
	}

	/*cout << "a=" << endl;
	print(a);
	cout << endl;
	cout << "b=" << endl;
	print(b);
	cout << endl;*/
	/*正规矩阵求最小二乘*/
	double aa[8][8], bb[8];
	vector<vector<double> > at(NN, vector<double>(a.size()));
	for (int i = 0; i < NN; i++)//求a的转置
	{
		for (int j = 0; j < 2 * count; j++)
		{
			at[i][j] = a[j][i];
		}
	}
	vector<vector<double>> ata(NN, vector<double>(NN, 0));

	vector<vector<double>> atb(NN, vector<double>(1, 0));
	if (MatrixMultiplyV(at, a, ata) && MatrixMultiplyV(at, b, atb))
	{
		for (int i = 0; i < NN; i++)
		{
			for (int j = 0; j < NN; j++)
			{
				aa[i][j] = ata[i][j];
			}
		}
		for (int i = 0; i < NN; i++)
		{
			bb[i] = atb[i][0];
		}
	}

	vector<double> transV;
	if (Inv_matrix(aa[0], NN))
	{
		double *trans = NULL;
		trans = MatrixMultiply(aa[0], NN, NN, bb, 1);
		//return trans;
		for (int i = 0; i < NN; i++)
		{
			transV.push_back(trans[i]);
		}
	}
	return transV;




	//cout << "a=" << endl;
	//print(a);
	//cout << endl;
	//cout << "b=" << endl;
	//print(b);
	//cout << endl;
	//vector<vector<double> > invA;
	//vector<vector<double> > _transV;
	//vector<double> transV;
	//if (Inv_svdMN(a, invA, 8))
	//{
	//	if (MatrixMultiplyV(invA, b, _transV))
	//	{
	//		for (int i = 0; i < 8; i++)
	//			transV.push_back(_transV[i][0]);
	//	}
	//}

	//return transV;

}
void resizeVV(vector<vector<double>> src, vector<vector<double>>& dst, int interpolation)
{
	//注意，src和dst表示正方形矩阵

	double dsize = (double)src.size() / dst.size();
	if (dst.size() && dst[0].size())
	{
		for (int i = 0; i < dst.size(); i++)
		{
			for (int j = 0; j < dst[0].size(); j++)
			{
				int _i = (int)i*dsize;//整数部分
				double _idec = i*dsize - _i;//小数部分
				int _j = (int)j*dsize;
				double _jdec = j*dsize - _j;
				if (_j >= 0 && _j < (src.size() - 1) && _i >= 0 && _i < (src.size() - 1))//双线性插值法
				{
					dst[i][j] = (1 - _idec)*(1 - _jdec)*src[_i][_j]
						+ _idec*(1 - _jdec)*src[_i + 1][_j]
						+ _jdec*(1 - _idec)*src[_j][_j + 1]
						+ _idec*_jdec*src[_i + 1][_j + 1];
				}
			}
		}

	}

}
double calcDistance(double* p1, double* p2, int n)
{
	double distance = 0;
	for (int i = 0; i < n; i++)
	{
		distance += fabs(p1[i] - p2[i]);
	}
	return distance;
}