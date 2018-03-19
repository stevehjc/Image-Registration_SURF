#define _CRT_SECURE_NO_WARNINGS
#include"mSURF.h"


IMAGE* toIMAGE(double* src, int row, int col)
{
	IMAGE* srcImg;
	srcImg = (IMAGE *)(malloc(sizeof(IMAGE)));
	srcImg->col = col;
	srcImg->row = row;
	srcImg->imagesize = row*col;
	srcImg->gray = (unsigned char*)src;

	return srcImg;
}


IMAGE* CreatImage(int row, int col)
{
	IMAGE *image;

	if (!(image = (IMAGE *)(malloc(sizeof(IMAGE)))))
	{
		return(NULL);
	}

	image->row = row;
	image->col = col;
	image->imagesize = (long)row*col;
	image->firstcol = image->firstrow = 0;
	image->lastrow = row;
	image->lastcol = col;
	image->scale = 1.0;
	strcpy(image->cmt, "?");
	strcpy(image->name, "?");

	image->gray = (unsigned char*)calloc(image->imagesize, sizeof(double));
	if (!image->gray)
	{
		free(image);
		return NULL;
	}

	image->previous = NULL;
	image->next = NULL;

	return (image);

}
void DeleteImage(IMAGE * image)
{
	if (image == NULL)
	{
		return;
	}
	if (image->gray != NULL) free(image->gray);
	image->gray = NULL;
	free(image);
	image = NULL;

}
void SURF_detect(IMAGE* srcImg, IMAGE* intergalSum, vector<mKeyPoint>&keypoints, float hessianThreshold)
{


	int nOctaves = 4;
	int nOctaveLayers = 2;
	mFastHessianDetector(intergalSum, keypoints, nOctaves, nOctaveLayers, hessianThreshold);

}
void SURF_descriptor(IMAGE* srcImg, IMAGE* intergalSum, vector<mKeyPoint>& keypoints, IMAGE* descriptors)
{
	int N = keypoints.size();
	if (N > 0)
	{

		mSURFInvoker(srcImg, intergalSum, keypoints, descriptors);


	}
}

void mFastHessianDetector(const IMAGE* sum, vector<mKeyPoint> &keypoints,
	int nOctaves, int nOctaveLayers, float hessianThreshold)
{
	int SAMPLE_SETEPO = 1;
	int nTotalLayers = (nOctaveLayers + 2)*nOctaves;
	int nMiddleLayers = nOctaveLayers*nOctaves;

	int SURF_HAAR_SIZE0 = 9;
	int SURF_HAAR_SIXE_INC = 6;

	vector<IMAGE*> dets(nTotalLayers);
	vector<IMAGE*> traces(nTotalLayers);
	vector<int> sizes(nTotalLayers);
	vector<int> sampleSteps(nTotalLayers);
	vector<int> middleIndices(nMiddleLayers);
	keypoints.clear();

	int index = 0, middleIndex = 0, step = SAMPLE_SETEPO;
	for (int octave = 0; octave < nOctaves; octave++)
	{
		for (int layer = 0; layer < nOctaveLayers + 2; layer++)
		{
			dets[index] = CreatImage((sum->row - 1) / step, (sum->col - 1) / step);
			traces[index] = CreatImage((sum->row - 1) / step, (sum->col - 1) / step);
			sizes[index] = (SURF_HAAR_SIZE0 + SURF_HAAR_SIXE_INC*layer) << octave;
			sampleSteps[index] = step;

			if (0 < layer&&layer <= nOctaveLayers)
				middleIndices[middleIndex++] = index;
			index++;
		}
		step *= 2;
	}


	mSURFBuildInvoker(sum, sizes, sampleSteps, dets, traces);

	mSURFFindInvoker(sum, dets, traces, sizes,
		sampleSteps, middleIndices, keypoints, nOctaveLayers, hessianThreshold);
	sort(keypoints.begin(), keypoints.end(), mKeypointGreater());


	for (int i = 0; i < nTotalLayers; i++)
	{
		DeleteImage(dets[i]);
		DeleteImage(traces[i]);
	}

}
void intergal(IMAGE *img, IMAGE* sum)
{
	int _row = img->row;
	int _col = img->col;
	int row = sum->row;
	int col = sum->col;

	double *inData = (double*)img->gray;
	double *outData = (double*)sum->gray;

	for (int j = 1; j < row; j++)
	{
		for (int i = 1; i < col; i++)
		{
			outData[j*col + i] = 0;
		}
	}
	for (int j = 1; j < row; j++)
	{
		for (int i = 1; i < col; i++)
		{
			outData[j*col + i] = inData[(j - 1)*_col + i - 1]
				+ outData[j*col + i - 1]
				+ outData[(j - 1)*col + i]
				- outData[(j - 1)*col + i - 1];
		}
	}

}
void mSURFBuildInvoker(const IMAGE* sum, vector<int>& sizes, vector<int>& sampleSteps, vector<IMAGE*>& dets, vector<IMAGE*> traces)
{
	int N = sizes.size();
	for (int i = 0; i < N; i++)
	{
		mcalcLayerDetAndTrace(sum, sizes[i], sampleSteps[i], dets[i], traces[i]);
	}
}

void mcalcLayerDetAndTrace(const IMAGE* sum, int size, int sampleStep, IMAGE* det, IMAGE* trace)
{
	const int NX = 3, NY = 3, NXY = 4;
	const int dx_s[NX][5] = { { 0, 2, 3, 7, 1 }, { 3, 2, 6, 7, -2 }, { 6, 2, 9, 7, 1 } };
	const int dy_s[NY][5] = { { 2, 0, 7, 3, 1 }, { 2, 3, 7, 6, -2 }, { 2, 6, 7, 9, 1 } };
	const int dxy_s[NXY][5] = { { 1, 1, 4, 4, 1 }, { 5, 1, 8, 4, -1 }, { 1, 5, 4, 8, -1 }, { 5, 5, 8, 8, 1 } };
	SurfHF Dx[NX], Dy[NY], Dxy[NXY];
	if (size > (sum->row - 1) || size > (sum->col - 1))
		return;
	mresizeHaarPattern(dx_s, Dx, NX, 9, size, sum->col);
	mresizeHaarPattern(dy_s, Dy, NY, 9, size, sum->col);
	mresizeHaarPattern(dxy_s, Dxy, NXY, 9, size, sum->col);


	int samples_i = 1 + (sum->row - 1 - size) / sampleStep;
	int samples_j = 1 + (sum->col - 1 - size) / sampleStep;


	int margin = (size / 2) / sampleStep;



	for (int i = 0; i < samples_i; i++)
	{


		const double* sum_ptr = (double*)(sum->gray + i*(sum->col)*sampleStep*sizeof(double));
		double* det_ptr = (double*)(det->gray + ((i + margin)*(det->col) + margin)*sizeof(double));
		double* trace_ptr = (double*)(trace->gray + ((i + margin)*(trace->col) + margin)*sizeof(double));
		for (int j = 0; j < samples_j; j++)
		{
			double dx = mcalcHarrPattern(sum_ptr, Dx, 3);
			double dy = mcalcHarrPattern(sum_ptr, Dy, 3);
			double dxy = mcalcHarrPattern(sum_ptr, Dxy, 4);
			sum_ptr += sampleStep;
			det_ptr[j] = dx*dy - 0.81f*dxy*dxy;
			trace_ptr[j] = dx + dy;
		}
	}



}
void mresizeHaarPattern(const int src[][5], SurfHF* dst, int n, int oldSize, int newSize, int widthStep)
{
	float ratio = (float)newSize / oldSize;
	for (int k = 0; k < n; k++)
	{
		int dx1 = int(ratio*src[k][0] + 0.5);
		int dy1 = int(ratio*src[k][1] + 0.5);
		int dx2 = int(ratio*src[k][2] + 0.5);
		int dy2 = int(ratio*src[k][3] + 0.5);
		dst[k].p0 = dy1*widthStep + dx1;
		dst[k].p1 = dy2*widthStep + dx1;
		dst[k].p2 = dy1*widthStep + dx2;
		dst[k].p3 = dy2*widthStep + dx2;
		dst[k].w = src[k][4] / ((float)(dx2 - dx1)*(dy2 - dy1));
	}
}

void mSURFFindInvoker(const IMAGE* sum, vector<IMAGE*>& dets, vector<IMAGE*>& traces, vector<int>& sizes,
	vector<int>& sampleSteps, vector<int>& middleIndices,
	vector<mKeyPoint>& keypoints, int nOctaveLayers, float hessianThreshold)
{
	int M = middleIndices.size();
	for (int i = 0; i < M; i++)
	{
		int layer = middleIndices[i];
		int octave = i / nOctaveLayers;
		mfindMaximaInLayer(sum, dets, traces, sizes, keypoints, octave, layer, hessianThreshold, sampleSteps[layer]);
	}

}






void mfindMaximaInLayer(const IMAGE* sum, const vector<IMAGE*> &dets, const vector<IMAGE*> &traces,
	const vector<int>&sizes, vector<mKeyPoint> &keypoints, int octave,
	int layer, float hessianThreshold, int sampleStep)
{
	const int NM = 1;
	const int dm[NM][5] = { { 0, 0, 9, 9, 1 } };
	SurfHF Dm;
	int size = sizes[layer];
	int layer_rows = (sum->row - 1) / sampleStep;
	int layer_cols = (sum->col - 1) / sampleStep;
	int margin = (sizes[layer + 1] / 2) / sampleStep + 1;


	int step = dets[layer]->col;
	for (int i = margin; i < layer_rows - margin; i++)
	{
		const double* det_ptr = (double*)(dets[layer]->gray + i*step*sizeof(double));
		const double* trace_ptr = (double*)(traces[layer]->gray + i*step*sizeof(double));
		for (int j = margin; j < layer_cols - margin; j++)
		{
			double val0 = det_ptr[j];
			if (val0 > hessianThreshold)
			{


				int sum_i = sampleStep*(i - (size / 2) / sampleStep);
				int sum_j = sampleStep*(j - (size / 2) / sampleStep);


				const double *det1 = (double*)(dets[layer - 1]->gray + (i*step + j)*sizeof(double));
				const double *det2 = (double*)(dets[layer]->gray + (i*step + j)*sizeof(double));
				const double *det3 = (double*)(dets[layer + 1]->gray + (i*step + j)*sizeof(double));
				float N9[3][9] = { { det1[-step - 1], det1[-step], det1[-step + 1],
					det1[-1], det1[0], det1[1],
					det1[step - 1], det1[step], det1[step + 1] },
					{ det2[-step - 1], det2[-step], det2[-step + 1],
					det2[-1], det2[0], det2[1],
					det2[step - 1], det2[step], det2[step + 1] },
					{ det3[-step - 1], det3[-step], det3[-step + 1],
					det3[-1], det3[0], det3[1],
					det3[step - 1], det3[step], det3[step + 1] } };

				if (val0 > N9[0][0] && val0 > N9[0][1] && val0 > N9[0][2] &&
					val0 > N9[0][3] && val0 > N9[0][4] && val0 > N9[0][5] &&
					val0 > N9[0][6] && val0 > N9[0][7] && val0 > N9[0][8] &&
					val0 > N9[1][0] && val0 > N9[1][1] && val0 > N9[1][2] &&
					val0 > N9[1][3] && val0 > N9[1][5] &&
					val0 > N9[1][6] && val0 > N9[1][7] && val0 > N9[1][8] &&
					val0 > N9[2][0] && val0 > N9[2][1] && val0 > N9[2][2] &&
					val0 > N9[2][3] && val0 > N9[2][4] && val0 > N9[2][5] &&
					val0 > N9[2][6] && val0 > N9[2][7] && val0 > N9[2][8])
				{


					float center_i = sum_i + (size - 1)*0.5f;
					float center_j = sum_j + (size - 1)*0.5f;
					mKeyPoint kpt(center_j, center_i, (float)sizes[layer],
						-1, val0, octave, (trace_ptr[j]) > 0);


					int ds = size - sizes[layer - 1];
					int interp_ok = interpolateKeypoint(N9, sampleStep, sampleStep, ds, kpt);


					if (interp_ok)
					{



						keypoints.push_back(kpt);
					}
				}
			}
		}
	}
}

static int interpolateKeypoint(float N9[3][9], int dx, int dy, int ds, mKeyPoint& kpt)
{

	double B[3] = { -(N9[1][5] - N9[1][3]) / 2,
		-(N9[1][7] - N9[1][1]) / 2,
		-(N9[2][4] - N9[0][4]) / 2 };
	double A[3][3] = { { N9[1][3] - 2 * N9[1][4] + N9[1][5],
		(N9[1][8] - N9[1][6] - N9[1][2] + N9[1][0]) / 4,
		(N9[2][5] - N9[2][3] - N9[0][5] + N9[0][3]) / 4 },
		{ (N9[1][8] - N9[1][6] - N9[1][2] + N9[1][0]) / 4,
		N9[1][1] - 2 * N9[1][4] + N9[1][7],
		(N9[2][7] - N9[2][1] - N9[0][7] + N9[0][1]) / 4 },
		{ (N9[2][5] - N9[2][3] - N9[0][5] + N9[0][3]) / 4,
		(N9[2][7] - N9[2][1] - N9[0][7] + N9[0][1]) / 4,
		N9[0][4] - 2 * N9[1][4] + N9[2][4] } };
	double x[3];
	double *temp = (double*)calloc(3, sizeof(double));
	if (Inv_matrix(A[0], 3))
	{
		temp = MatrixMultiply(A[0], 3, 3, B, 1);
		x[0] = temp[0];
		x[1] = temp[1];
		x[2] = temp[2];
	}
	else
		return 0;

	free(temp);

	bool ok = (x[0] != 0 || x[1] != 0 || x[2] != 0) &&
		std::abs(x[0]) <= 1 && std::abs(x[1]) <= 1 && std::abs(x[2]) <= 1;
	if (ok)
	{
		kpt.pt.x += x[0] * dx;
		kpt.pt.y += x[1] * dy;
		kpt.size = (float)(int)(kpt.size + x[2] * ds + 0.5);
	}
	return ok;
}


void mSURFInvoker(const IMAGE* img, const IMAGE* sum, vector<mKeyPoint> & keypoints, IMAGE* descriptors)
{
	enum { ORI_RADIUS = 6, ORI_WIN = 60, PATCH_SZ = 20 };
	const int nOriSampleBound = (2 * ORI_RADIUS + 1)*(2 * ORI_RADIUS + 1);
	vector<Point2F> apt;
	vector<float> aptw;
	vector<float> DW;
	apt.resize(nOriSampleBound);
	aptw.resize(nOriSampleBound);
	DW.resize(PATCH_SZ*PATCH_SZ);
	float SURF_ORI_SIGMA = 2.5f;

	double* G_ori = mgetGaussianKernel(2 * ORI_RADIUS + 1, SURF_ORI_SIGMA);

	int nOriSamples = 0;
	for (int i = -ORI_RADIUS; i <= ORI_RADIUS; i++)
	{
		for (int j = -ORI_RADIUS; j <= ORI_RADIUS; j++)
		{
			if (i*i + j*j <= ORI_RADIUS*ORI_RADIUS)
			{

				apt[nOriSamples].x = i;
				apt[nOriSamples].y = j;

				aptw[nOriSamples++] = G_ori[i + ORI_RADIUS] * G_ori[j + ORI_RADIUS];
			}
		}
	}
	int SURF_DESC_SIGMA = 3.3f;

	double* G_desc = mgetGaussianKernel(PATCH_SZ, SURF_DESC_SIGMA);
	for (int i = 0; i < PATCH_SZ; i++)
	{
		for (int j = 0; j < PATCH_SZ; j++)
			DW[i*PATCH_SZ + j] = G_desc[i] * G_desc[j];
	}





	const int NX = 2, NY = 2;
	const int dx_s[NX][5] = { { 0, 0, 2, 4, -1 }, { 2, 0, 4, 4, 1 } };
	const int dy_s[NY][5] = { { 0, 0, 4, 2, 1 }, { 0, 2, 4, 4, -1 } };

	float X[nOriSampleBound], Y[nOriSampleBound], angle[nOriSampleBound];
	vector<vector<double>> mPATCH(PATCH_SZ + 1, vector<double>(PATCH_SZ + 1, 0));
	float DX[PATCH_SZ][PATCH_SZ], DY[PATCH_SZ][PATCH_SZ];
	int dsize = 64;
	int k, k1 = 0, k2 = keypoints.size();
	float maxSize = 0;
	for (k = k1; k < k2; k++)
	{
		maxSize = std::max(maxSize, keypoints[k].size);
	}
	int imaxSize = std::max((int)ceil((PATCH_SZ + 1)*maxSize*1.2f / 9.0f), 1);

	for (k = k1; k < k2; k++)
	{
		int i, j, kk, nangle;
		double* vec;
		SurfHF dx_t[NX], dy_t[NY];
		mKeyPoint& kp = keypoints[k];
		float size = kp.size;
		Point2F center = kp.pt;
		float s = size*1.2f / 9.0f;
		int grad_wav_size = 2 * round(2 * s);
		if (sum->row < grad_wav_size || sum->col < grad_wav_size)
		{


			kp.size = -1;
			continue;
		}
		float descriptor_dir = 360.f - 90.f;
		mresizeHaarPattern(dx_s, dx_t, NX, 4, grad_wav_size, sum->col);
		mresizeHaarPattern(dy_s, dy_t, NY, 4, grad_wav_size, sum->col);
		for (kk = 0, nangle = 0; kk < nOriSamples; kk++)
		{
			int x = round(center.x + apt[kk].x*s - (float)(grad_wav_size - 1) / 2);
			int y = round(center.y + apt[kk].y*s - (float)(grad_wav_size - 1) / 2);
			if (y < 0 || y >= sum->row - grad_wav_size ||
				x < 0 || x >= sum->col - grad_wav_size)
				continue;
			const double* ptr = (double*)(sum->gray + (y*sum->col + x)*sizeof(double));
			float vx = mcalcHarrPattern(ptr, dx_t, 2);
			float vy = mcalcHarrPattern(ptr, dy_t, 2);
			X[nangle] = vx*aptw[kk];
			Y[nangle] = vy*aptw[kk];
			nangle++;
		}
		if (nangle == 0)
		{
			kp.size = -1;
			continue;
		}


		for (int i = 0; i < nangle; i++)
		{
			float temp = atan2(Y[i], X[i])*(180 / PI);
			if (temp < 0)
				angle[i] = temp + 360;
			else
				angle[i] = temp;

		}

		float bestx = 0, besty = 0, descriptor_mod = 0;
		int SURF_ORI_SEARCH_INC = 5;
		for (i = 0; i < 360; i += SURF_ORI_SEARCH_INC)
		{
			float sumx = 0, sumy = 0, temp_mod;
			for (j = 0; j < nangle; j++)
			{
				int d = std::abs(round(angle[j]) - i);
				if (d < ORI_WIN / 2 || d > 360 - ORI_WIN / 2)
				{
					sumx += X[j];
					sumy += Y[j];
				}
			}
			temp_mod = sumx*sumx + sumy*sumy;
			if (temp_mod > descriptor_mod)
			{
				descriptor_mod = temp_mod;
				bestx = sumx;
				besty = sumy;
			}
		}
		descriptor_dir = atan2(-besty, bestx);

		kp.angle = descriptor_dir;


		int win_size = (int)((PATCH_SZ + 1)*s);
		vector<vector<double> > mwin(win_size, vector<double>(win_size));
		descriptor_dir *= (float)(PI / 180);
		float sin_dir = -std::sin(descriptor_dir);
		float cos_dir = std::cos(descriptor_dir);




		float win_offset = -(float)(win_size - 1) / 2;
		float start_x = center.x + win_offset*cos_dir + win_offset*sin_dir;
		float start_y = center.y - win_offset*sin_dir + win_offset*cos_dir;
		int ncols1 = img->col - 1, nrows1 = img->row - 1;
		size_t imgstep = img->col;
		for (i = 0; i < win_size; i++, start_x += sin_dir, start_y += cos_dir)
		{
			double pixel_x = start_x;
			double pixel_y = start_y;
			for (j = 0; j < win_size; j++, pixel_x += cos_dir, pixel_y -= sin_dir)
			{
				int ix = floor(pixel_x), iy = floor(pixel_y);
				if ((unsigned)ix < (unsigned)ncols1 &&
					(unsigned)iy < (unsigned)nrows1)
				{
					float a = (float)(pixel_x - ix), b = (float)(pixel_y - iy);
					const double* imgptr = (double*)(img->gray + (iy*img->col + ix)*sizeof(double));
					mwin[i][j] = (double)
						round(imgptr[0] * (1.f - a)*(1.f - b) +
						imgptr[1] * a*(1.f - b) +
						imgptr[imgstep] * (1.f - a)*b +
						imgptr[imgstep + 1] * a*b);
				}
				else
				{
					int x = std::min(std::max((int)round(pixel_x), 0), ncols1);
					int y = std::min(std::max((int)round(pixel_y), 0), nrows1);
					mwin[i][j] = (double)(img->gray[y*img->col + x]);
				}
			}
		}

		resizeVV(mwin, mPATCH, 0);
		for (i = 0; i < PATCH_SZ; i++)
		for (j = 0; j < PATCH_SZ; j++)
		{
			float dw = DW[i*PATCH_SZ + j];
			float vx = (mPATCH[i][j + 1] - mPATCH[i][j] + mPATCH[i + 1][j + 1] - mPATCH[i + 1][j])*dw;
			float vy = (mPATCH[i + 1][j] - mPATCH[i][j] + mPATCH[i + 1][j + 1] - mPATCH[i][j + 1])*dw;
			DX[i][j] = vx;
			DY[i][j] = vy;
		}
		vec = (double*)(descriptors->gray + (k*descriptors->col)*sizeof(double));
		for (kk = 0; kk < dsize; kk++)
			vec[kk] = 0;
		double square_mag = 0;


		for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			for (int y = i * 5; y < i * 5 + 5; y++)
			{
				for (int x = j * 5; x < j * 5 + 5; x++)
				{
					float tx = DX[y][x], ty = DY[y][x];
					vec[0] += tx; vec[1] += ty;
					vec[2] += (float)fabs(tx);
					vec[3] += (float)fabs(ty);
				}
			}
			for (kk = 0; kk < 4; kk++)
				square_mag += vec[kk] * vec[kk];
			vec += 4;
		}
		vec = (double*)(descriptors->gray + (k*descriptors->col)*sizeof(double));
		float scale = (float)(1. / (sqrt(square_mag) + DBL_EPSILON));
		for (kk = 0; kk < dsize; kk++)
			vec[kk] *= scale;
	}



}




