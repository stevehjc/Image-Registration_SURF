#ifndef _PROC_H_0123456_
#define _PROC_H_0123456_

#include<vector>
#include<iostream>
#include<algorithm>


using namespace std;


#ifndef _BITMAPFILEHEADER_
#define _BITMAPFILEHEADER_
typedef struct  _tagBITMAPFILEHEADER{
	unsigned short    bfType;
	unsigned long     bfSize;
	unsigned short    bfReserved1;
	unsigned short    bfReserved2;
	unsigned long     bfOffBits;
}_BITMAPFILEHEADER;
#endif

#ifndef _BITMAPINFOHEADER_
#define _BITMAPINFOHEADER_
typedef struct _tagBITMAPINFOHEADER{
	unsigned long       biSize;
	signed long         biWidth;
	signed long         biHeight;
	unsigned short      biPlanes;
	unsigned short      biBitCount;
	unsigned long       biCompression;
	unsigned long       biSizeImage;
	signed long         biXPelsPerMeter;
	signed long         biYPelsPerMeter;
	unsigned long       biClrUsed;
	unsigned long       biClrImportant;
} _BITMAPINFOHEADER;
#endif

#ifndef _RGBQUAD_
#define _RGBQUAD_
typedef struct _tagRGBQUAD {
	unsigned char    rgbBlue;
	unsigned char    rgbGreen;
	unsigned char    rgbRed;
	unsigned char    rgbReserved;
} _RGBQUAD;
#endif


class Point2D
{
public:
	double x;
	double y;
	Point2D(double _x, double _y)
		:x(_x), y(_y){}
	Point2D() :x(0), y(0) {}
};

typedef struct tagIMAGE {
	int row;        /* Number of rows (dy) */
	int col;        /* Number of columns (dx) */
	int imagesize;   /* Size allocated (in bytes) for the gray plane */
	unsigned char *gray;     /* The Gray level plane (may be NULL) */

	float scale;     /* Scale of the picture (should be 1 for original pict.) */
	char cmt[128]; /* Comments */
	char name[64]; /* Name of the image */

	/* Defines the signifiant part of the picture : */
	int firstcol;    /* index of the first col not affected by left side effect*/
	int lastcol;     /* index of the last col not affected by right side effect*/
	int firstrow;    /* index of the first row not aff. by upper side effect */
	int lastrow;     /* index of the last row not aff. by lower side effect */

	/* For use in Movies only */
	struct tagIMAGE *previous; /* Pointer to the previous image (may be NULL) */
	struct tagIMAGE *next; /* Pointer to the next image (may be NULL) */

} IMAGE;

void * fspace_1d(int col, int length);



double calcDistance(double* p1, double* p2, int n);

int WriteBMPFile(unsigned char *image, int Row, int Col, char *FileName);
int WriteBMPFile(double *image, int Row, int Col, char *FileName);
unsigned char * ReadBMPFile(int *Row, int *Col, char *FileName);

bool mSURF_Detection(double *srcImg1, double *srcImg2, int Rows, int Cols,
	vector<Point2D> &matchpoints1, vector<Point2D> &matchpoints2, int matches_num);
void mWarpPerspective(double* srcImg1, double* imageTransform1, vector<double>para, int Rows, int Cols);
double* getPerspectiveTransformIata(vector<Point2D> src, vector<Point2D> dst);
vector<double> getPerspectiveTransformLSM(vector<Point2D> src, vector<Point2D> dst);
vector<double> getPerspectiveTransformLSM2(vector<Point2D> src, vector<Point2D> dst);
vector<double> getRANSAC2(vector<Point2D>p1, vector<Point2D>p2, double threshold, int count);
void resizeVV(vector<vector<double>> src, vector<vector<double>>& dst, int interpolation);





#endif