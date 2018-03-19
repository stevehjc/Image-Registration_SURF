
#ifndef _MSURF_H_0123456_
#define _MSURF_H_0123456_

#include"proc.h"
#include"math_m.h"


struct SurfHF
{
	int p0, p1, p2, p3;
	float w;

	SurfHF() : p0(0), p1(0), p2(0), p3(0), w(0) {}
};


class Point2F
{
public:
	float x;
	float y;
	Point2F() :x(0), y(0) {}
	Point2F(float _x, float _y)
		:x(_x), y(_y)
	{}
	//float setX(float _x)
	//{
	//	x = _x;
	//}


};



class mKeyPoint
{
public:
	//mKeyPoint() : pt(0, 0), size(0), angle(-1), response(0), octave(0), class_id(-1) {};

	mKeyPoint(Point2F _pt, float _size, float _angle = -1,
		float _response = 0, int _octave = 0, int _class_id = -1)
		: pt(_pt), size(_size), angle(_angle),
		response(_response), octave(_octave), class_id(_class_id) {}

	mKeyPoint(float x, float y, float _size, float _angle = -1,
		float _response = 0, int _octave = 0, int _class_id = -1)
		: pt(x, y), size(_size), angle(_angle),
		response(_response), octave(_octave), class_id(_class_id) {}

	Point2F pt;
	float size;
	float angle;
	float response;
	int octave;
	int class_id;

	//~mKeyPoint();


};
struct mKeypointGreater
{
	inline bool operator()(const mKeyPoint& kp1, const mKeyPoint& kp2) const
	{
		if (kp1.response > kp2.response) return true;
		if (kp1.response < kp2.response) return false;
		if (kp1.size > kp2.size) return true;
		if (kp1.size < kp2.size) return false;
		if (kp1.octave > kp2.octave) return true;
		if (kp1.octave < kp2.octave) return false;
		if (kp1.pt.y < kp2.pt.y) return false;
		if (kp1.pt.y > kp2.pt.y) return true;
		return kp1.pt.x < kp2.pt.x;
	}
};

class mDMatch
{
public:
	mDMatch() : queryIdx(-1), trainIdx(-1), imgIdx(-1), distance(DBL_MAX) {}
	mDMatch(int _queryIdx, int _trainIdx, double _distance) :
		queryIdx(_queryIdx), trainIdx(_trainIdx), imgIdx(-1), distance(_distance) {}
	mDMatch(int _queryIdx, int _trainIdx, int _imgIdx, double _distance) :
		queryIdx(_queryIdx), trainIdx(_trainIdx), imgIdx(_imgIdx), distance(_distance) {}

	int queryIdx; // query descriptor index
	int trainIdx; // train descriptor index
	int imgIdx;   // train image index

	double distance;

	// less is better
	bool operator<(const mDMatch &m) const
	{
		return distance < m.distance;
	}
};
struct mathcePointsBetter
{
	inline bool operator()(const mDMatch& mp1, const mDMatch& mp2) const
	{

		if (mp1.distance < mp2.distance) return true;
		if (mp1.distance > mp2.distance) return false;
		//if (mp1.queryIdx < mp2.queryIdx) return true;
		//if (mp1.queryIdx > mp2.queryIdx) return false;
		/*if (mp1.octave > mp2.octave) return true;
		if (mp1.octave < mp2.octave) return false;
		if (mp1.pt.y < mp2.pt.y) return false;
		if (mp1.pt.y > mp2.pt.y) return true;*/
		return mp1.queryIdx < mp2.queryIdx;
	}
};


inline float mcalcHarrPattern(const double* origin, const SurfHF* f, int n)
{
	double d = 0;
	for (int k = 0; k < n; k++)
		d += (origin[f[k].p0] + origin[f[k].p3] - origin[f[k].p1] - origin[f[k].p2])*f[k].w;
	return (float)d;
}




IMAGE* toIMAGE(double* src, int row, int col);
IMAGE* CreatImage(int row, int col);
void DeleteImage(IMAGE * image);
void SURF_detect(IMAGE* srcImg, IMAGE* intergalSum, vector<mKeyPoint>&keypoints, float hessianThreshold);
void SURF_descriptor(IMAGE* srcImg, IMAGE* intergalSum, vector<mKeyPoint>& keypoints, IMAGE* descriptors);
void mFastHessianDetector(const IMAGE* sum, vector<mKeyPoint> &keypoints,
	int nOctaves, int nOctaveLayers, float hessianThreshold);
void mSURFBuildInvoker(const IMAGE* sum, vector<int>&sizes, vector<int>& sampleSteps, vector<IMAGE*>& dests, vector<IMAGE*>traces);
void mSURFFindInvoker(const IMAGE* sum, vector<IMAGE*>& dets, vector<IMAGE*>& traces, vector<int>& sizes,
	vector<int>& sampleSteps, vector<int>& middleIndices, vector<mKeyPoint>& keypoints, int nOctaveLayers, float hessianThreshold);
void mSURFInvoker(const IMAGE* img, const IMAGE* sum, vector<mKeyPoint> & keypoints, IMAGE* descriptors);
void intergal(IMAGE *img, IMAGE* sum);
static int interpolateKeypoint(float N9[3][9], int dx, int dy, int ds, mKeyPoint& kpt);
void mcalcLayerDetAndTrace(const IMAGE* sum, int size, int sampleStep, IMAGE* det, IMAGE* trace);
void mresizeHaarPattern(const int src[][5], SurfHF* dst, int n, int oldSize, int newSize, int widthStep);
void mfindMaximaInLayer(const IMAGE* sum, const vector<IMAGE*> &dets, const vector<IMAGE*> &traces,
	const vector<int>&sizes, vector<mKeyPoint> &keypoints, int octave,
	int layer, float hessianThreshold, int sampleStep);

#endif