#pragma once
#include <vector>
#include <ipp.h>
#include <mkl.h>
#include "ViType.h"

#define PI 3.14159265386
#define fitting_sigma 20

using namespace std;

enum ArcDirection//圆弧方向
{
	CLOCKWISE,//顺时针
	ANTICLOCLWISE//逆时针
};
typedef struct
{
	IMG_RCOORD CirCen;//圆心
	IMG_REAL Radius;//半径
}StructCircle;//圆

typedef struct
{
	IMG_RCOORD CirCen;//圆心
	IMG_REAL Radius;//半径
	IMG_RCOORD startPoint;//起点
	IMG_RCOORD endPoint;//终点
	IMG_REAL startAngle;
	IMG_REAL sweepAngle;
	ArcDirection direct;//方向
}StructArc;//圆弧

typedef struct tagEllipseStruct		//椭圆结构
{
	Ipp32f a;
	Ipp32f b;
	Ipp32f c;
	Ipp32f d;
	Ipp32f e;
	Ipp32f f;
	Ipp32f center_x;	//中心x坐标
	Ipp32f center_y;	//中心y坐标
	Ipp32f axis_a;		//长轴？
	Ipp32f axis_b;		//短轴？
	Ipp32f rotateTheta;	//旋转角度
}EllipseStruct;

typedef struct tagRegionPeakPoint		//连通域四个极值点
{
	int xMin;
	int xMax;
	int yMin;
	int yMax;
	tagRegionPeakPoint()
	{
		xMin = 0;
		xMax = 0;
		yMin = 0;
		yMax = 0;
	}

}RegionPeakPoint;

typedef struct tagRegionWidthHeight  //连通域宽高
{
	int width;
	int height;
	tagRegionWidthHeight()
	{
		width = 0;
		height = 0;
	}
}RegionWidthHeight;

typedef struct
{
	IMG_COORD xyInteger; //像素点
	IMG_RCOORD xyDecimal;//亚像素点
	int gradient;
	float angle;
}edgeInformation;//边缘点

class CVisAlg2DBase
{
public:
	CVisAlg2DBase();
	~CVisAlg2DBase();

	
	///////////////////Pyramid////////////////////////
	//Author:Liu Ping
	int VisPyramid(IMG_UBBUF src, unsigned char * pDst, int & pyramid_width, int & pyramid_height, int level);

	//Author:Shen Jiancheng
	IppStatus VisPyramid2(Ipp8u * pSrc, IppiSize roiSize, IppiPyramid *& pPyrStruct, Ipp8u **& pPyrImage, int level);



	/////////////////// Filter/////////////////////////
	
	//Author:  Jiang He/20170407
	int VisFilterMean(const unsigned char *src, const int srcHeight, const int srcWidth, unsigned char *dst, const unsigned char kernelSize, unsigned int divisor);

	int VisFilterGaussian(const unsigned char *src, const int srcHeight, const int srcWidth, unsigned char *dst, const unsigned char winWidth);
	
	int VisFilterMedian(const unsigned char *src, const int srcHeight, const int srcWidth, unsigned char *dst, const unsigned char winWidth);

	
	////////////////// Fitting////////////////////////////

	// Author: 申健成/20170227
	int VisFittingCircular(const float *point_pos, const int m, StructCircle &circular_fit, const int iteration_times);
	int VisFittingArc(const float *point_pos, const int m, StructArc &arc_fit, const int iteration_times);

	
	//Author:申健成/詹铭泽/20170412
	bool VisFitting_Ellipse(Ipp32f *xy, int m, EllipseStruct &resu);

	//Author:申健成/詹铭泽/20170412
	bool  VisFitting_Ellipse_Str(Ipp32f *xy, int m, EllipseStruct &resu);
	IMG_VVOID inverse(IMG_LREAL* A, IMG_INT N);

	
	// Author:Tan Ling/20170412
	IMG_INT GNERAL_LINE_FITTING(vector<IMG_RCOORD>PointCor, IMG_REAL &Slope, IMG_REAL &B, IMG_REAL Sigma, IMG_INT FLAG);
	


	//////////////// Histogram related////////////////
	//Author:Liu Ping
	int VisHistogram(IMG_UBBUF src, Ipp32u * pHist, int nBins);

	


	////////////Segment////////////
	int VisSegmentManual(const unsigned char * pSrc, const unsigned int nW, const unsigned int nH, unsigned char * pDst, const unsigned char nT);

	int VisSegmentOtsu(const unsigned char * pSrc, const unsigned int nW, const unsigned int nH, unsigned char * pDst);

	int VisCalcOtsu(const unsigned char * pSrc, const unsigned int nW, const unsigned int nH, unsigned char & nThres);

	int VisSegmentDynamic(const unsigned char * pSrc, const unsigned int nW, const unsigned int nH, unsigned char * pDst, const unsigned int avgWinWidth);

	int VisSegmentCntDots(unsigned char * pSrc, const unsigned int nW, const unsigned int nH, unsigned char * pDst, const double fWhiteRatio);

	//Author:Huang Yige
	void VisGetBlobThreshold(Ipp8u *srcRoi, IppiSize Roi, int BlobThreshold);

	///////////////////Morphology//////////////////

	int VisMorphologyErode(unsigned char * src1, int srcWidth, int srcHeight, unsigned char * dst);

	int VisMorphologyDilation(unsigned char * src1, int srcWidth, int srcHeight, unsigned char * dst);


	int VisMorphologyOpen(IMG_UBBUF src, IMG_UBBUF dst);

	int VisMorphologyClose(IMG_UBBUF src, IMG_UBBUF dst);

	
	// Author: 詹铭泽/20170412
	int getRegionInfo(IMG_UBBUF src,int regionNum,int *regionArea,RegionPeakPoint *regionPeakPoint,RegionWidthHeight *regionWH,
	IppiPoint *regionGraCenPoint);

	int VisHoleFill(unsigned char * srcSeg, int srcWidth, int srcHeight, unsigned char * dst);
	
	//Author:Huang Yige
	void VisLabelMarker(Ipp8u *srcRoi, IppiSize Roi, int &markersNum);

	//Author:Huang Yige
	void VisMoment(Ipp8u *srcRoi, IppiSize Roi, double &hu);

	
	


	///////// Gradient///////////Sobel等
	
	int GradientCompute(unsigned char * src, int srcWidth, int srcHeight, int rowNum, int colNum, int * dst, IMG_COORD * edgePoint);
    int SobelFilter_8u16s_C1_5x5or3x3(IMG_UBYTE * pSrc, IppiSize roiSize, Ipp16s * pDst, Ipp32f * pAngle, int flag);
	
	////Edge////////
	//Canny、own
	int GradientCompute(unsigned char * src, int srcWidth, int srcHeight, int rowNum, int colNum, float threshold, int * dst, vector<edgeInformation>& edgePoint);


	IMG_INT VisEdge_detection(IMG_UBYTE * srcRoi, IMG_SIZE roiSize, int threshold, IMG_WORD * dstRoi, IMG_UBYTE * dstRoiE, edgeInformation *& edgeArray, IMG_INT & eNum);

	
	
private:
	IMG_REAL getIntraClassVariance(Ipp16s * src, int srcRows, int srcCols, int & varTh);

};

