//#include "stdafx.h"
#include "VisAlg2DBase.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <stdexcept>
#include "BmErrorCode.h"
#include <Windows.h>

#define fitting_sigma 10
#define PI 3.14159265386

using namespace std;

CVisAlg2DBase::CVisAlg2DBase()
{
}


CVisAlg2DBase::~CVisAlg2DBase()
{
}


///////Pyramid
//Author:Liu Ping
int CVisAlg2DBase::VisPyramid(IMG_UBBUF src, unsigned char * pDst, int & pyramid_width, int & pyramid_height, int level)
{
	__try
	{
		IppStatus   status = ippStsNoErr;
		//ofstream outfile("pyramidData.txt");

		IppiPyramid*  pPyrStruct = NULL;
		unsigned char **  pPyrImage = NULL;

		//init var
		//double sigma = 3;
		unsigned char* pSrc = src.ptr;
		IMG_SIZE roiSize = src.size;

		if (!pSrc) { status = ippStsNoMemErr; goto exit; }

		int srcStep = roiSize.width * sizeof(unsigned char);
		float      rate = 2.f;                  // Neighbour levels ratio
		signed short kernel[3] = { 1,1,1 };		// Separable symmetric kernel of odd length

												//signed short *kernel = (signed short *)malloc(3 * sizeof(signed short));
												//__GetGaussianKernel_dim1(kernel, 3, sigma);	// preserved

		int pyrBufferSize = 0;
		int pyrStructSize = 0;
		unsigned char   *pPyrBuffer = NULL;
		unsigned char   *pPyrStrBuffer = NULL;

		int      pyrLStateSize = 0;
		int      pyrLBufferSize = 0;
		unsigned char   *pPyrLStateBuf = NULL;
		unsigned char   *pPyrLBuffer = NULL;

		// Computes the temporary work buffer size
		status = ippiPyramidGetSize(&pyrStructSize, &pyrBufferSize, level, { roiSize.width,roiSize.height }, rate);
		if (status < 0) goto exit;

		//pPyrBuffer = ippsMalloc_8u(pyrBufferSize);
		//pPyrStrBuffer = ippsMalloc_8u(pyrStructSize);
		pPyrBuffer = (unsigned char*)malloc(pyrBufferSize * sizeof(unsigned char));
		pPyrStrBuffer = (unsigned char*)malloc(pyrStructSize * sizeof(unsigned char));	//not pop
		if ((pyrBufferSize && !pPyrBuffer) || (pyrStructSize && !pPyrStrBuffer)) {
			status = ippStsNoMemErr; goto exit_proc;
		}

		// Initializes Gaussian structure for pyramids
		//pPyrStruct = (IppiPyramid*)malloc(level * sizeof(IppiPyramid));	
		status = ippiPyramidInit(&pPyrStruct, level, { roiSize.width,roiSize.height }, rate, pPyrStrBuffer, pPyrBuffer);
		if (status < 0) goto exit_proc;

		// ????????????????Correct maximum scale level 
		level = pPyrStruct->level;

		// Allocate structures to calculate pyramid layers
		status = ippiPyramidLayerDownGetSize_8u_C1R({ roiSize.width,roiSize.height }, rate, 3, &pyrLStateSize, &pyrLBufferSize);
		if (status < 0) goto exit_proc;

		//pPyrLStateBuf = ippsMalloc_8u(pyrLStateSize);
		//pPyrLBuffer = ippsMalloc_8u(pyrLBufferSize);
		pPyrLStateBuf = (unsigned char*)malloc(pyrLStateSize * sizeof(unsigned char));
		pPyrLBuffer = (unsigned char*)malloc(pyrLBufferSize * sizeof(unsigned char));
		if ((pyrLStateSize && !pPyrLStateBuf) || (pyrLBufferSize && !pPyrLBuffer)) { status = ippStsNoMemErr; goto exit; }

		// Initialize the structure for creating a lower pyramid layer
		status = ippiPyramidLayerDownInit_8u_C1R((IppiPyramidDownState_8u_C1R**)&pPyrStruct->pState, { roiSize.width,roiSize.height }, rate, kernel, 3, IPPI_INTER_LINEAR, pPyrLStateBuf, pPyrLBuffer);
		if (status < 0) goto exit_proc;

		// Allocate pyramid layers
		pPyrImage = pPyrStruct->pImage;
		pPyrImage[0] = pSrc;
		pPyrStruct->pStep[0] = srcStep;

		for (int i = 1; i <= level; i++)
		{
			//pPyrImage[i] = ippiMalloc_8u_C1(pPyrStruct->pRoi[i].width, pPyrStruct->pRoi[i].height, &pPyrStruct->pStep[i]);
			pPyrImage[i] = (unsigned char*)malloc((pPyrStruct->pRoi[i].width) * (pPyrStruct->pRoi[i].height) * sizeof(unsigned char));
			pPyrStruct->pStep[i] = (pPyrStruct->pRoi[i].width) * sizeof(unsigned char);
			if (!pPyrImage[i]) { status = ippStsNoMemErr; goto exit_proc; }
		}

		// Perform downsampling of the image with 5x5 Gaussian kernel
		for (int i = 1; i <= level; i++)
		{
			status = ippiPyramidLayerDown_8u_C1R(pPyrImage[i - 1], pPyrStruct->pStep[i - 1], pPyrStruct->pRoi[i - 1], pPyrImage[i], pPyrStruct->pStep[i], pPyrStruct->pRoi[i], (IppiPyramidDownState_8u_C1R*)pPyrStruct->pState);
			if (status < 0) goto exit_proc;

		}

		for (int i = 0; i < pPyrStruct->pRoi[level].height; i++)
		{
			for (int j = 0; j < pPyrStruct->pRoi[level].width; j++)
			{
				pDst[i * pPyrStruct->pRoi[level].width + j] = pPyrImage[level][i * pPyrStruct->pRoi[level].width + j];
			}
		}
		pyramid_width = pPyrStruct->pRoi[level].width;
		pyramid_height = pPyrStruct->pRoi[level].height;

		//	test //
		//for (int i = 0; i < pPyrStruct->pRoi[level].height; i++)
		//{
		//	for (int j = 0; j < pPyrStruct->pRoi[level].width; j++)
		//	{
		//		outfile << (float)pPyrImage[level][i * pPyrStruct->pRoi[level].width + j] << " ";
		//	}
		//	outfile << endl;
		//}
		//outfile.close();



	exit_proc:
		for (int i = 1; i <= level; i++)
			free(pPyrImage[i]);
		free(pPyrStrBuffer);
		free(pPyrBuffer);
		free(pPyrLBuffer);
		free(pPyrLStateBuf);
		//for (int i = 1; i <= level; i++)
		//	ippiFree(pPyrImage[i]);
		//ippiFree(pPyrLStateBuf);
		//ippiFree(pPyrLBuffer);
		//ippiFree(pPyrStrBuffer);
		//ippiFree(pPyrBuffer);
		goto exit;

	exit:
		//printf("pyramidStatus: %s\n", ippGetStatusString(status));
		if (!status)
		{
		}
		return status;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
		return _CODE_THROW;
	}
}


//Author:Shen Jiancheng
IppStatus CVisAlg2DBase::VisPyramid2(Ipp8u* pSrc, IppiSize roiSize, IppiPyramid *&pPyrStruct, Ipp8u **&pPyrImage, int level)
{
	//pool.Push();
	//char savefile[20];
	IppStatus   status = ippStsNoErr;
	//if (!pSrc) { status = ippStsNoMemErr; goto exit; }
	int srcStep = roiSize.width * sizeof(Ipp8u);
	Ipp32f      rate = 2.f;                  // Neighbour levels ratio
	Ipp16s kernel[3] = { 1,1,1 };		// Separable symmetric kernel of odd length
										//__GetGaussianKernel(kernel, 3, sigma);	// preserved

	int i = 0;

	int pyrBufferSize = 0;
	int pyrStructSize = 0;
	Ipp8u       *pPyrBuffer = NULL;
	Ipp8u       *pPyrStrBuffer = NULL;

	int      pyrLStateSize = 0;
	int      pyrLBufferSize = 0;
	Ipp8u   *pPyrLStateBuf = NULL;
	Ipp8u   *pPyrLBuffer = NULL;

	// Computes the temporary work buffer size
	status = ippiPyramidGetSize(&pyrStructSize, &pyrBufferSize, level, roiSize, rate);
	if (status !=0) goto exit;

	pPyrBuffer = (Ipp8u*)malloc(pyrBufferSize * sizeof(Ipp8u));
	pPyrStrBuffer = (Ipp8u*)malloc(pyrStructSize * sizeof(Ipp8u));
	/*if ((pyrBufferSize && !pPyrBuffer) || (pyrStructSize && !pPyrStrBuffer)) {
	status = ippStsNoMemErr; goto exit;
	}*/

	// Initializes Gaussian structure for pyramids
	//pPyrStruct = (IppiPyramid*)pool.Malloc(level * sizeof(IppiPyramid));	
	status = ippiPyramidInit(&pPyrStruct, level, roiSize, rate, pPyrStrBuffer, pPyrBuffer);
	if (status != 0) goto exit;

	// ????????????????Correct maximum scale level 
	level = pPyrStruct->level;

	// Allocate structures to calculate pyramid layers
	status = ippiPyramidLayerDownGetSize_8u_C1R(roiSize, rate, 3, &pyrLStateSize, &pyrLBufferSize);
	if (status != 0) goto exit;

	pPyrLStateBuf = (Ipp8u*)malloc(pyrLStateSize * sizeof(Ipp8u));
	pPyrLBuffer = (Ipp8u*)malloc(pyrLBufferSize * sizeof(Ipp8u));
	//if ((pyrLStateSize && !pPyrLStateBuf) || (pyrLBufferSize && !pPyrLBuffer)) { status = ippStsNoMemErr; goto exit; }

	// Initialize the structure for creating a lower pyramid layer
	status = ippiPyramidLayerDownInit_8u_C1R((IppiPyramidDownState_8u_C1R**)&pPyrStruct->pState, roiSize, rate, kernel, 3, IPPI_INTER_LINEAR, pPyrLStateBuf, pPyrLBuffer);
	if (status != 0) goto exit;

	// Allocate pyramid layers
	pPyrImage = pPyrStruct->pImage;
	pPyrImage[0] = pSrc;
	pPyrStruct->pStep[0] = srcStep;
	for (i = 1; i <= level; i++)
	{
		pPyrImage[i] = (Ipp8u*)malloc((pPyrStruct->pRoi[i].width) * (pPyrStruct->pRoi[i].height) * sizeof(Ipp8u));
		pPyrStruct->pStep[i] = (pPyrStruct->pRoi[i].width) * sizeof(Ipp8u);
		if (!pPyrImage[i]) { status = ippStsNoMemErr; goto exit; }
	}

	// Perform downsampling of the image with 5x5 Gaussian kernel
	for (i = 1; i <= level; i++)
	{
		status = ippiPyramidLayerDown_8u_C1R(pPyrImage[i - 1], pPyrStruct->pStep[i - 1], pPyrStruct->pRoi[i - 1], pPyrImage[i], pPyrStruct->pStep[i], pPyrStruct->pRoi[i], (IppiPyramidDownState_8u_C1R*)pPyrStruct->pState);
		if (status != 0) goto exit;

	}


	if (pPyrBuffer != NULL)
		free(pPyrBuffer);
	/*if (pPyrStrBuffer != NULL)
	free(pPyrStrBuffer);*/
	if (pPyrLBuffer != NULL)
		free(pPyrLBuffer);
	if (pPyrLStateBuf != NULL)
		free(pPyrLStateBuf);
	exit:	//pool.Pop();			//do not pop!!!!!!!!!!!!!!!!!
	return status;
}



/******************************  Filter  *****************************************/

//均值滤波
//VisMeanFilter功能说明：基于ipp卷积函数ippiConv_8u_C1R，
//滤波后的图像大小为(srcWidth - kernelSize + 1)*(srcHeight - kernelSize + 1)，边缘灰度值与原图相同。
//Input:
//unsigned char *src, int srcHeight, int srcWidth   输入图像
//unsigned int kernelSize  核的尺寸，大于等于3的奇数
//int divisor   可取为kernelSize*kernelSize，若大于该值即为对图像整体灰度值进行了调暗处理，若大于该值即为对图像进行了调亮处理。
//output
//unsigned char *dst 输出图像，与输入图像宽高一致。
//Return
// 0   正常
// -1  异常
//Author:  Jiang He/20170407
int CVisAlg2DBase::VisFilterMean(const unsigned char *src, const int srcHeight, const int srcWidth, unsigned char *dst, 
	const unsigned char kernelSize, unsigned int divisor)
{
	//如果kernelSize不是大于等于3的奇数，返回-1.
	if (kernelSize < 3)
	{
		return -1;
	}
	else
	{
		if (!(kernelSize % 2))
		{
			return -1;
		}
	}

	if (divisor == 0)
	{
		divisor = kernelSize*kernelSize;
	}

	unsigned char *src2;
	src2 = (Ipp8u*)malloc(kernelSize * kernelSize);// *sizeof(Ipp8u));
	for (int i = 0; i < kernelSize*kernelSize; i++)
	{
		src2[i] = 1;
	}


	IppStatus status = ippStsNoErr;
	IppiSize src1Size = { srcWidth ,srcHeight };
	int src1Step = srcWidth * sizeof(Ipp8u);

	IppiSize src2Size = { kernelSize,kernelSize };
	int src2Step = kernelSize * sizeof(Ipp8u);
	//IppEnum algType = (IppEnum)(ippAlgAuto | ippiROIFull | ippiNormNone);
	IppEnum algType = (IppEnum)(ippAlgAuto | ippiROIValid | ippiNormNone);
	int buffSize = 0;
	status = ippiConvGetBufferSize(src1Size, src2Size, ipp8u, 1, algType, &buffSize);
	if (status != ippStsNoErr)
	{
		free(src2);
		return -1;
	}


	//int divisor = 15;
	Ipp8u *buffer;
	buffer = (Ipp8u*)ippsMalloc_8u(buffSize);

	//unsigned char *dstt = (unsigned char*)malloc((srcWidth + kernelSize - 1)*(srcHeight + kernelSize - 1)*sizeof(unsigned char));
	//int dstStep = (srcWidth + kernelSize - 1)*sizeof(Ipp8u);
	unsigned char *dstt = (unsigned char*)malloc((srcWidth - kernelSize + 1)*(srcHeight - kernelSize + 1) * sizeof(unsigned char));
	int dstStep = (srcWidth - kernelSize + 1) * sizeof(Ipp8u);

	status = ippiConv_8u_C1R(src, src1Step, src1Size, src2, src2Step, src2Size, dstt, dstStep, kernelSize*kernelSize, algType, buffer);

	
	if (status != ippStsNoErr)
	{
		free(src2);
		ippsFree(buffer);
		free(dstt);
		return -1;
	}

	//for (int i = 0; i < srcHeight; i++)
	//{
	//	for (int j = 0; j < srcWidth; j++)
	//	{
	//		dst[j + i*srcWidth] = src[j + i*srcWidth];
	//	}
	//}
	memcpy(dst, src, srcWidth*srcHeight*sizeof(unsigned char));

	int hs = kernelSize / 2;
	unsigned int pos1 = hs*srcWidth;
	unsigned int pos2 = 0;
	int w1 = srcWidth - kernelSize + 1;
	for (int i = hs; i < srcHeight - hs; i++)
	{
		memcpy(&dst[pos1 + hs], &dstt[pos2], w1);
		//for (int j = hs; j < srcWidth - hs; j++)
		//{
		//	dst[j + pos1] = dstt[(j - hs) + pos2];
		//}
		pos1 += srcWidth;
		pos2 += w1;
	}

	free(src2);
	ippsFree(buffer);
	free(dstt);
	return 0;
}



//高斯滤波
//VisGaussianFilter功能说明：基于ipp高斯滤波函数ippiFilterGaussianBorder_8u_C1R，
//Input
//unsigned char *src, int srcHeight, int srcWidth   输入图像
//float sigma  Standard deviation of the Gaussian kernel.
//unsigned int winWidth  核的尺寸，大于等于3的奇数
//output
//unsigned char *dst 输出图像，与输入图像宽高一致。
//
//Return
// 0   正常
// -1  异常
//Author:  Jiang He/20170407
int CVisAlg2DBase::VisFilterGaussian(const unsigned char *src, const int srcHeight, const int srcWidth, unsigned char *dst, const unsigned char winWidth)
{
	//如果kernelSize不是大于等于3的奇数，返回-1.
	if (winWidth < 3)
	{
		return -1;
	}
	else
	{
		if (!(winWidth % 2))
		{
			return -1;
		}
	}

	IppiSize roiSize = { srcWidth,srcHeight };
	Ipp32u kernelSize = winWidth; //must be odd and greater or equal to 3.
	int tmpBufSize = 0, specSize = 0;
	Ipp32f sigma = float(kernelSize);
	IppiBorderType borderType = ippBorderRepl;

	ippcvFilterGaussianSpec *spec = NULL;
	Ipp8u *buffer;

	//1、get buffer size
	IppStatus status1 = ippiFilterGaussianGetBufferSize(roiSize, kernelSize, ipp8u, 1, &specSize, &tmpBufSize);

	if (status1 != ippStsNoErr)
	{
		return -1;
	}

	//2、init
	spec = (IppFilterGaussianSpec*)ippsMalloc_8u(specSize);
	buffer = ippsMalloc_8u(tmpBufSize);
	IppStatus status2 = ippiFilterGaussianInit(roiSize, kernelSize, sigma, borderType, ipp8u, 1, spec, buffer);

	if (status2 != ippStsNoErr)
	{
		if (buffer != NULL)
		{
			ippsFree(buffer);
		}
		if (spec != NULL)
		{
			ippsFree(spec);
		}
		return -1;
	}

	//3、filter
	int srcStep = srcWidth * sizeof(unsigned char);
	Ipp8u borderValue = 0;
	IppStatus status3 = ippiFilterGaussianBorder_8u_C1R(src, srcStep, dst, srcStep, roiSize, borderValue, spec, buffer);


	if (status3 != ippStsNoErr)
	{
		if (buffer != NULL)
		{
			ippsFree(buffer);
		}
		if (spec != NULL)
		{
			ippsFree(spec);
		}
		return -1;
	}


	if (buffer != NULL)
	{
		ippsFree(buffer);
	}
	if (spec != NULL)
	{
		ippsFree(spec);
	}
	return 0;
}

//中值滤波
//Input
//unsigned char *src, int srcHeight, int srcWidth   输入图像
//const unsigned char winWidth  核的尺寸，大于等于3的奇数
//output
//unsigned char *dst 输出图像，与输入图像宽高一致。
//
//Return
// 0   正常
// -1  异常
//Author:  Jiang He/20170407
int CVisAlg2DBase::VisFilterMedian(const unsigned char *src, const int srcHeight, const int srcWidth, unsigned char *dst, const unsigned char winWidth)
{
	//如果kernelSize不是大于等于3的奇数，返回 - 1.
		if (winWidth < 3)
		{
			return -1;
		}
		else
		{
			if (!(winWidth % 2))
			{
				return -1;
			}
		}


	IppStatus status;

	int srcStep = srcWidth * sizeof(unsigned char);

	IppiSize dstSize = { srcWidth - 1,srcHeight - 1 };
	IppiSize maskSize = { winWidth,winWidth };
	IppDataType dataType = ipp8u;
	IppiBorderType borderType = ippBorderRepl;
	Ipp8u borderValue = NULL;

	int buffSize = 0;
	Ipp8u *buffer = NULL;

	status = ippiFilterMedianBorderGetBufferSize(dstSize, maskSize, dataType, 1, &buffSize);
	if (status != ippStsNoErr) return -1;

	buffer = ippsMalloc_8u(buffSize);
	status = ippiFilterMedianBorder_8u_C1R(src + srcStep, srcStep, dst, srcStep, dstSize, maskSize, borderType, borderValue, buffer);

	ippsFree(buffer);
	if (status != ippStsNoErr) return -1;

	return 0;
}


/**********************************************  Fitting  *******************************************/

/**********************************************/
//VisFittingCircular,拟合圆
// Input:
//    *point_pos,输入点坐标，默认第一个为起点，最后一个为终止点
//     m，输入点对的个数
//	   iteration_times 迭代次数，输入0就行
// Output:
//    StructCircle, 圆心坐标和半径
// Return:
//     0 - 正常
//     1 - 点数少于3
//     2 - 无法拟合成圆
//     ...
// Author: 申健成/20170227
/**********************************************/
int CVisAlg2DBase::VisFittingCircular(const float *point_pos, const int m, StructCircle &circular_fit, const int iteration_times)
{
	int status = 0;
	if (m < 3)
	{
		printf("point number <3,can't fitting circular\n");
		return 1;
	}
	int tempxnum = 0, tempynum = 0;
	for (int i = 0; i < m * 2; i = i + 2)
	{
		if (point_pos[0] == point_pos[i])
		{
			tempxnum++;
		}
		if (point_pos[1] == point_pos[i + 1])
		{
			tempynum++;
		}
	}
	if (tempxnum == m || tempynum == m)
	{
		printf("It is a line,can't fitting circular\n");
		return 2;
	}

	float lastsum = 0, nowsum = 0;
	float *A;
	A = (float *)malloc(3 * sizeof(float));
	A[0] = A[1] = A[2] = 0;
	float *B, *C, *recordB, *recordC;
	B = (float *)malloc(3 * m * sizeof(float));
	C = (float *)malloc(m * sizeof(float));
	recordB = (float *)malloc(3 * m * sizeof(float));
	recordC = (float *)malloc(m * sizeof(float));
	float *error;
	float *error_weight;
	error_weight = (float *)malloc(m * sizeof(float));
	error = (float *)malloc(m * sizeof(float));
	for (int i = 0; i < m; i++)
	{
		error_weight[i] = 1.0;
		error[i] = 0;
	}
	int temp_time = 0;
	do {
		lastsum = nowsum;
		nowsum = 0;

		for (int i = 0, j = 0, k = 0; i < 3 * m; i = i + 3, j = j + 2, k++)
		{
			recordB[i] = B[i] = -2 * point_pos[j] * error_weight[k];
			recordB[i + 1] = B[i + 1] = -2 * point_pos[j + 1] * error_weight[k];
			recordB[i + 2] = B[i + 2] = 1 * error_weight[k];
			recordC[k] = C[k] = -(point_pos[j] * point_pos[j] + point_pos[j + 1] * point_pos[j + 1])* error_weight[k];
		}

		int info = LAPACKE_sgels(LAPACK_ROW_MAJOR, 'N', m, 3, 1, B, 3, C, 1);
		if (info != 0)
		{
			status = info;
			goto end;
		}
		A[0] = C[0];
		A[1] = C[1];
		A[2] = C[2];
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, 1, 3, 1, recordB, 3, A, 1, -1, recordC, 1);

		for (int i = 0; i < (m); i++)
		{
			error[i] = recordC[i] * recordC[i];
			nowsum = error[i] + nowsum;
		}
		for (int i = 0; i < m; i++)
		{
			error_weight[i] = float(1.0f / (1 + exp((error[i]) / (fitting_sigma*fitting_sigma))));
		}
		temp_time++;
	} while (fabs(nowsum - lastsum) > 0.0002 && temp_time < iteration_times);

	if (nowsum >= 1000000)
	{
		printf("Too big error\n");
		goto end;
		return 3;
	}
	circular_fit.CirCen.x = C[0];
	circular_fit.CirCen.y = C[1];
	circular_fit.Radius = (float)sqrt(C[0] * C[0] + C[1] * C[1] - C[2]);

end:
	if (error_weight != nullptr)
		free(error_weight);
	if (A != nullptr)
		free(A);
	if (B != nullptr)
		free(B);
	if (C != nullptr)
		free(C);
	if (recordB != nullptr)
		free(recordB);
	if (recordC != nullptr)
		free(recordC);
	if (error != nullptr)
		free(error);
	return status;
}

int CVisAlg2DBase::VisFittingArc(const float *point_pos, const int m, StructArc &arc_fit, const int iteration_times)
{
	int status = 0;
	if (m < 3)
	{
		//printf("point number <3,can't fitting circular\n");
		return 1;
	}
	int tempxnum = 0, tempynum = 0;
	for (int i = 0; i < m * 2; i = i + 2)
	{
		if (point_pos[0] == point_pos[i])
		{
			tempxnum++;
		}
		if (point_pos[1] == point_pos[i + 1])
		{
			tempynum++;
		}
	}
	if (tempxnum == m || tempynum == m)
	{
		//printf("It is a line,can't fitting circular\n");
		return 2;
	}

	float lastsum = 0, nowsum = 0;
	float *A;
	A = (float *)malloc(3 * sizeof(float));
	A[0] = A[1] = A[2] = 0;
	float *B, *C, *recordB, *recordC;
	B = (float *)malloc(3 * m * sizeof(float));
	C = (float *)malloc(m * sizeof(float));
	recordB = (float *)malloc(3 * m * sizeof(float));
	recordC = (float *)malloc(m * sizeof(float));
	float *error;
	float *error_weight;
	error_weight = (float *)malloc(m * sizeof(float));
	error = (float *)malloc(m * sizeof(float));
	for (int i = 0; i < m; i++)
	{
		error_weight[i] = 1.0;
		error[i] = 0;
	}
	int temp_time = 0;
	do {
		lastsum = nowsum;
		nowsum = 0;

		for (int i = 0, j = 0, k = 0; i < 3 * m; i = i + 3, j = j + 2, k++)
		{
			recordB[i] = B[i] = -2 * point_pos[j] * error_weight[k];
			recordB[i + 1] = B[i + 1] = -2 * point_pos[j + 1] * error_weight[k];
			recordB[i + 2] = B[i + 2] = 1 * error_weight[k];
			recordC[k] = C[k] = -(point_pos[j] * point_pos[j] + point_pos[j + 1] * point_pos[j + 1])* error_weight[k];
		}

		int info = LAPACKE_sgels(LAPACK_ROW_MAJOR, 'N', m, 3, 1, B, 3, C, 1);
		if (info != 0)
		{
			status = info;
			goto end;
		}
		A[0] = C[0];
		A[1] = C[1];
		A[2] = C[2];
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, 1, 3, 1, recordB, 3, A, 1, -1, recordC, 1);

		for (int i = 0; i < (m); i++)
		{
			error[i] = recordC[i] * recordC[i];
			nowsum = error[i] + nowsum;
		}
		for (int i = 0; i < m; i++)
		{
			error_weight[i] = float(1.0f / (1 + exp((error[i]) / (fitting_sigma*fitting_sigma))));
		}
		temp_time++;
	} while (fabs(nowsum - lastsum) > 0.0002 && temp_time < iteration_times);

	if (nowsum >= 1000000)
	{
		//printf("Too big error\n");
		return 3;
	}
	arc_fit.CirCen.x = C[0];
	arc_fit.CirCen.y = C[1];
	arc_fit.Radius = (float)sqrt(C[0] * C[0] + C[1] * C[1] - C[2]);
	float r = arc_fit.Radius;


	//优劣弧判断 
	float d1 = (point_pos[0] - point_pos[2 * m - 2])*(point_pos[0] - point_pos[2 * m - 2]) +
		(point_pos[1] - point_pos[2 * m - 1])*(point_pos[1] - point_pos[2 * m - 1]);
	float d2 = (point_pos[0] - point_pos[m - 1])*(point_pos[0] - point_pos[m - 1]) +
		(point_pos[1] - point_pos[m])*(point_pos[1] - point_pos[m]);
	float d3 = (point_pos[2 * m - 2] - point_pos[m - 1])*(point_pos[2 * m - 2] - point_pos[m - 1]) +
		(point_pos[2 * m - 1] - point_pos[m])*(point_pos[2 * m - 1] - point_pos[m]);

	float costh = (d2 + d3 - d1) / (2 * sqrt(d2)*sqrt(d3));
	//costh=acos(costh)*180/PI;
	costh = 90 - costh * 90;

	int flag = 3;//yuan costh=180
	if (costh > 90 && costh < 180)
	{
		flag = 0;//劣弧
	}
	else if (costh == 90)
	{
		flag = 2;//半弧
	}
	else if (costh < 90)
	{
		flag = 1;
	}

	// 计算起始角度
	float angle_th = 0;
	float record_x = point_pos[0] - arc_fit.CirCen.x;
	float record_y = arc_fit.CirCen.y - point_pos[1];
	if (record_y >= 0)
	{
		angle_th = 360 - atan2(record_y, record_x) * 180 / PI;
	}
	else
	{
		angle_th = -atan2(record_y, record_x) * 180 / PI;//逆时针
	}
	arc_fit.startAngle = angle_th;

	float angle_th_end = 0;
	float record_x_end = point_pos[2 * m - 2] - arc_fit.CirCen.x;
	float record_y_end = arc_fit.CirCen.y - point_pos[2 * m - 1];
	if (record_y_end >= 0)
	{
		angle_th_end = 360 - atan2(record_y_end, record_x_end) * 180 / PI;
	}
	else
	{
		angle_th_end = -atan2(record_y_end, record_x_end) * 180 / PI;//逆时针
	}

	float sweepAngle = 0;

	{
		if (angle_th_end > angle_th)
		{
			sweepAngle = angle_th_end - angle_th;
			if (sweepAngle < 180)
			{
				if (flag == 0)//劣弧 
				{
					sweepAngle = sweepAngle;
				}
				else if (flag == 2)//半弧
				{
					sweepAngle = 180;
				}
				else if (flag == 1)//优弧
				{
					sweepAngle = 360 - sweepAngle;
				}
				else
				{
					sweepAngle = 360;
				}
			}
			else if (sweepAngle > 180)
			{
				if (flag == 0)//劣弧 
				{
					sweepAngle = 360 - sweepAngle;
				}
				else if (flag == 2)//半弧
				{
					sweepAngle = 180;
				}
				else if (flag == 1)//优弧
				{
					sweepAngle = sweepAngle;
				}
				else
				{
					sweepAngle = 360;
				}
			}

		}
		else
		{
			sweepAngle = 360 + angle_th_end - angle_th;
			if (sweepAngle < 180)
			{
				if (flag == 0)//劣弧 
				{
					sweepAngle = sweepAngle;
				}
				else if (flag == 2)//半弧
				{
					sweepAngle = 180;
				}
				else if (flag == 1)//优弧
				{
					sweepAngle = 360 - sweepAngle;
				}
				else
				{
					sweepAngle = 360;
				}
			}
			else if (sweepAngle > 180)
			{
				if (flag == 0)//劣弧 
				{
					sweepAngle = 360 - sweepAngle;
				}
				else if (flag == 2)//半弧
				{
					sweepAngle = 180;
				}
				else if (flag == 1)//优弧
				{
					sweepAngle = sweepAngle;
				}
				else
				{
					sweepAngle = 360;
				}
			}
		}
	}
	int flagclockwise = 0;//////////////////////////////////////////////////////////0表示顺时针，1表示逆时针
						  //坐标系 用户角度
	float record_x_mid = point_pos[m - 1] - arc_fit.CirCen.x;
	float record_y_mid = arc_fit.CirCen.y - point_pos[m];


	if (point_pos[0] == point_pos[2 * m - 2])//垂直
	{
		if (point_pos[1] < point_pos[2 * m - 1])
		{
			if (point_pos[m - 1] > point_pos[0])
			{
				flagclockwise = 1;
			}
		}
		else if (point_pos[1] > point_pos[2 * m - 1])
		{
			if (point_pos[m - 1] < point_pos[0])
			{
				flagclockwise = 1;
			}
		}

	}
	else if (point_pos[1] == point_pos[2 * m - 1])//水平
	{
		if (point_pos[0] < point_pos[2 * m - 2])
		{
			if (point_pos[m] < point_pos[1])
			{
				flagclockwise = 1;
			}
		}
		else if (point_pos[0] > point_pos[2 * m - 2])
		{
			if (point_pos[m] > point_pos[1])
			{
				flagclockwise = 1;
			}
		}
	}
	else if (point_pos[0] != point_pos[2 * m - 2] && point_pos[1] != point_pos[2 * m - 1])
	{

		if ((point_pos[0] > point_pos[2 * m - 2] &&
			0 > (record_y_end + (record_x_mid - record_x_end)*(record_y - record_y_end) / (record_x - record_x_end) - record_y_mid))
			|| (point_pos[0] < point_pos[2 * m - 2] &&
				0 < (record_y_end + (record_x_mid - record_x_end)*(record_y - record_y_end) / (record_x - record_x_end) - record_y_mid))) {
			flagclockwise = 1;
		}

	}
	else
	{
		status = -1;
		goto end;
	}
	if (flagclockwise == 1)
	{
		arc_fit.sweepAngle = -sweepAngle;
	}
	else
	{
		arc_fit.sweepAngle = sweepAngle;
	}
	//arc_fit.sweepAngle = sweepAngle;



end:
	if (error_weight != nullptr)
		free(error_weight);
	if (A != nullptr)
		free(A);
	if (B != nullptr)
		free(B);
	if (C != nullptr)
		free(C);
	if (recordB != nullptr)
		free(recordB);
	if (recordC != nullptr)
		free(recordC);
	if (error != nullptr)
		free(error);
	return status;
}


//////////////////////////////////////////////////
//VisFitting_Ellipse功能说明：输入边缘点坐标和个数，拟合椭圆
//Input
//	Ipp32f *xy  输入边缘点坐标，x0,y0,x1,y1,x2,y2......
//	int m  边缘点个数

//Output
//	EllipseStruct &resu  输出椭圆结构参数

//Return
//	true - 正常
//	false - 失败
//Author:申健成/詹铭泽/20170412
//////////////////////////////////////////////////
bool CVisAlg2DBase::VisFitting_Ellipse(Ipp32f *xy, int m, EllipseStruct &resu)
{
	if (m < 5)
		return false;
	bool status = true;
	Ipp32f lastsum = 0, nowsum = 0;
	Ipp32f *error_weight;
	error_weight = (Ipp32f *)malloc(5 * m * sizeof(Ipp32f));
	Ipp32f *A;
	A = (Ipp32f *)malloc(5 * sizeof(Ipp32f));
	A[0] = A[1] = A[2] = A[3] = A[4] = 0;
	Ipp32f *B, *C, *recordB, *recordC;
	B = (Ipp32f *)malloc(25 * m * sizeof(Ipp32f));
	C = (Ipp32f *)malloc(5 * m * sizeof(Ipp32f));
	recordB = (Ipp32f *)malloc(25 * m * sizeof(Ipp32f));
	recordC = (Ipp32f *)malloc(5 * m * sizeof(Ipp32f));
	Ipp32f *error;
	error = (Ipp32f *)malloc(5 * m * sizeof(Ipp32f));
	for (int i = 0; i < 5 * m; i++)
	{
		error_weight[i] = 1.0;
		error[i] = 0;
	}
	Ipp32f tempx2 = 0;
	Ipp32f tempx3 = 0;
	Ipp32f tempy2 = 0;
	Ipp32f tempy3 = 0;
	Ipp32f tempy4 = 0;
	int temp_time = 1;
	do {
		lastsum = nowsum;
		nowsum = 0;

		for (int i = 0, j = 0, k = 0; i < 25 * m; i = i + 25, j = j + 2, k = k + 5)
		{
			tempx2 = xy[j] * xy[j];
			tempx3 = xy[j] * xy[j] * xy[j];
			tempy2 = xy[j + 1] * xy[j + 1];
			tempy3 = xy[j + 1] * xy[j + 1] * xy[j + 1];
			tempy4 = xy[j + 1] * xy[j + 1] * xy[j + 1] * xy[j + 1];

			recordB[i] = B[i] = tempx2 * tempy2*error_weight[k];
			recordB[i + 1] = B[i + 1] = xy[j] * tempy3*error_weight[k];
			recordB[i + 2] = B[i + 2] = tempx2 * xy[j + 1] * error_weight[k];
			recordB[i + 3] = B[i + 3] = xy[j] * tempy2*error_weight[k];
			recordB[i + 4] = B[i + 4] = xy[j] * xy[j + 1] * error_weight[k];

			recordB[i + 5] = B[i + 5] = xy[j] * tempy3*error_weight[k + 1];
			recordB[i + 6] = B[i + 6] = tempy4*error_weight[k + 1];
			recordB[i + 7] = B[i + 7] = xy[j] * tempy2*error_weight[k + 1];
			recordB[i + 8] = B[i + 8] = tempy3*error_weight[k + 1];
			recordB[i + 9] = B[i + 9] = tempy2*error_weight[k + 1];

			recordB[i + 10] = B[i + 10] = tempx2 * xy[j + 1] * error_weight[k + 2];
			recordB[i + 11] = B[i + 11] = xy[j] * tempy2*error_weight[k + 2];
			recordB[i + 12] = B[i + 12] = tempx2*error_weight[k + 2];
			recordB[i + 13] = B[i + 13] = xy[j] * xy[j + 1] * error_weight[k + 2];
			recordB[i + 14] = B[i + 14] = xy[j] * error_weight[k + 2];

			recordB[i + 15] = B[i + 15] = xy[j] * tempy2*error_weight[k + 3];
			recordB[i + 16] = B[i + 16] = tempy3*error_weight[k + 3];
			recordB[i + 17] = B[i + 17] = xy[j] * xy[j + 1] * error_weight[k + 3];
			recordB[i + 18] = B[i + 18] = tempy2*error_weight[k + 3];
			recordB[i + 19] = B[i + 19] = xy[j + 1] * error_weight[k + 3];

			recordB[i + 20] = B[i + 20] = xy[j] * xy[j + 1] * error_weight[k + 4];
			recordB[i + 21] = B[i + 21] = tempy2*error_weight[k + 4];
			recordB[i + 22] = B[i + 22] = xy[j] * error_weight[k + 4];
			recordB[i + 23] = B[i + 23] = xy[j + 1] * error_weight[k + 4];
			recordB[i + 24] = B[i + 24] = m * error_weight[k + 4];

			recordC[k] = C[k] = -tempx3 * xy[j + 1] * error_weight[k];
			recordC[k + 1] = C[k + 1] = -tempx2 * tempy2*error_weight[k + 1];
			recordC[k + 2] = C[k + 2] = -tempx3*error_weight[k + 2];
			recordC[k + 3] = C[k + 3] = -tempx2*xy[j + 1] * error_weight[k + 3];
			recordC[k + 4] = C[k + 4] = -tempx2*error_weight[k + 4];
		}

		int info = LAPACKE_sgels(LAPACK_ROW_MAJOR, 'N', 5 * m, 5, 1, B, 5, C, 1);
		if (info != 0)
		{
			status = false;
			goto end;
		}
		A[0] = C[0];
		A[1] = C[1];
		A[2] = C[2];
		A[3] = C[3];
		A[4] = C[4];

		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 5 * m, 1, 5, 1, recordB, 5, A, 1, 1, recordC, 1);

		for (int i = 0; i < (5 * m); i++)
		{
			error[i] = recordC[i];//*recordC[i];
			nowsum = error[i] + nowsum;
		}
		for (int i = 0; i < 5 * m; i++)
		{
			error_weight[i] = 1.0 / (1 + exp((error[i]) / (fitting_sigma*fitting_sigma * 100000)));
		}
		temp_time++;

	} while (fabs(nowsum - lastsum) > 0.02 && temp_time < 0);
	Ipp32f a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
	a = 1;
	b = C[0];
	c = C[1];
	d = C[2];
	e = C[3];
	f = C[4];

	/*cout << C[0] << endl;
	cout << C[1] << endl;
	cout << C[2] << endl;
	cout << C[3] << endl;
	cout << C[4] << endl;*/
	Ipp32f x0 = (2 * C[1] * C[2] - C[0] * C[3]) / (C[0] * C[0] - 4 * C[1]);
	Ipp32f y0 = (2 * C[3] - C[0] * C[3]) / (C[0] * C[0] - 4 * C[1]);
	Ipp32f fenzi = 2 * (C[0] * C[2] * C[3] - C[1] * C[2] * C[2] - C[3] * C[3] + 4 * C[4] * C[1] - C[0] * C[0] * C[4]);

	Ipp32f fenmu = (C[0] * C[0] - 4 * C[1]) * (C[1] - sqrt(C[0] * C[0] + (1 - C[1]) * (1 - C[1])) + 1);

	Ipp32f femmu2 = (C[0] * C[0] - 4 * C[1]) * (C[1] + sqrt(C[0] * C[0] + (1 - C[1]) * (1 - C[1])) + 1);

	//long temp_axis1= sqrt(fabs(fenzi / fenmu));

	//short temp_axis2= sqrt(fabs(fenzi / femmu2));

	//Ipp32f theta = atan(sqrt((x0 *x0 - y0*y0*C[1]) / (x0 *x0 *C[1] - y0*y0)) + 0.0001) * 180 / PI;

	Ipp32f center_x = 0;
	Ipp32f center_y = 0;
	if (0 != (C[0] * C[0] - 4 * C[1]))
	{
		center_x = (2 * C[1] * C[2] - C[0] * C[3]) / (C[0] * C[0] - 4 * C[1]);
		center_y = (2 * C[3] - C[0] * C[2]) / (C[0] * C[0] - 4 * C[1]);
	}

	Ipp32f axis_a = sqrt(fabs(fenzi / fenmu));
	Ipp32f axis_b = sqrt(fabs(fenzi / femmu2));
	Ipp32f rotateTheta = 0;
	float temp_d_x_max = 0; float temp_d_y_max = 0;
	float temp_d_x_min = 0; float temp_d_y_min = 0;
	if (axis_a > axis_b)
	{
		double temp_d_max = 0; double temp_d_min = sqrt((xy[0] - center_x)*(xy[0] - center_x) + (xy[1] - center_y)*(xy[1] - center_y));
		for (int i = 0; i < 2 * m; i = i + 2)
		{
			double d = sqrt((xy[i] - center_x)*(xy[i] - center_x) + (xy[i + 1] - center_y)*(xy[i + 1] - center_y));
			if (d > temp_d_max)
			{
				temp_d_max = d;
				temp_d_x_max = xy[i];
				temp_d_y_max = xy[i + 1];
			}
			if (d < temp_d_min)
			{
				temp_d_min = d;
				temp_d_x_min = xy[i];
				temp_d_y_min = xy[i + 1];
			}
			/*if (d >= axis_a&&d <= (axis_a + 2))
			{
			if (center_x <= xy[i])
			{
			rotateTheta = atan2(xy[i + 1] - center_y, xy[i] - center_x) * 180 / PI;
			}
			else
			{
			rotateTheta = atan2(center_y - xy[i + 1], center_x - xy[i]) * 180 / PI;
			}
			break;
			}*/
		}
		if (rotateTheta == 0)
		{
			rotateTheta = atan2(temp_d_y_max - center_y, temp_d_x_max - center_x) * 180 / PI;
		}
	}
	else {
		double temp_d_max = 0; double temp_d_min = sqrt((xy[0] - center_x)*(xy[0] - center_x) + (xy[1] - center_y)*(xy[1] - center_y));

		for (int i = 0; i < 2 * m; i = i + 2)
		{
			double d = sqrt((xy[i] - center_x)*(xy[i] - center_x) + (xy[i + 1] - center_y)*(xy[i + 1] - center_y));
			if (d > temp_d_max)
			{
				temp_d_max = d;
				temp_d_x_max = xy[i];
				temp_d_y_max = xy[i + 1];
			}
			if (d < temp_d_min)
			{
				temp_d_min = d;
				temp_d_x_min = xy[i];
				temp_d_y_min = xy[i + 1];
			}
			/*if (d >= axis_b&&d <= (axis_b + 2))
			{
			if (center_x <= xy[i])
			{
			rotateTheta = atan2(xy[i + 1] - center_y, xy[i] - center_x) * 180 / PI;
			}
			else
			{
			rotateTheta = atan2(center_y - xy[i + 1], center_x - xy[i]) * 180 / PI;
			}
			break;
			}*/
		}
		if (rotateTheta == 0)
		{
			rotateTheta = atan2(temp_d_y_max - center_y, temp_d_x_max - center_x) * 180 / PI;
		}
	}
	axis_a = sqrt((temp_d_x_max - center_x)*(temp_d_x_max - center_x) + (temp_d_y_max - center_y)*(temp_d_y_max - center_y));
	axis_b = sqrt((temp_d_x_min - center_x)*(temp_d_x_min - center_x) + (temp_d_y_min - center_y)*(temp_d_y_min - center_y));
	resu.a = a;
	resu.b = b;
	resu.c = c;
	resu.d = d;
	resu.e = e;
	resu.f = f;
	resu.center_x = center_x;
	resu.center_y = center_y;
	resu.axis_a = axis_a;
	resu.axis_b = axis_b;
	resu.rotateTheta = rotateTheta;

end:
	if (error_weight != nullptr) free(error_weight);
	if (A != nullptr) free(A);
	if (B != nullptr) free(B);
	if (C != nullptr) free(C);
	if (recordB != nullptr) free(recordB);
	if (recordC != nullptr) free(recordC);
	if (error != nullptr) free(error);
	return status;
}

//////////////////////////////////////////////////
//VisOmmET_Fitting_Ellipse_Str功能说明：输入边缘点坐标和个数，拟合椭圆
//Input
//	Ipp32f *xy  输入边缘点坐标，x0,y0,x1,y1,x2,y2......
//	int m  边缘点个数

//Output
//	EllipseStruct &resu  输出椭圆结构参数

//Return
//	true - 正常
//	false - 失败
//Author:申健成/詹铭泽/20170412
//////////////////////////////////////////////////
bool  CVisAlg2DBase::VisFitting_Ellipse_Str(Ipp32f *xy, int m, EllipseStruct &resu)
{
	if (m < 5)
		return false;
	float *A, *X;
	A = (float*)malloc(5 * m * sizeof(float));
	X = (float*)malloc(m * sizeof(float));
	for (int i = 0, j = 0, k = 0; i < 5 * m; i = i + 5, j = j + 2, k++)
	{
		A[i] = xy[j + 1] * xy[j + 1];
		A[i + 1] = xy[j] * xy[j + 1];
		A[i + 2] = xy[j];
		A[i + 3] = xy[j + 1];
		A[i + 4] = 1;
		X[k] = -xy[j] * xy[j];
	}
	float *AT;
	AT = (float*)malloc(5 * m * sizeof(float));
	MKL_Somatcopy('r', 't', m, 5, 1, A, 5, AT, m);
	float *ATA;
	ATA = (float*)malloc(5 * 5 * sizeof(float));
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 5, 5, m, 1, AT, m, A, 5, 0, ATA, 5);
	double temp_ata[25];
	for (int i = 0; i < 25; i++)
	{
		temp_ata[i] = (double)ATA[i];
	}
	inverse(temp_ata, 5);
	for (int i = 0; i < 25; i++)
	{
		ATA[i] = temp_ata[i];
	}

	float *tempmat;
	tempmat = (float*)malloc(5 * m * sizeof(float));
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 5, m, 5, 1, ATA, 5, AT, m, 0, tempmat, m);
	float *P;
	P = (float*)malloc(5 * sizeof(float));
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 5, 1, m, 1, tempmat, m, X, 1, 0, P, 1);

	double P2[4] = { -2,-P[1],-P[1],-2 * P[0] };
	inverse(P2, 2);
	float ceter_xy[2] = { 0,0 };
	float P34[2] = { P[2],P[3] };
	float tempp2[4] = { P2[0],P2[1],P2[2],P2[3] };
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 1, 2, 1, tempp2, 2, P34, 1, 0, ceter_xy, 1);

	double u = 1.0 / (ceter_xy[0] * ceter_xy[0] + P[0] * ceter_xy[1] * ceter_xy[1] + P[1] * ceter_xy[0] * ceter_xy[1] - P[4]);
	double v = u*P[0];
	double w = u*P[1];
	double x = u*P[2];
	double y = u*P[3];
	double z = u*P[4] + 1;

	double tanth = (u - v) / w - sqrt(((u - v) / w)*((u - v) / w) + 1);
	double a_axis = cos(atan(tanth)) / sqrt((u - tanth*tanth*v) / (1 - tanth*tanth*tanth*tanth));
	double b_axis = cos(atan(tanth)) / sqrt((v - tanth*tanth*u) / (1 - tanth*tanth*tanth*tanth));
	double tanth1 = (u - v) / w + sqrt(((u - v) / w)*((u - v) / w) + 1);
	double actan = atan(tanth1) * 180 / PI;
	resu.center_x = ceter_xy[0];
	resu.center_y = ceter_xy[1];



	float temp_d_x_max = 0; float temp_d_y_max = 0;
	float temp_d_x_min = 0; float temp_d_y_min = 0;

	if (a_axis > b_axis)
	{
		double temp_d_max = 0; double temp_d_min = sqrt((xy[0] - ceter_xy[0])*(xy[0] - ceter_xy[0]) + (xy[1] - ceter_xy[1])*(xy[1] - ceter_xy[1]));
		for (int i = 0; i < 2 * m; i = i + 2)
		{
			double d = sqrt((xy[i] - ceter_xy[0])*(xy[i] - ceter_xy[0]) + (xy[i + 1] - ceter_xy[1])*(xy[i + 1] - ceter_xy[1]));
			if (d > temp_d_max)
			{
				temp_d_max = d;
				temp_d_x_max = xy[i];
				temp_d_y_max = xy[i + 1];
			}
			if (d < temp_d_min)
			{
				temp_d_min = d;
				temp_d_x_min = xy[i];
				temp_d_y_min = xy[i + 1];
			}
			/*if (d >= a_axis&&d <= (a_axis + 2))
			{
			if (ceter_xy[0] <= xy[i])
			{
			resu.rotateTheta = atan2(xy[i + 1] - ceter_xy[1], xy[i] - ceter_xy[0]) * 180 / PI;
			}
			else
			{
			resu.rotateTheta = atan2(ceter_xy[1] - xy[i + 1], ceter_xy[0] - xy[i]) * 180 / PI;
			}
			break;
			}*/
		}
		//if (resu.rotateTheta == 0)
		{
			resu.rotateTheta = atan2(temp_d_y_max - ceter_xy[1], temp_d_x_max - ceter_xy[0]) * 180 / PI;
			//resu.rotateTheta = atan2(temp_d_y - ceter_xy[1], temp_d_x - ceter_xy[0]) * 180 / PI;
		}

	}
	else
	{
		double temp_d_max = 0; double temp_d_min = sqrt((xy[0] - ceter_xy[0])*(xy[0] - ceter_xy[0]) + (xy[1] - ceter_xy[1])*(xy[1] - ceter_xy[1]));
		for (int i = 0; i < 2 * m; i = i + 2)
		{
			double d = sqrt((xy[i] - ceter_xy[0])*(xy[i] - ceter_xy[0]) + (xy[i + 1] - ceter_xy[1])*(xy[i + 1] - ceter_xy[1]));
			if (d > temp_d_max)
			{
				temp_d_max = d;
				temp_d_x_max = xy[i];
				temp_d_y_max = xy[i + 1];
			}
			if (d < temp_d_min)
			{
				temp_d_min = d;
				temp_d_x_min = xy[i];
				temp_d_y_min = xy[i + 1];
			}
			/*if (d >= b_axis&&d <= (b_axis + 2))
			{
			if (ceter_xy[0] <= xy[i])
			{
			resu.rotateTheta = atan2(xy[i + 1] - ceter_xy[1], xy[i] - ceter_xy[0]) * 180 / PI;
			}
			else
			{
			resu.rotateTheta = atan2(ceter_xy[1] - xy[i + 1], ceter_xy[0] - xy[i]) * 180 / PI;
			}
			break;
			}*/
		}
		//if (resu.rotateTheta == 0)
		{
			resu.rotateTheta = atan2(temp_d_y_max - ceter_xy[1], temp_d_x_max - ceter_xy[0]) * 180 / PI;
			//resu.rotateTheta = atan2(temp_d_y - ceter_xy[1], temp_d_x - ceter_xy[0]) * 180 / PI;
		}
	}
	resu.axis_a = sqrt((temp_d_x_max - ceter_xy[0])*(temp_d_x_max - ceter_xy[0]) + (temp_d_y_max - ceter_xy[1])*(temp_d_y_max - ceter_xy[1]));
	resu.axis_b = sqrt((temp_d_x_min - ceter_xy[0])*(temp_d_x_min - ceter_xy[0]) + (temp_d_y_min - ceter_xy[1])*(temp_d_y_min - ceter_xy[1]));

	//resu.rotateTheta = 90 - actan;
	//resu.rotateTheta = (atan(tanth)+0.0001)*180/PI;
	//atan(sqrt((x0 *x0 - y0*y0*C[1]) / (x0 *x0 *C[1] - y0*y0)) + 0.0001) * 180 / PI;
	if (A != nullptr) free(A);
	if (X != nullptr) free(X);
	if (P != nullptr) free(P);
	if (AT != nullptr) free(AT);
	if (ATA != nullptr) free(ATA);
	if (tempmat != nullptr) free(tempmat);

	return true;
}
IMG_VVOID CVisAlg2DBase::inverse(IMG_LREAL* A, IMG_INT N)
{
	int *IPIV = new int[N + 1];
	int LWORK = N*N;
	double *WORK = new double[LWORK];
	int INFO;

	dgetrf_(&N, &N, A, &N, IPIV, &INFO);
	dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);

	delete IPIV;
	delete WORK;
}

/**********************************************/
// GNERAL_LINE_FITTING, 通用直线拟合
// Input:
//     vector<IMG_RCOORD>PointCor, 点集合
//     IMG_REAL &Slope, 初始斜率输入，FLAG = 0 的时候需输入
//     IMG_REAL &B, 初始截距输入，FLAG = 0 的时候需输入
//     IMG_REAL Sigma, 权重sigma
//     IMG_INT FLAG，FLAG=0的时候需要输入初始直线斜率以及截距，FLAG=1时不需要输入
// Output:
//     IMG_REAL &Slope, 拟合后的斜率
//     IMG_REAL &B, 拟合后的截距
// Return:
//     -1 - 异常
//     0  - 正常
// Author: Tan Ling/20170412
/**********************************************/
IMG_INT CVisAlg2DBase::GNERAL_LINE_FITTING(vector<IMG_RCOORD>PointCor, IMG_REAL &Slope, IMG_REAL &B, IMG_REAL Sigma, IMG_INT FLAG)
{
	if (PointCor.size() == 0)
		return -1;
	IMG_REAL *MatA, *MatA_Copy, *MatB, *MatC, *MatC_Copy, *Weight, *ErrorRecord;
	IMG_INT ErrorMaxTimes = 20, ErrorTimes, Length = PointCor.size();//length means efficient point numbers subject to sigma

	MatA = (IMG_REAL*)malloc(sizeof(IMG_REAL) * 2 * PointCor.size()); //column: 2, rows: number of points
	MatB = (IMG_REAL*)malloc(sizeof(IMG_REAL) * 2);
	MatC = (IMG_REAL*)malloc(sizeof(IMG_REAL) * PointCor.size());
	MatA_Copy = (IMG_REAL*)malloc(sizeof(IMG_REAL) * 2 * PointCor.size()); //column: 2, rows: number of points
	MatC_Copy = (IMG_REAL*)malloc(sizeof(IMG_REAL) * PointCor.size());
	ErrorRecord = (IMG_REAL*)malloc(sizeof(IMG_REAL) * ErrorMaxTimes * 3);
	Weight = (IMG_REAL*)malloc(sizeof(IMG_REAL) * PointCor.size());
	for (IMG_INT j = 0; j < (IMG_INT)PointCor.size(); j++)
	{
		MatA[j * 2 + 0] = PointCor[j].x;
		MatA[j * 2 + 1] = 1;

		MatA_Copy[j * 2 + 0] = MatA[j * 2 + 0];
		MatA_Copy[j * 2 + 1] = MatA[j * 2 + 1];

		MatC[j] = PointCor[j].y;
		MatC_Copy[j] = MatC[j];
		Weight[j] = 1;
	}

	switch (FLAG)
	{
	case 0:
		MatB[0] = Slope;
		MatB[1] = B;
		break;
	case 1:
		//step: initial fitting
		LAPACKE_sgels(LAPACK_ROW_MAJOR, 'N', PointCor.size(), 2, 1, MatA_Copy, 2, MatC_Copy, 1);//fitting slope and b
		MatB[0] = MatC_Copy[0];
		MatB[1] = MatC_Copy[1];
	default:
		break;
	}

	////step: choose points subject to sigma
	for (IMG_INT j = 0; j < Length; j++)
	{
		MatA_Copy[j * 2 + 0] = MatA[j * 2 + 0];
		MatA_Copy[j * 2 + 1] = MatA[j * 2 + 1];
		MatC_Copy[j] = MatC[j];
	}



	//step: fitting 
	ErrorTimes = 0;
	while (1)
	{
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Length, 1, 2, 1, MatA, 2, MatB, 1, -1, MatC_Copy, 1); //calculate error
		ErrorRecord[ErrorTimes * 3 + 1] = MatB[0];
		ErrorRecord[ErrorTimes * 3 + 2] = MatB[1];
		for (IMG_INT j = 0; j < Length; j++)
		{
			ErrorRecord[ErrorTimes * 3] += MatC_Copy[j];
		}
		ErrorRecord[ErrorTimes * 3] /= (IMG_REAL)Length;
		if (ErrorTimes >= 1 && abs(ErrorRecord[ErrorTimes * 3] - ErrorRecord[(ErrorTimes - 1) * 3]) <= 0.001) //error sum covergence to 0.001
		{
			Slope = MatB[0];
			B = MatB[1];
			break;
		}
		if (ErrorTimes == ErrorMaxTimes - 1) // runing time > max time, choose lowest error sum to return
		{
			IMG_INT FLAG = 0;
			for (IMG_INT i = 0; i < ErrorMaxTimes - 1; i++)
			{
				if (ErrorRecord[3 * i] < ErrorRecord[3 * (i + 1)])
					ErrorRecord[3 * (i + 1)] = ErrorRecord[3 * i];
				else
					FLAG = i + 1;
			}

			Slope = ErrorRecord[FLAG * 3 + 1];
			B = ErrorRecord[FLAG * 3 + 2];
			break;
		}
		ErrorTimes++;
		for (IMG_INT j = 0; j < Length; j++)
		{
			MatA_Copy[j * 2 + 0] = MatA[j * 2 + 0];
			MatA_Copy[j * 2 + 1] = MatA[j * 2 + 1];

			Weight[j] = (IMG_REAL)1.0 / (1.f + exp(MatC_Copy[j] * MatC_Copy[j] / (Sigma * Sigma)));

			MatA_Copy[j * 2 + 0] = Weight[j] * MatA_Copy[j * 2 + 0];
			MatA_Copy[j * 2 + 1] = Weight[j] * MatA_Copy[j * 2 + 1];


			MatC_Copy[j] = MatC[j] * Weight[j]; //reset MatCCopy after calculating error
		}
		LAPACKE_sgels(LAPACK_ROW_MAJOR, 'N', Length, 2, 1, MatA_Copy, 2, MatC_Copy, 1);//fitting slope and b
		MatB[0] = MatC_Copy[0];
		MatB[1] = MatC_Copy[1];
		for (IMG_INT j = 0; j < Length; j++) //reset MatCCopy after fitting line
		{
			MatC_Copy[j] = MatC[j];
		}
	}
	free(MatC_Copy);
	free(Weight);
	free(MatA_Copy);
	free(MatA);
	free(MatB);
	free(MatC);
	free(ErrorRecord);
	return 0;
}



/**********************************************  Histogram related  *******************************************/

/**********************************************/
// VisHistogram, 统计直方图
// Input:
//		pSrc ， 输入图像指针
//		nW,nH	图像宽高
// Output:
//		pDst ， 输出图像指针
// Return:
//		ippstatus
// Liu Ping, 2017.05.09
/**********************************************/
int CVisAlg2DBase::VisHistogram(IMG_UBBUF src, Ipp32u *pHist, int nBins)
{
	__try
	{
		//1、统计每个灰度值出现的次数
		IppStatus status = ippStsNoErr;

		IppDataType dataType = ipp8u;
		IppiSize srcSize = { src.size.width,src.size.height };
		int srcStep = src.size.width * sizeof(unsigned char);

		//int nBins = 255;
		int nLevels[] = { nBins + 1 };
		int uniform = 1;
		int specSize = 0;
		int bufferSize = 0;

		status = ippiHistogramGetBufferSize(dataType, srcSize, nLevels, 1, uniform, &specSize, &bufferSize);
		if (status != ippStsNoErr) return -1;

		Ipp32f lowerLevel[] = { 0 };
		Ipp32f upperLevel[] = { 255 };
		IppiHistogramSpec *Spec = NULL;
		Spec = (IppiHistogramSpec*)ippsMalloc_8u(specSize);
		status = ippiHistogramUniformInit(dataType, lowerLevel, upperLevel, nLevels, 1, Spec);
		if (status != ippStsNoErr)
		{
			ippsFree(Spec);
			return -1;
		}

		Ipp32f *pLevels, *ppLevels[1];
		pLevels = (Ipp32f*)malloc((nBins + 1) * sizeof(Ipp32f));

		ppLevels[0] = pLevels;
		status = ippiHistogramGetLevels(Spec, ppLevels);

		//Ipp32u *pHist = (Ipp32u*)malloc((nBins) * sizeof(Ipp32u));
		Ipp8u *buffer = NULL;
		buffer = (Ipp8u*)malloc(bufferSize * sizeof(Ipp8u));

		status = ippiHistogram_8u_C1R(src.ptr, srcStep, srcSize, pHist, Spec, buffer);

		//free
		free(pLevels);
		free(buffer);
		ippsFree(Spec);

		if (status != ippStsNoErr)
		{
			return -1;
		}

		return 0;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
		return _CODE_THROW;
	}
}



/**********************************************  Segment *******************************************/

/**********************************************/
// threshold_manual, 手动阈值分割
// Input:
//		pSrc ， 输入图像指针
//		nW,nH	图像宽高
//		nT	分割阈值
// Output:
//		pDst ， 输出图像指针
// Return:
//		ippstatus
// Liu Ping, 2017.05.06
/**********************************************/
int CVisAlg2DBase::VisSegmentManual(const unsigned char * pSrc, const unsigned int nW, const unsigned int nH, unsigned char * pDst, const unsigned char nT)
{
	__try
	{
		IppiSize imgSize;
		imgSize.width = nW;
		imgSize.height = nH;
		int status = ippStsNoErr;
		int width = imgSize.width;
		int height = imgSize.height;

		//	do threshold	//
		int srcStep = width * sizeof(unsigned char);
		int dstStep = srcStep;
		status = ippiThreshold_LTValGTVal_8u_C1R(pSrc, srcStep, pDst, dstStep, imgSize, nT, 0, nT, 255);

		return status;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
		return _CODE_THROW;
	}
}

/**********************************************/
// threshold_OTSU, 自动阈值分割
// Input:
//		pSrc ， 输入图像指针
//		nW,nH	图像宽高
// Output:
//		pDst ， 输出图像指针
// Return:
//		ippstatus
// Liu Ping, 2017.05.07
/**********************************************/
int CVisAlg2DBase::VisSegmentOtsu(const unsigned char * pSrc, const unsigned int nW, const unsigned int nH, unsigned char * pDst)
{
	__try
	{
		IppiSize imgSize;
		imgSize.width = nW;
		imgSize.height = nH;
		int status = ippStsNoErr;
		int width = imgSize.width;
		int height = imgSize.height;

		//	calc threshold	//
		int srcStep = width * sizeof(unsigned char);
		unsigned char nThreshold = 0;
		status = ippiComputeThreshold_Otsu_8u_C1R(pSrc, srcStep, imgSize, &nThreshold);
		if (status != 0)
			throw 1;

		//	do segment	//
		int dstStep = srcStep;
		status = ippiThreshold_LTValGTVal_8u_C1R(pSrc, srcStep, pDst, dstStep, imgSize, nThreshold, 0, nThreshold, 255);

		return status;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
		return _CODE_THROW;
	}
}

/**********************************************/
// threshold_OTSU, 自动阈值计算，不分割
// Input:
//		pSrc ， 输入图像指针
//		nW,nH	图像宽高
// Output:
//		pDst ， 输出图像指针
// Return:
//		ippstatus
// Liu Ping, 2017.05.07
/**********************************************/
int CVisAlg2DBase::VisCalcOtsu(const unsigned char * pSrc, const unsigned int nW, const unsigned int nH, unsigned char & nThres)
{
	__try
	{
		IppiSize imgSize;
		imgSize.width = nW;
		imgSize.height = nH;
		int status = ippStsNoErr;
		int width = imgSize.width;
		int height = imgSize.height;

		//	calc threshold	//
		int srcStep = width * sizeof(unsigned char);
		status = ippiComputeThreshold_Otsu_8u_C1R(pSrc, srcStep, imgSize, &nThres);
		if (status != 0)
			throw 1;

		return status;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
		return _CODE_THROW;
	}
}

/**********************************************/
// threshold_Dynamic, 动态阈值分割
// Input:
//		pSrc ， 输入图像指针
//		nW,nH	图像宽高
// Output:
//		pDst ， 输出图像指针
// Return:
//		ippstatus
// Liu Ping, 2017.05.08
/**********************************************/
int CVisAlg2DBase::VisSegmentDynamic(const unsigned char * pSrc, const unsigned int nW, const unsigned int nH, unsigned char * pDst, const unsigned int avgWinWidth)
{
	__try
	{
		IppiSize imgSize;
		imgSize.width = nW;
		imgSize.height = nH;
		int status = ippStsNoErr;
		int width = imgSize.width;
		int height = imgSize.height;
		int nS = width*height;

		// calc average image
		unsigned char * pThre = new unsigned char[nS];

		VisFilterMedian(pSrc, height, width, pThre, avgWinWidth);

		// do segment
		for (int j = 0; j < nS; j++)
		{
			if ((pSrc[j]) <= pThre[j])
			{
				pDst[j] = 0;
			}
			else
			{
				pDst[j] = 255;
			}
		}

		delete[] pThre;
		return status;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
		return _CODE_THROW;
	}
}

/**********************************************/
// threshold_CountDots, 数点阈值分割
// Input:
//		pSrc ， 输入图像指针
//		nW,nH	图像宽高
// Output:
//		pDst ， 输出图像指针
// Return:
//		ippstatus
// Liu Ping, 2017.05.09
/**********************************************/
int CVisAlg2DBase::VisSegmentCntDots(unsigned char * pSrc, const unsigned int nW, const unsigned int nH, unsigned char * pDst, const double fWhiteRatio)
{
	__try
	{
		int nBins = 255;

		int status = ippStsNoErr;
		unsigned char nThreshold = 0;

		//	calc threshold	//
		Ipp32u *pHist = new Ipp32u[nBins];
		IMG_UBBUF src;
		src.ptr = pSrc;
		src.size.width = nW;
		src.size.height = nH;
		src.linestep = nW;

		//VisHistogram(src, pHist, nBins);
		//IppStatus status = ippStsNoErr;

		IppDataType dataType = ipp8u;
		IppiSize srcSize = { src.size.width,src.size.height };
		int srcStep = src.size.width * sizeof(unsigned char);

		int nLevels[] = { nBins + 1 };
		int uniform = 1;
		int specSize = 0;
		int bufferSize = 0;

		status = ippiHistogramGetBufferSize(dataType, srcSize, nLevels, 1, uniform, &specSize, &bufferSize);
		if (status != ippStsNoErr) return -1;

		Ipp32f lowerLevel[] = { 0 };
		Ipp32f upperLevel[] = { 255 };
		IppiHistogramSpec *Spec = NULL;
		Spec = (IppiHistogramSpec*)ippsMalloc_8u(specSize);
		status = ippiHistogramUniformInit(dataType, lowerLevel, upperLevel, nLevels, 1, Spec);
		if (status != ippStsNoErr)
		{
			ippsFree(Spec);
			return -1;
		}

		Ipp32f *pLevels, *ppLevels[1];
		pLevels = (Ipp32f*)malloc((nBins + 1) * sizeof(Ipp32f));

		ppLevels[0] = pLevels;
		status = ippiHistogramGetLevels(Spec, ppLevels);

		//Ipp32u *pHist = (Ipp32u*)malloc((nBins) * sizeof(Ipp32u));
		Ipp8u *buffer = NULL;
		buffer = (Ipp8u*)malloc(bufferSize * sizeof(Ipp8u));

		status = ippiHistogram_8u_C1R(src.ptr, srcStep, srcSize, pHist, Spec, buffer);

		//free
		free(pLevels);
		free(buffer);
		ippsFree(Spec);

		unsigned int i, nCnt = 0;
		unsigned int nWhite = unsigned int(nW*nH*fWhiteRatio + 0.5);
		for (i = nBins - 1; i > 0; i--)
		{
			nCnt += pHist[i];
			if (nCnt > nWhite)
			{
				nThreshold = i;
				break;
			}
		}
		delete[] pHist;

		//	do segment	//
		//int srcStep = nW * sizeof(unsigned char);
		int dstStep = srcStep;
		IppiSize imgSize;
		imgSize.width = nW;
		imgSize.height = nH;
		status = ippiThreshold_LTValGTVal_8u_C1R(pSrc, srcStep, pDst, dstStep, imgSize, nThreshold, 0, nThreshold, 255);

		return 0;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
		return _CODE_THROW;
	}
}


/**********************************************/
// VisSS_GetBlobThreshold, 功能说明:获取动态阈值对图像进行二值化。
// Input:
// 	Ipp8u *srcRoi, 输入图片
//	IppiSize Roi,感兴趣区域大小
//  int BlobThreshold,将图像分为几段，按行分段（可以分为5段，6段，7段）
//     ...
// Output:
//     ...
// Return:
//     0 - 正常
// Author: Huang Yige/05/02/2017
/**********************************************/
void CVisAlg2DBase::VisGetBlobThreshold(Ipp8u *srcRoi, IppiSize Roi, int BlobThreshold)
{

	Ipp8u* dstRoi;
	dstRoi = (Ipp8u*)malloc(Roi.width * sizeof(Ipp8u));

	int SmallRoiY = (int)((float)Roi.height / (float)BlobThreshold);
	int srcDstStep = Roi.width * sizeof(Ipp8u);

	Ipp8u* SmallRoi1;
	SmallRoi1 = (Ipp8u*)malloc(Roi.width*SmallRoiY * sizeof(Ipp8u));
	IppiSize ThresholdRoiFirst = { Roi.width ,SmallRoiY };

	Ipp8u* SmallRoiLast;
	SmallRoiLast = (Ipp8u*)malloc(Roi.width*(Roi.height - (BlobThreshold - 1) * SmallRoiY) * sizeof(Ipp8u));
	IppiSize ThresholdRoiLast = { Roi.width ,Roi.height - (BlobThreshold - 1) * SmallRoiY };

	Ipp64f* pMinThreshold;
	pMinThreshold = (Ipp64f*)malloc(BlobThreshold * sizeof(Ipp64f));

	Ipp64f* pMinFilter;
	pMinFilter = (Ipp64f*)malloc(Roi.height * sizeof(Ipp64f));

	IppiSize FilterRoi = { Roi.width,1 };

	for (int i = 0; i < (BlobThreshold - 1); i++)
	{
		for (int j = 0; j < Roi.width*SmallRoiY; j++)
		{
			SmallRoi1[j] = srcRoi[i*(Roi.width*SmallRoiY) + j];
		}
		ippiMean_8u_C1R(SmallRoi1, srcDstStep, ThresholdRoiFirst, &pMinThreshold[i]);
	}

	for (int j = 0; j < ThresholdRoiLast.width*ThresholdRoiLast.height; j++)
	{
		SmallRoiLast[j] = srcRoi[(Roi.width*SmallRoiY) * (BlobThreshold - 1) + j];
	}

	ippiMean_8u_C1R(SmallRoiLast, srcDstStep, ThresholdRoiLast, &pMinThreshold[(BlobThreshold - 1)]);

	//最小二乘法求取直线
	float A = 0.0, B = 0.0;
	if (BlobThreshold == 5)
	{
		A = (float)((5.0 * (0.5 * (float)pMinThreshold[0] + 1.5 * (float)pMinThreshold[1] + 2.5 * (float)pMinThreshold[2] + 3.5 * (float)pMinThreshold[3] + 4.5 * (float)pMinThreshold[4]) -
			(0.5 + 1.5 + 2.5 + 3.5 + 4.5)*(float)(pMinThreshold[0] + pMinThreshold[1] + pMinThreshold[2] + pMinThreshold[3] + pMinThreshold[4])) /
			(5.0*(0.5*0.5 + 1.5*1.5 + 2.5*2.5 + 3.5*3.5 + 4.5*4.5) - (0.5 + 1.5 + 2.5 + 3.5 + 4.5)*(0.5 + 1.5 + 2.5 + 3.5 + 4.5)));
		B = (float)((pMinThreshold[0] + pMinThreshold[1] + pMinThreshold[2] + pMinThreshold[3] + pMinThreshold[4]) / 5.0 -
			A * (0.5 + 1.5 + 2.5 + 3.5 + 4.5) / 5.0);
	}
	if (BlobThreshold == 6)
	{
		A = (float)((6.0 * (0.5 * (float)pMinThreshold[0] + 1.5 * (float)pMinThreshold[1] + 2.5 * (float)pMinThreshold[2] + 3.5 * (float)pMinThreshold[3] + 4.5 * (float)pMinThreshold[4] + 5.5 * (float)pMinThreshold[5]) -
			(0.5 + 1.5 + 2.5 + 3.5 + 4.5 + 5.5)*((float)pMinThreshold[0] + (float)pMinThreshold[1] + (float)pMinThreshold[2] + (float)pMinThreshold[3] + (float)pMinThreshold[4] + (float)pMinThreshold[5])) /
			(6.0*(0.5*0.5 + 1.5*1.5 + 2.5*2.5 + 3.5*3.5 + 4.5*4.5 + 5.5*5.5) - (0.5 + 1.5 + 2.5 + 3.5 + 4.5 + 5.5)*(0.5 + 1.5 + 2.5 + 3.5 + 4.5 + 5.5)));
		B = (float)(((float)pMinThreshold[0] + (float)pMinThreshold[1] + (float)pMinThreshold[2] + (float)pMinThreshold[3] + (float)pMinThreshold[4] + (float)pMinThreshold[5]) / 6.0 -
			A * (0.5 + 1.5 + 2.5 + 3.5 + 4.5 + 5.5) / 6.0);
	}
	if (BlobThreshold == 7)
	{
		A = (float)((7.0 * (0.5 * (float)pMinThreshold[0] + 1.5 * (float)pMinThreshold[1] + 2.5 * (float)pMinThreshold[2] + 3.5 * (float)pMinThreshold[3] + 4.5 * (float)pMinThreshold[4] + 5.5 * (float)pMinThreshold[5] + 6.5 * (float)pMinThreshold[6]) -
			(0.5 + 1.5 + 2.5 + 3.5 + 4.5 + 5.5 + 6.5)*((float)pMinThreshold[0] + (float)pMinThreshold[1] + (float)pMinThreshold[2] + (float)pMinThreshold[3] + (float)pMinThreshold[4] + (float)pMinThreshold[5] + (float)pMinThreshold[6])) /
			(7.0*(0.5*0.5 + 1.5*1.5 + 2.5*2.5 + 3.5*3.5 + 4.5*4.5 + 5.5*5.5 + 6.5*6.5) - (0.5 + 1.5 + 2.5 + 3.5 + 4.5 + 5.5 + 6.5)*(0.5 + 1.5 + 2.5 + 3.5 + 4.5 + 5.5 + 6.5)));
		B = (float)(((float)pMinThreshold[0] + (float)pMinThreshold[1] + (float)pMinThreshold[2] + (float)pMinThreshold[3] + (float)pMinThreshold[4] + (float)pMinThreshold[5] + (float)pMinThreshold[6]) / 7.0 -
			A * (0.5 + 1.5 + 2.5 + 3.5 + 4.5 + 5.5 + 6.5) / 7.0);
	}


	float coodrate_x = 0.0;

	for (int i = 0; i < Roi.height; i++)
	{
		coodrate_x = (float)i / (float)SmallRoiY;

		pMinFilter[i] = B + A*coodrate_x;
	}

	for (int i = 0; i < Roi.height; i++)
	{
		for (int j = 0; j < Roi.width; j++)
		{
			dstRoi[j] = srcRoi[i*Roi.width + j];
		}

		if ((int)pMinFilter[i]  < 0)
		{
			ippiThreshold_LTValGTVal_8u_C1IR(dstRoi, srcDstStep, FilterRoi, 0, 255, 0, 0);
		}
		else
		{
			ippiThreshold_LTValGTVal_8u_C1IR(dstRoi, srcDstStep, FilterRoi, (int)pMinFilter[i], 255, (int)pMinFilter[i], 0);

		}

		for (int j = 0; j < Roi.width; j++)
		{
			srcRoi[i*Roi.width + j] = dstRoi[j];
		}
	}

	for (int i = 0; i < Roi.height; i++)
	{
		for (int j = 0; j < Roi.width; j++)
		{
			if (srcRoi[i*Roi.width + j] != 0 && srcRoi[i*Roi.width + j] != 255)
			{
				srcRoi[i*Roi.width + j] = 0;
			}
		}
	}

	free(dstRoi);
	free(SmallRoiLast);
	free(SmallRoi1);
	free(pMinThreshold);
	free(pMinFilter);
}


/**********************************************  Morphology *******************************************/
/**********************************************/
// VisMorphologyErode, 形态学腐蚀, 窗宽(3 或 5)
// Input:
//		src		输入图像
//  Size        3或者5或者7
// Output:
//		dst		输出图像
// Return:
//		ippstatus
// Liu Ping, 2017.05.09
/**********************************************/
int CVisAlg2DBase::VisMorphologyErode(unsigned char* src1, int srcWidth, int srcHeight, unsigned char*dst)
{
	unsigned char*src = (unsigned char*)malloc(srcHeight*srcWidth*sizeof(unsigned char));
	memcpy(src, src1, srcHeight*srcWidth*sizeof(unsigned char));

	IppStatus status = ippStsNoErr;

	Ipp8u pMask[3 * 3] = { 1,1,1,1,0,1,1,1,1 };
	//Ipp8u pMask[5 * 5] = { 1,1, 1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1 };
	
	IppiSize maskSize = { 3,3 };
	//IppiSize maskSize = { 5,5 };



	

	int pSpecSize;
	int pBufferSize;

	int srcStep = srcWidth*sizeof(unsigned char);
	int dstStep = srcStep;
	IppiSize roiSize = { srcWidth,srcHeight };
	IppiBorderType borderType = ippBorderRepl;
	IppiMorphState *pMorphSpec;
	Ipp8u *pBuffer = NULL;

	status = ippiMorphologyBorderGetSize_8u_C1R(roiSize, maskSize, &pSpecSize, &pBufferSize);
	if (status != ippStsNoErr)
	{
		free(src);
		return -1;
	}

	pMorphSpec = (IppiMorphState*)ippsMalloc_8u(pSpecSize);
	pBuffer = (Ipp8u*)ippsMalloc_8u(pBufferSize);

	status = ippiMorphologyBorderInit_8u_C1R(roiSize, pMask, maskSize, pMorphSpec, pBuffer);
	if (status != ippStsNoErr)
	{
		free(src);
		ippsFree(pMorphSpec);
		ippsFree(pBuffer);
		return -1;
	}

	status = ippiErodeBorder_8u_C1R(src, srcStep, dst, dstStep, roiSize, borderType, 0, pMorphSpec, pBuffer);
	if (status != ippStsNoErr)
	{
		free(src);
		ippsFree(pMorphSpec);
		ippsFree(pBuffer);
		return -1;
	}

	free(src);
	ippsFree(pMorphSpec);
	ippsFree(pBuffer);

	return 0;
}

/**********************************************/
// VisMorphologyDilation, 形态学膨胀, 窗宽(3 或 5)
// Input:
//		src		输入图像	
// Output:
//		dst		输出图像
// Return:
//		ippstatus
// Liu Ping, 2017.05.09
/**********************************************/
int CVisAlg2DBase::VisMorphologyDilation(unsigned char* src1, int srcWidth, int srcHeight, unsigned char*dst)
{
	unsigned char*src = (unsigned char*)malloc(srcHeight*srcWidth*sizeof(unsigned char));
	memcpy(src, src1, srcHeight*srcWidth*sizeof(unsigned char));

	IppStatus status = ippStsNoErr;

	Ipp8u pMask[3 * 3] = {1,1,1,1,0,1,1,1,1};
	//Ipp8u pMask[5 * 5] = { 1,1, 1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1 };
	
	IppiSize maskSize = { 3,3 };
	//IppiSize maskSize = { 5,5 };

	int pSpecSize;
	int pBufferSize;

	int srcStep = srcWidth*sizeof(unsigned char);
	int dstStep = srcStep;
	IppiSize roiSize = { srcWidth,srcHeight };
	IppiBorderType borderType = ippBorderRepl;
	IppiMorphState *pMorphSpec;
	Ipp8u *pBuffer = NULL;

	status = ippiMorphologyBorderGetSize_8u_C1R(roiSize, maskSize, &pSpecSize, &pBufferSize);
	if (status != ippStsNoErr)
	{
		free(src);
		return -1;
	}

	pMorphSpec = (IppiMorphState*)ippsMalloc_8u(pSpecSize);
	pBuffer = (Ipp8u*)ippsMalloc_8u(pBufferSize);

	status = ippiMorphologyBorderInit_8u_C1R(roiSize, pMask, maskSize, pMorphSpec, pBuffer);
	if (status != ippStsNoErr)
	{
		free(src);
		ippsFree(pMorphSpec);
		ippsFree(pBuffer);
		return -1;
	}

	status = ippiDilateBorder_8u_C1R(src, srcStep, dst, dstStep, roiSize, borderType, 0, pMorphSpec, pBuffer);
	if (status != ippStsNoErr)
	{
		free(src);
		ippsFree(pMorphSpec);
		ippsFree(pBuffer);
		return -1;
	}

	free(src);
	ippsFree(pMorphSpec);
	ippsFree(pBuffer);
	return 0;
}

/**********************************************/
// VisMorphologyDilation, 形态学open, 窗宽(3 或 5)
// Input:
//		src		输入图像	
// Output:
//		dst		输出图像
// Return:
//		ippstatus
// Liu Ping, 2017.05.09
/**********************************************/
int CVisAlg2DBase::VisMorphologyOpen(IMG_UBBUF src, IMG_UBBUF dst)
{
	__try//try//
	{
		IppStatus status = ippStsNoErr;

		IppiMorphAdvState *pSpec = NULL;
		Ipp8u* pBuffer = NULL;

		Ipp8u pMask[9] = { 1,1,1,1,0,1,1,1,1 };
		IppiSize maskSize = { 3,3 };
		//Ipp8u pMask[5 * 5] = { 1, 1, 1,1,1,1, 1, 1,1,1,
		//	1,1, 0, 1,1,
		//	1, 1, 1,1,1,1, 1, 1,1,1 };
		//IppiSize maskSize = { 5, 5 };

		Ipp8u *srcRoi = src.ptr;
		IppiSize Roi;
		Roi.width = src.size.width;
		Roi.height = src.size.height;

		Ipp8u *dstRoi;
		dstRoi = dst.ptr;

		int srcStep = Roi.width * sizeof(Ipp8u);
		int dstStep = Roi.width * sizeof(Ipp8u);
		int dstSize = Roi.width;

		int specSize = 0, bufferSize = 0;
		IppiBorderType borderType = ippBorderRepl;
		Ipp8u borderValue = 0;

		status = ippiMorphAdvGetSize_8u_C1R(Roi, maskSize, &specSize, &bufferSize);
		if (status != ippStsNoErr)
			return status;

		pSpec = (IppiMorphAdvState*)ippsMalloc_8u(specSize);
		pBuffer = (Ipp8u*)ippsMalloc_8u(bufferSize);

		status = ippiMorphAdvInit_8u_C1R(Roi, pMask, maskSize, pSpec, pBuffer);
		if (status != ippStsNoErr)
		{
			ippsFree(pBuffer);
			ippsFree(pSpec);
			return status;
		}

		status = ippiMorphOpenBorder_8u_C1R(dstRoi, srcStep, srcRoi, dstStep, Roi, borderType, borderValue, pSpec, pBuffer);

		ippsFree(pBuffer);
		ippsFree(pSpec);

		return status;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
		return _CODE_THROW;
	}
}

/**********************************************/
// VisMorphologyDilation, 形态学close, 窗宽(3 或 5)
// Input:
//		src		输入图像	
// Output:
//		dst		输出图像
// Return:
//		ippstatus
// Liu Ping, 2017.05.09
/**********************************************/
int CVisAlg2DBase::VisMorphologyClose(IMG_UBBUF src, IMG_UBBUF dst)
{
	__try//try//
	{
		IppStatus status = ippStsNoErr;

		IppiMorphAdvState *pSpec = NULL;
		Ipp8u* pBuffer = NULL;

		Ipp8u pMask[9] = { 1,1,1,1,0,1,1,1,1 };
		IppiSize maskSize = { 3,3 };
		//Ipp8u pMask[5 * 5] = { 1, 1, 1,1,1,1, 1, 1,1,1,
		//	1,1, 0, 1,1,
		//	1, 1, 1,1,1,1, 1, 1,1,1 };
		//IppiSize maskSize = { 5, 5 };

		Ipp8u *srcRoi = src.ptr;
		IppiSize Roi;
		Roi.width = src.size.width;
		Roi.height = src.size.height;

		Ipp8u *dstRoi;
		dstRoi = dst.ptr;

		int srcStep = Roi.width * sizeof(Ipp8u);
		int dstStep = Roi.width * sizeof(Ipp8u);
		int dstSize = Roi.width;

		int specSize = 0, bufferSize = 0;
		IppiBorderType borderType = ippBorderRepl;
		Ipp8u borderValue = 0;

		status = ippiMorphAdvGetSize_8u_C1R(Roi, maskSize, &specSize, &bufferSize);
		if (status != ippStsNoErr)
			return status;

		pSpec = (IppiMorphAdvState*)ippsMalloc_8u(specSize);
		pBuffer = (Ipp8u*)ippsMalloc_8u(bufferSize);

		status = ippiMorphAdvInit_8u_C1R(Roi, pMask, maskSize, pSpec, pBuffer);
		if (status != ippStsNoErr)
		{
			ippsFree(pBuffer);
			ippsFree(pSpec);
			return status;
		}

		status = ippiMorphCloseBorder_8u_C1R(srcRoi, srcStep, dstRoi, dstStep, Roi, borderType, borderValue, pSpec, pBuffer);

		ippsFree(pBuffer);
		ippsFree(pSpec);

		return status;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
		return _CODE_THROW;
	}
}



/**********************************************/
// getRegionInfo, 功能说明：计算连通域面积，长宽，极值点，重心等信息
// Input:
//     IMG_UBBUF src, 输入图像（标记完的图）
//     int regionNum, 连通域个数
//	   
// Output:
//     int *regionArea, 保存连通域的面积
//	   RegionPeakPoint *regionPeakPoint, 保存连通域x，y的极值
//     RegionWidthHeight *regionWH,  保存连通域长宽
//	   IppiPoint *regionGraCenPoint , 保存连通域重心
//
// Return:
//     0 - 正常
//     -1 - 越界异常
// Author: Jimmy Zhan 2017/3/20
/**********************************************/
int CVisAlg2DBase::getRegionInfo(IMG_UBBUF src,
	int regionNum,
	int *regionArea,
	RegionPeakPoint *regionPeakPoint,
	RegionWidthHeight *regionWH,
	IppiPoint *regionGraCenPoint)
{
	unsigned char *pMarker = src.ptr;
	IMG_SIZE roiSize = src.size;
	//init
	memset(regionArea, 0, sizeof(int) * regionNum);
	for (int i = 0; i < regionNum; i++)
	{
		regionPeakPoint[i].xMax = 0;
		regionPeakPoint[i].xMin = 0;
		regionPeakPoint[i].yMax = 0;
		regionPeakPoint[i].yMin = 0;
		regionWH[i].height = 0;
		regionWH[i].width = 0;
		regionGraCenPoint[i].x = 0;
		regionGraCenPoint[i].y = 0;
	}


	//	calculate region area and region width / height
	for (int val = 1; val < regionNum + 1; val++)
	{
		int Xmin = roiSize.width - 1;
		int Xmax = 0;
		int Ymin = roiSize.height - 1;
		int Ymax = 0;
		int xSum = 0;
		int ySum = 0;
		for (int i = 0; i < roiSize.height; i++)
		{
			for (int j = 0; j < roiSize.width; j++)
			{
				if (pMarker[i * roiSize.width + j] == val)
				{
					regionArea[val]++;
					xSum += j;
					ySum += i;
					if (j < Xmin)	Xmin = j;
					if (j > Xmax) Xmax = j;
					if (i < Ymin)	Ymin = i;
					if (i > Ymax)	Ymax = i;
				}
			}
		}

		if (Xmin < 0 || Xmax > roiSize.width - 1 || Ymin < 0 || Ymax > roiSize.height)
			return -1;

		regionPeakPoint[val].xMin = Xmin;
		regionPeakPoint[val].xMax = Xmax;
		regionPeakPoint[val].yMin = Ymin;
		regionPeakPoint[val].yMax = Ymax;

		if (Xmax > Xmin)
			regionWH[val].width = Xmax - Xmin;
		if (Ymax > Ymin)
			regionWH[val].height = Ymax - Ymin;

		if (regionArea[val] > 0)
		{
			regionGraCenPoint[val].x = xSum / regionArea[val];
			regionGraCenPoint[val].y = ySum / regionArea[val];
		}

	}

	return 0;
}


//VisHoleFill功能：填充二值图封闭孔洞
//Return 
//0  正常
//-1 异常
//Author:Jiang He/20170601
int CVisAlg2DBase::VisHoleFill(unsigned char* srcSeg, int srcWidth, int srcHeight, unsigned char* dst)
{
	//延展两个像素扩大一周
	int width = srcWidth + 4;
	int height = srcHeight + 4;
	unsigned char* DstF;
	DstF = (unsigned char*)malloc(width*height*sizeof(unsigned char));

	//四周赋值为0
	for (int j = 0; j < width; j++)
	{
		for (int i = 0; i < 2; i++)
		{
			DstF[j + i*width] = 0;
		}
		for (int t = height - 2; t < height; t++)
		{
			DstF[j + t*width] = 0;
		}
	}

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			DstF[j + i*width] = 0;
		}
		for (int t = width - 2; t < width; t++)
		{
			DstF[t + i*width] = 0;
		}
	}

	for (int i = 2; i < height - 2; i++)
	{
		for (int j = 2; j < width - 2; j++)
		{
			DstF[j + i*width] = srcSeg[(j - 2) + (i - 2)*srcWidth];
		}
	}

	IppStatus status;
	//以（0，0）为种子点，找一个连通域，赋值为255
	IppiPoint seed;
	seed.x = 0;
	seed.y = 0;

	int pBufSize = 0;
	int srcStep = width*sizeof(Ipp8u);
	IppiSize roiSize = { width,height };
	status = ippiFloodFillGetSize(roiSize, &pBufSize);
	if (status != ippStsNoErr)
	{
		free(DstF);
		return -1;
	}
	Ipp8u *pBuffer = NULL;
	pBuffer = (Ipp8u*)malloc(pBufSize*sizeof(Ipp8u));
	IppiConnectedComp pRegion;
	status = ippiFloodFill_8Con_8u_C1IR(DstF, srcStep, roiSize, seed, 255, &pRegion, pBuffer);
	free(pBuffer);
	if (status != ippStsNoErr)
	{
		free(DstF);
		return -1;
	}
	///////////////////////////////////

	//////////
	Ipp8u* DstS = (unsigned char*)malloc(srcWidth*srcHeight*sizeof(unsigned char));
	for (int i = 0; i < srcHeight; i++)
	{
		for (int j = 0; j < srcWidth; j++)
		{
			DstS[j + i*srcWidth] = DstF[j + 2 + (i + 2)*width];
		}
	}

	//取反
	int Step = srcWidth*sizeof(unsigned char);
	Ipp8u* pDst = (unsigned char*)malloc(srcWidth*srcHeight*sizeof(unsigned char));
	IppiSize Size = { srcWidth,srcHeight };
	status = ippiNot_8u_C1R(DstS, Step, pDst, Step, Size);
	if (status != ippStsNoErr)
	{
		free(DstS);
		free(DstF);
		free(pDst);
		return -1;
	}

	////pDst+srcSeg
	status = ippiAdd_8u_C1RSfs(pDst, Step, srcSeg, Step, dst, Step, Size, 0);
	if (status != ippStsNoErr)
	{
		free(DstF);
		free(pDst);
		free(DstS);
		return -1;
	}

	free(DstF);
	free(pDst);
	free(DstS);

	return 0;
}


/**********************************************/
// VisSS_LabelMarker, 功能说明:对图像中的Blob块进行提取(经过二值化的图片，背景为0，特征区域为255)。
// Input:
// 	Ipp8u *srcRoi, 输入图片
//	IppiSize Roi,感兴趣区域大小
//     ...
// Output:
// 	int &markersNum, 输出Blob块的个数
//     ...
// Return:
//     0 - 正常
//     -1 - 暂无待定
//     -2 - 暂无待定
// Author: Huang Yige/05/02/2017
/**********************************************/
void CVisAlg2DBase::VisLabelMarker(Ipp8u *srcRoi, IppiSize Roi, int &markersNum)
{
	int markerStep = Roi.width * sizeof(Ipp8u);
	int minLabel = 1;
	int maxLabel = 200;
	Ipp8u *pBuffer = NULL;
	int bufferSize = 0;

	ippiLabelMarkersGetBufferSize_8u_C1R(Roi, &bufferSize);

	pBuffer = ippsMalloc_8u(bufferSize);

	ippiLabelMarkers_8u_C1IR(srcRoi, markerStep, Roi, minLabel, maxLabel, ippiNormInf, &markersNum, pBuffer);

	ippsFree(pBuffer);
}

/**********************************************/
// VisMoment, 功能说明:计算图像的矩。
// Input:
// 	Ipp8u *srcRoi, 输入图片
//	IppiSize Roi,感兴趣区域大小
//     ...
// Output:
// 	double &hu, 图像的矩
//     ...
// Return:
//     0 - 正常
//     -1 - 暂无待定
//     -2 - 暂无待定
// Author: Huang Yige/03/18/2017
/**********************************************/
void CVisAlg2DBase::VisMoment(Ipp8u *srcRoi, IppiSize Roi, double &hu)
{
	//IppiSize Roi = { src.size.width,src.size.height };
	int sizeState = 0;
	IppiHuMoment_64f hmH = { 0 };
	IppiMomentState_64f* pState = NULL;

	ippiMomentGetStateSize_64f(ippAlgHintNone, &sizeState);

	pState = (IppiMomentState_64f*)ippMalloc(sizeState);

	ippiMomentInit_64f(pState, ippAlgHintNone);

	ippiMoments64f_8u_C1R(srcRoi, Roi.width * sizeof(Ipp8u), Roi, pState);

	ippiGetHuMoments_64f(pState, 0, hmH);

	hu = hmH[0];

	ippsFree(pState);
}





/**********************************************  Gradient *******************************************/
//计算Y方向的梯度强度,
//rowNum行和colNum列的矩形块的灰度值的和相减
//edgePoint 大小为（srcWidth-colNum+1）
//Return 
//0    正常
//-1   异常
int CVisAlg2DBase::GradientCompute(unsigned char *src, int srcWidth, int srcHeight, int rowNum, int colNum, int *dst, IMG_COORD *edgePoint)
{
	//////0、异常判断
	int w = 2 * rowNum;
	int h = colNum;
	if (src == NULL || srcWidth < w || srcHeight < h)
	{
		return -1;
	}

	//////1、边界处理（与原图相等）
	for (int i = 0; i < srcHeight*srcWidth; i++)
	{
		dst[i] = (int)src[i];
	}
	//memcpy(dst,src,srcWidth*srcHeight*sizeof(unsigned char));

	//////3、计算gradient
	int len = colNum / 2;
	int a = 0;

	IMG_COORD point = { 0,0 };
	for (int j = len; j < srcWidth - len; j++)
	{
		int maxGrad = 0;
		for (int i = rowNum; i < srcHeight - rowNum; i++)
		{
			int sum1 = 0;
			int sum2 = 0;
			for (int t = 0; t < rowNum; t++)
			{
				for (int r = 0; r < colNum; r++)
				{
					sum1 = sum1 + src[(i - t - 1)*srcWidth + j - len + r];
					sum2 = sum2 + src[(i + t + 1)*srcWidth + j - len + r];
				}
			}

			dst[j + i*srcWidth] = sum1 - sum2;
			if (sum1 - sum2 > maxGrad)
			{
				maxGrad = sum1 - sum2;
				point.x = j;
				point.y = i;
			}
		}
		edgePoint[a].x = point.x;
		edgePoint[a].y = point.y;
		a++;
	}


	return 0;
}



/**********************************************  Edge *******************************************/
int CVisAlg2DBase::GradientCompute(unsigned char *src, int srcWidth, int srcHeight, int rowNum, int colNum, float threshold, int *dst, vector<edgeInformation> &edgePoint)
{
	//////0、异常判断
	int h = 2 * rowNum;
	int w = colNum;
	if (src == NULL || srcWidth < w || srcHeight < h)
	{
		return -1;
	}

	//////1、边界处理（与原图相等）
	for (int i = 0; i < srcHeight*srcWidth; i++)
	{
		dst[i] = (int)src[i];
	}
	//memcpy(dst,src,srcWidth*srcHeight*sizeof(unsigned char));

	//////3、计算gradient
	int len = colNum / 2;
	int a = 0;
	int tempGray;

	IMG_COORD point = { 0,0 };
	for (int j = len; j < srcWidth - len; j++)
	{
		for (int i = rowNum; i < srcHeight - rowNum; i++)
		{
			int sum1 = 0;
			int sum2 = 0;
			for (int t = 0; t < rowNum; t++)
			{
				for (int r = 0; r < colNum; r++)
				{
					sum1 = sum1 + src[(i - t - 1)*srcWidth + j - len + r];
					sum2 = sum2 + src[(i + t + 1)*srcWidth + j - len + r];
				}
			}
			tempGray = sum2 - sum1;
			dst[j + i*srcWidth] = tempGray;
			/*if (tempGray > threshold)
			{
			dst[j + i*srcWidth] = tempGray;
			temp.xyDecimal = { (float)j,(float)i };
			edgePoint.push_back(temp);
			}*/
		}
	}

	/////////
	edgeInformation temp;
	int k1 = 0;
	int k2 = 0;
	int k3 = 0;
	float deci = 0;
	for (int j = 0; j < srcWidth; j++)
	{
		for (int i = 1; i < srcHeight - 1; i++)
		{
			if (dst[j + i*srcWidth]>threshold && dst[j + i*srcWidth]>dst[j + (i - 1)*srcWidth] &&
				dst[j + i*srcWidth] >= dst[j + (i + 1)*srcWidth])
			{
				if (dst[j + (i - 1)*srcWidth]>0 && dst[j + (i + 1)*srcWidth] > 0)
				{
					k1 = dst[j + (i + 1)*srcWidth];
					k2 = dst[j + i*srcWidth];
					k3 = dst[j + (i - 1)*srcWidth];
					deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));
					temp.angle = 0;
					temp.gradient = dst[j + i*srcWidth];
					temp.xyInteger.x = j;
					temp.xyInteger.y = i;
					temp.xyDecimal.x = j;
					temp.xyDecimal.y = i - deci;
					edgePoint.push_back(temp);
				}
			}
		}
	}

	return 0;
}

//边缘检测
///////////////////////////////////////////////////////////////////////////
//VisEdge_detection功能说明
//Input
//srcRoi   输入图像
//roiSize  输入图像的尺寸
//threshold  梯度强度的阈值.如果用户不知道阈值设为多少合适，可以输入0（算法自动获取合适的阈值）。
//
//output
//dstRoi  梯度强度
//edgeInformation *&edgeArray  边缘点信息，包括像素坐标、亚像素坐标、梯度强度、角度
//
//函数返回
//正常情况下返回1；
//如果用户输入阈值小于0，函数返回-1；
//如果输入参数不正确，包括图像尺寸不正确、srcRoi大小与尺寸不符合，函数返回-1。
//Author：Jiang He/20170227
////////////////////////////////////////////////
int CVisAlg2DBase::SobelFilter_8u16s_C1_5x5or3x3(IMG_UBYTE * pSrc, IppiSize roiSize, Ipp16s *pDst, Ipp32f *pAngle, int flag)
{

	IppiMaskSize mask = ippMskSize5x5;
	if (flag == 3)
	{
		mask = ippMskSize3x3;
	}
	else
	{
		if (flag == 5)
		{
			mask = ippMskSize5x5;
		}
		else
		{
			return -1;
		}
	}

	IppiBorderType bordertype = ippBorderRepl; //Border is replicated from the edge pixels
	Ipp16s *pHoriz, *pVert;

	int srcStep = roiSize.width * sizeof(Ipp8u);

	int dstStep = roiSize.width * sizeof(Ipp16s);
	int angleStep = roiSize.width * sizeof(Ipp32f);
	int bufferSize;
	int bufLen = roiSize.width * roiSize.height;
	IppStatus statusVert, statusHoriz, status;
	Ipp8u *pBuffer;
	IppNormType normType = ippNormL2;//input gradient magnitude

	pVert = (Ipp16s *)malloc(sizeof(Ipp16s)*bufLen);
	pHoriz = (Ipp16s *)malloc(sizeof(Ipp16s)*bufLen);

	ippiGradientVectorGetBufferSize(roiSize, mask, ipp16s, 1, &bufferSize);
	pBuffer = (Ipp8u *)malloc(bufferSize);
	ippiGradientVectorSobel_8u16s_C1R(pSrc, srcStep, pVert, dstStep, pHoriz, dstStep, pDst, dstStep, pAngle, angleStep, roiSize, mask, normType, bordertype, NULL, pBuffer);

	free(pVert);
	free(pHoriz);
	free(pBuffer);

	return 0;
}

//给定一个值threshold，计算最大类间方差
IMG_REAL CVisAlg2DBase::getIntraClassVariance(Ipp16s* src, int srcRows, int srcCols, int &varTh)//int &varian)
{
	//intra-class variance
	float varian = 0;

	int sumPixel = srcRows*srcCols;
	int sumGrayValue = 0;
	int average = 0;

	int sumApixel = 0;
	double PA = 0;
	int sumAgray = 0;
	int averageA = 0;

	int sumBpixel = 0;
	double PB = 0;
	int sumBgray = 0;
	int averageB = 0;

	for (int i = 0; i < sumPixel; i++)
	{
		sumGrayValue = sumGrayValue + src[i];
		if (src[i] < varTh)
		{
			sumApixel++;
			sumAgray = sumAgray + src[i];
		}
	}

	average = sumGrayValue / sumPixel;
	PA = (double)sumApixel / sumPixel;
	if (sumApixel > 0)
	{
		averageA = sumAgray / sumApixel;
	}
	else
	{
		averageA = 0;
	}

	sumBpixel = sumPixel - sumApixel;
	PB = 1.0 - PA;
	sumBgray = sumGrayValue - sumAgray;
	if (sumBpixel > 0)
	{
		averageB = sumBgray / sumBpixel;
	}
	else
	{
		averageB = 0;
	}

	//ICV = PA?(MA?M)2 + PB?(MB?M)2
	varian = PA * (averageA - average) * (averageA - average) + PB * (averageB - average) * (averageB - average);

	return varian;
}

IMG_INT CVisAlg2DBase::VisEdge_detection(IMG_UBYTE *srcRoi, IMG_SIZE roiSize, int threshold, IMG_WORD *dstRoi, IMG_UBYTE *dstRoiE, edgeInformation *&edgeArray, IMG_INT &eNum)
{
	//如果阈值小于0，函数直接返回-1
	if (threshold < 0)
	{
		return -1;
	}

	int roiRows = roiSize.height;
	int roiCols = roiSize.width;

	//角度信息
	Ipp32f *angAll;
	angAll = (Ipp32f*)malloc(roiRows*roiCols*sizeof(Ipp32f));

	std::vector<edgeInformation> edgeInfor;
	edgeInformation edInf;

	int k = 0;//记录边缘点的个数
	Ipp16u k1;//抛物线拟合的三个已知点
	Ipp16u k2;
	Ipp16u k3;
	float deci;//抛物线拟合顶点的小数部分，即对应的亚像素
	float sumx = 0;//边缘点的x坐标之和
	float sumy = 0;
	int numberChannels = 1; //the source image is single channel

	IppiSize dstRoiSize = { roiCols,roiRows };

	int sta = SobelFilter_8u16s_C1_5x5or3x3(srcRoi, dstRoiSize, dstRoi, angAll, 3);
	if (sta != 0)
	{
		free(angAll);
		return -1;
	}


	//把角度由[-PI，PI]变为[0，360]
	/*for (int i = 0; i < roiRows; i++)
	{
	for (int j = 0; j < roiCols; j++)
	{
	angAll[j + i * roiCols] = (float)(180 - angAll[j + i * roiCols] / PI * 180);
	}
	}*/

	/*
	FILE *sx;
	sx = fopen("E:\\ProjectMy\\a1sobel.txt", "w");
	FILE *ang;
	ang = fopen("E:\\ProjectMy\\a1ang.txt", "w");
	for (int i = 0; i<roiRows; i++)
	{
	for (int j = 0; j < roiCols; j++)
	{
	fprintf(sx, "%d   ", dstRoi[j+i*roiCols]);
	fprintf(ang, "%f   ", angAll[j+i*roiCols]);
	}
	fprintf(sx,"\n");
	fprintf(ang, "\n");
	}
	fclose(sx);
	fclose(ang);
	*/
	//自动获取梯度强度的阈值
	//什么情况下算是没有输入阈值呢？
	if (threshold == 0)
	{
		//Otsu法，遍历所有的灰度值，从1到255，使intra-class invariance最大的那个值，即为要求的阈值
		int varian = 0;
		int temp = 0;
		for (int p = 1; p < 800; p++)
		{
			temp = getIntraClassVariance(dstRoi, roiRows, roiCols, p);
			if (varian < temp)
			{
				varian = temp;
				threshold = p;
			}
		}
	}
	//printf("%d\n",threshold);

	//到亚像素
	for (int i = 1; i<roiRows - 1; i++)
	{
		for (int j = 1; j<roiCols - 1; j++)
		{
			if (dstRoi[j + i*roiCols] > threshold)
			{
				angAll[j + i * roiCols] = (float)180 - angAll[j + i * roiCols] / PI * 180;
				if ((angAll[j + i*roiCols]>22.5) && (angAll[j + i*roiCols] < 67.5))
				{
					if ((dstRoi[j + i*roiCols] > dstRoi[j - 1 + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + (i + 1)*roiCols]))
					{
						k1 = dstRoi[j - 1 + (i - 1)*roiCols];
						k2 = dstRoi[j + i*roiCols];
						k3 = dstRoi[(j + 1) + (i + 1)*roiCols];
						deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));


						edInf.xyInteger.x = j;
						edInf.xyInteger.y = i;
						edInf.xyDecimal.x = j + deci;
						edInf.xyDecimal.y = i + deci;
						edInf.gradient = dstRoi[j + i*roiCols];
						edInf.angle = angAll[j + i*roiCols];
						edgeInfor.push_back(edInf);
						k++;
					}
				}
				else
				{
					if ((angAll[j + i*roiCols] > 202.5) && (angAll[j + i*roiCols] < 247.5))
					{
						if ((dstRoi[j + i*roiCols] > dstRoi[j - 1 + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + (i + 1)*roiCols]))
						{
							k3 = dstRoi[j - 1 + (i - 1)*roiCols];
							k2 = dstRoi[j + i*roiCols];
							k1 = dstRoi[(j + 1) + (i + 1)*roiCols];
							deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

							edInf.xyInteger.x = j;
							edInf.xyInteger.y = i;
							edInf.xyDecimal.x = j - deci;
							edInf.xyDecimal.y = i - deci;
							edInf.gradient = dstRoi[j + i*roiCols];
							edInf.angle = angAll[j + i*roiCols];
							edgeInfor.push_back(edInf);
							k++;
						}
					}
					else
					{
						if ((angAll[j + i*roiCols] > 112.5) && (angAll[j + i*roiCols] < 157.5))
						{

							if ((dstRoi[j + i*roiCols] > dstRoi[(j + 1) + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j - 1) + (i + 1)*roiCols]))
							{
								k1 = dstRoi[(j + 1) + (i - 1)*roiCols];
								k2 = dstRoi[j + i*roiCols];
								k3 = dstRoi[(j - 1) + (i + 1)*roiCols];
								deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

								edInf.xyInteger.x = j;
								edInf.xyInteger.y = i;
								edInf.xyDecimal.x = j - deci;
								edInf.xyDecimal.y = i + deci;
								edInf.gradient = dstRoi[j + i*roiCols];
								edInf.angle = angAll[j + i*roiCols];
								edgeInfor.push_back(edInf);
								k++;
							}
						}
						else
						{
							if ((angAll[j + i*roiCols] > 292.5) && (angAll[j + i*roiCols] < 337.5))
							{
								if ((dstRoi[j + i*roiCols] > dstRoi[(j + 1) + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j - 1) + (i + 1)*roiCols]))
								{
									k3 = dstRoi[(j + 1) + (i - 1)*roiCols];
									k2 = dstRoi[j + i*roiCols];
									k1 = dstRoi[(j - 1) + (i + 1)*roiCols];
									deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

									edInf.xyInteger.x = j;
									edInf.xyInteger.y = i;
									edInf.xyDecimal.x = j + deci;
									edInf.xyDecimal.y = i - deci;
									edInf.gradient = dstRoi[j + i*roiCols];
									edInf.angle = angAll[j + i*roiCols];
									edgeInfor.push_back(edInf);
									k++;
								}
							}
							else
							{
								if (((angAll[j + i*roiCols] >= -1) && (angAll[j + i*roiCols] <= 22.5)) || ((angAll[j + i*roiCols] >= 337.5) && (angAll[j + i*roiCols] <= 361)))
								{
									if ((dstRoi[j + i*roiCols] > dstRoi[(j - 1) + i*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + i*roiCols]))
									{
										k1 = dstRoi[(j - 1) + i*roiCols];
										k2 = dstRoi[j + i*roiCols];
										k3 = dstRoi[(j + 1) + i*roiCols];
										deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

										edInf.xyInteger.x = j;
										edInf.xyInteger.y = i;
										edInf.xyDecimal.x = j + deci;
										edInf.xyDecimal.y = i;
										edInf.gradient = dstRoi[j + i*roiCols];
										edInf.angle = angAll[j + i*roiCols];
										edgeInfor.push_back(edInf);
										k++;
									}
								}
								else
								{
									if ((angAll[j + i*roiCols] <= 202.5) && (angAll[j + i*roiCols] >= 157.5))
									{
										if ((dstRoi[j + i*roiCols] > dstRoi[(j - 1) + i*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + i*roiCols]))
										{
											k3 = dstRoi[(j - 1) + i*roiCols];
											k2 = dstRoi[j + i*roiCols];
											k1 = dstRoi[(j + 1) + i*roiCols];
											deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

											edInf.xyInteger.x = j;
											edInf.xyInteger.y = i;
											edInf.xyDecimal.x = j - deci;
											edInf.xyDecimal.y = i;
											edInf.gradient = dstRoi[j + i*roiCols];
											edInf.angle = angAll[j + i*roiCols];
											edgeInfor.push_back(edInf);
											k++;
										}
									}
									else
									{
										if ((angAll[j + i*roiCols] >= 67.5) && (angAll[j + i*roiCols] <= 112.5))
										{

											if ((dstRoi[j + i*roiCols] > dstRoi[j + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[j + (i + 1)*roiCols]))
											{
												k1 = dstRoi[j + (i - 1)*roiCols];
												k2 = dstRoi[j + i*roiCols];
												k3 = dstRoi[j + (i + 1)*roiCols];
												deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

												edInf.xyInteger.x = j;
												edInf.xyInteger.y = i;
												edInf.xyDecimal.x = j;
												edInf.xyDecimal.y = i + deci;
												edInf.gradient = dstRoi[j + i*roiCols];
												edInf.angle = angAll[j + i*roiCols];
												edgeInfor.push_back(edInf);
												k++;
											}
										}
										else
										{
											if ((angAll[j + i*roiCols] >= 247.5) && (angAll[j + i*roiCols] <= 292.5))
											{
												if ((dstRoi[j + i*roiCols] > dstRoi[j + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[j + (i + 1)*roiCols]))
												{
													k3 = dstRoi[j + (i - 1)*roiCols];
													k2 = dstRoi[j + i*roiCols];
													k1 = dstRoi[j + (i + 1)*roiCols];
													deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

													edInf.xyInteger.x = j;
													edInf.xyInteger.y = i;
													edInf.xyDecimal.x = j;
													edInf.xyDecimal.y = i - deci;
													edInf.gradient = dstRoi[j + i*roiCols];
													edInf.angle = angAll[j + i*roiCols];
													edgeInfor.push_back(edInf);
													k++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}


	eNum = k;
	//二值图
	for (int t = 0; t < roiCols*roiRows; t++)//二值图像，所有像素先都赋值为0，边缘点赋值255
	{
		dstRoiE[t] = 0;
	}
	for (int q = 0; q < k; q++)
	{
		dstRoiE[edgeInfor[q].xyInteger.x + edgeInfor[q].xyInteger.y * roiCols] = 255;
	}
	/*
	FILE *Binary;
	Binary = fopen("E:\\ProjectMy\\Binary.txt", "w");
	for (int i = 0; i<roiRows; i++)
	{
	for (int j = 0; j < roiCols; j++)
	{
	fprintf(Binary, "%d   ", dstRoiE[j + i*roiCols]);
	}
	fprintf(Binary, "\n");
	}
	fclose(Binary);
	*/
	//以数组的方式传出边缘信息
	edgeArray = (edgeInformation*)malloc(k*sizeof(edgeInformation));
	for (int i = 0; i < k; i++)
	{
		edgeArray[i] = edgeInfor[i];
	}


	/*FILE *e;
	e = fopen("E:\\project03\\e.txt", "w");
	for (int i = 0; i < k; i++)
	{
	fprintf(e,"%d   %d   \n", edgeInfor[i].xInteger, edgeInfor[i].yInteger);
	}
	*/

	free(angAll);
	return 1;
}


