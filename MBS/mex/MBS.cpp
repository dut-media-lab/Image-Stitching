/*****************************************************************************
*	Implemetation of the saliency detction method described in paper
*	"Minimum Barrier Salient Object Detection at 80 FPS", Jianming Zhang, 
*	Stan Sclaroff, Zhe Lin, Xiaohui Shen, Brian Price, Radomir Mech, ICCV, 
*       2015
*	
*	Copyright (C) 2015 Jianming Zhang
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*	If you have problems about this software, please contact: 
*       jimmie33@gmail.com
*******************************************************************************/

#include <iostream>
#include <ctime>
#include "opencv2/opencv.hpp"
#include "mexopencv.hpp"
#include "MBS.hpp"
#include <vector>
#include <cmath>

using namespace std;
using namespace cv;

#define MAX_IMG_DIM 300
#define TOLERANCE 0.01
#define FRAME_MAX 20
#define SOBEL_THRESH 0.4

MBS::MBS(const Mat& src)
:mAttMapCount(0)
{
	mSrc=src.clone();
	mSaliencyMap = Mat::zeros(src.size(), CV_32FC1);

	split(mSrc, mFeatureMaps);

	for (int i = 0; i < mFeatureMaps.size(); i++)
	{
		//normalize(mFeatureMaps[i], mFeatureMaps[i], 255.0, 0.0, NORM_MINMAX);
		medianBlur(mFeatureMaps[i], mFeatureMaps[i], 5);
	}
}

void MBS::computeSaliency(bool use_geodesic)
{

        if (use_geodesic)
		mMBSMap = fastGeodesic(mFeatureMaps);
	else
		mMBSMap = fastMBS(mFeatureMaps);
	normalize(mMBSMap, mMBSMap, 0.0, 1.0, NORM_MINMAX);
	mSaliencyMap = mMBSMap;
}

Mat MBS::getSaliencyMap()
{
	Mat ret;
	normalize(mSaliencyMap, ret, 0.0, 255.0, NORM_MINMAX);
	ret.convertTo(ret, CV_8UC1);
	return ret;
}


void rasterScan(const Mat& featMap, Mat& map, Mat& lb, Mat& ub)
{
	Size sz = featMap.size();
	float *pMapup = (float*)map.data + 1;
	float *pMap = pMapup + sz.width;
	uchar *pFeatup = featMap.data + 1;
	uchar *pFeat = pFeatup + sz.width;
	uchar *pLBup = lb.data + 1;
	uchar *pLB = pLBup + sz.width;
	uchar *pUBup = ub.data + 1;
	uchar *pUB = pUBup + sz.width;

	float mapPrev;
	float featPrev;
	uchar lbPrev, ubPrev;

	float lfV, upV;
	int flag;
	for (int r = 1; r < sz.height - 1; r++)
	{
		mapPrev = *(pMap - 1);
		featPrev = *(pFeat - 1);
		lbPrev = *(pLB - 1);
		ubPrev = *(pUB - 1);


		for (int c = 1; c < sz.width - 1; c++)
		{
			lfV = MAX(*pFeat, ubPrev) - MIN(*pFeat, lbPrev);//(*pFeat >= lbPrev && *pFeat <= ubPrev) ? mapPrev : mapPrev + abs((float)(*pFeat) - featPrev);
			upV = MAX(*pFeat, *pUBup) - MIN(*pFeat, *pLBup);//(*pFeat >= *pLBup && *pFeat <= *pUBup) ? *pMapup : *pMapup + abs((float)(*pFeat) - (float)(*pFeatup));

			flag = 0;
			if (lfV < *pMap)
			{
				*pMap = lfV;
				flag = 1;
			}
			if (upV < *pMap)
			{
				*pMap = upV;
				flag = 2;
			}

			switch (flag)
			{
			case 0:		// no update
				break;
			case 1:		// update from left
				*pLB = MIN(*pFeat, lbPrev);
				*pUB = MAX(*pFeat, ubPrev);
				break;
			case 2:		// update from up
				*pLB = MIN(*pFeat, *pLBup);
				*pUB = MAX(*pFeat, *pUBup);
				break;
			default:
				break;
			}

			mapPrev = *pMap;
			pMap++; pMapup++;
			featPrev = *pFeat;
			pFeat++; pFeatup++;
			lbPrev = *pLB;
			pLB++; pLBup++;
			ubPrev = *pUB;
			pUB++; pUBup++;
		}
		pMapup += 2; pMap += 2;
		pFeat += 2; pFeatup += 2;
		pLBup += 2; pLB += 2;
		pUBup += 2; pUB += 2;
	}
}

void invRasterScan(const Mat& featMap, Mat& map, Mat& lb, Mat& ub)
{
	Size sz = featMap.size();
	int datalen = sz.width*sz.height;
	float *pMapdn = (float*)map.data + datalen - 2;
	float *pMap = pMapdn - sz.width;
	uchar *pFeatdn = featMap.data + datalen - 2;
	uchar *pFeat = pFeatdn - sz.width;
	uchar *pLBdn = lb.data + datalen - 2;
	uchar *pLB = pLBdn - sz.width;
	uchar *pUBdn = ub.data + datalen - 2;
	uchar *pUB = pUBdn - sz.width;
	
	float mapPrev;
	float featPrev;
	uchar lbPrev, ubPrev;

	float rtV, dnV;
	int flag;
	for (int r = 1; r < sz.height - 1; r++)
	{
		mapPrev = *(pMap + 1);
		featPrev = *(pFeat + 1);
		lbPrev = *(pLB + 1);
		ubPrev = *(pUB + 1);

		for (int c = 1; c < sz.width - 1; c++)
		{
			rtV = MAX(*pFeat, ubPrev) - MIN(*pFeat, lbPrev);//(*pFeat >= lbPrev && *pFeat <= ubPrev) ? mapPrev : mapPrev + abs((float)(*pFeat) - featPrev);
			dnV = MAX(*pFeat, *pUBdn) - MIN(*pFeat, *pLBdn);//(*pFeat >= *pLBdn && *pFeat <= *pUBdn) ? *pMapdn : *pMapdn + abs((float)(*pFeat) - (float)(*pFeatdn));

			flag = 0;
			if (rtV < *pMap)
			{
				*pMap = rtV;
				flag = 1;
			}
			if (dnV < *pMap)
			{
				*pMap = dnV;
				flag = 2;
			}

			switch (flag)
			{
			case 0:		// no update
				break;
			case 1:		// update from right
				*pLB = MIN(*pFeat, lbPrev);
				*pUB = MAX(*pFeat, ubPrev);
				break;
			case 2:		// update from down
				*pLB = MIN(*pFeat, *pLBdn);
				*pUB = MAX(*pFeat, *pUBdn);
				break;
			default:
				break;
			}

			mapPrev = *pMap;
			pMap--; pMapdn--;
			featPrev = *pFeat;
			pFeat--; pFeatdn--;
			lbPrev = *pLB;
			pLB--; pLBdn--;
			ubPrev = *pUB;
			pUB--; pUBdn--;
		}


		pMapdn -= 2; pMap -= 2;
		pFeatdn -= 2; pFeat -= 2;
		pLBdn -= 2; pLB -= 2;
		pUBdn -= 2; pUB -= 2;
	}
}

cv::Mat fastMBS(const std::vector<cv::Mat> featureMaps)
{
	assert(featureMaps[0].type() == CV_8UC1);

	Size sz = featureMaps[0].size();
	Mat ret = Mat::zeros(sz, CV_32FC1);
	if (sz.width < 3 || sz.height < 3)
		return ret;

	for (int i = 0; i < featureMaps.size(); i++)
	{
		Mat map = Mat::zeros(sz, CV_32FC1);
		Mat mapROI(map, Rect(1, 1, sz.width - 2, sz.height - 2));
		mapROI.setTo(Scalar(100000));
		Mat lb = featureMaps[i].clone();
		Mat ub = featureMaps[i].clone();

		rasterScan(featureMaps[i], map, lb, ub);
		invRasterScan(featureMaps[i], map, lb, ub);
		rasterScan(featureMaps[i], map, lb, ub);
		
		ret += map;
	}

	return ret;
	
}

float getThreshForGeo(const Mat& src)
{
	float ret;
	Size sz = src.size();

	uchar *pFeatup = src.data + 1;
	uchar *pFeat = pFeatup + sz.width;
	uchar *pfeatdn = pFeat + sz.width;

	float featPrev;

	for (int r = 1; r < sz.height - 1; r++)
	{
		featPrev = *(pFeat - 1);

		for (int c = 1; c < sz.width - 1; c++)
		{
			float temp = MIN(abs(*pFeat-featPrev),abs(*pFeat-*(pFeat+1)));
			temp = MIN(temp,abs(*pFeat-*pFeatup));
			temp = MIN(temp,abs(*pFeat-*pfeatdn));
			ret += temp;

			featPrev = *pFeat;
			pFeat++; pFeatup++; pfeatdn++;
		}
		pFeat += 2; pFeatup += 2; pfeatdn += 2;
	}
	return ret / ((sz.width - 2)*(sz.height - 2));
}

void rasterScanGeo(const Mat& featMap, Mat& map, float thresh)
{
	Size sz = featMap.size();
	float *pMapup = (float*)map.data + 1;
	float *pMap = pMapup + sz.width;
	uchar *pFeatup = featMap.data + 1;
	uchar *pFeat = pFeatup + sz.width;

	float mapPrev;
	float featPrev;

	float lfV, upV;
	int flag;
	for (int r = 1; r < sz.height - 1; r++)
	{
		mapPrev = *(pMap - 1);
		featPrev = *(pFeat - 1);


		for (int c = 1; c < sz.width - 1; c++)
		{
			lfV = (abs(featPrev - *pFeat)>thresh ? abs(featPrev - *pFeat):0.0f) + mapPrev;
			upV = (abs(*pFeatup - *pFeat)>thresh ? abs(*pFeatup - *pFeat):0.0f) + *pMapup;
			
			if (lfV < *pMap)
			{
				*pMap = lfV;
			}
			if (upV < *pMap)
			{
				*pMap = upV;
			}

			mapPrev = *pMap;
			pMap++; pMapup++;
			featPrev = *pFeat;
			pFeat++; pFeatup++;
		}
		pMapup += 2; pMap += 2;
		pFeat += 2; pFeatup += 2;
	}
}

void invRasterScanGeo(const Mat& featMap, Mat& map, float thresh)
{
	Size sz = featMap.size();
	int datalen = sz.width*sz.height;
	float *pMapdn = (float*)map.data + datalen - 2;
	float *pMap = pMapdn - sz.width;
	uchar *pFeatdn = featMap.data + datalen - 2;
	uchar *pFeat = pFeatdn - sz.width;

	float mapPrev;
	float featPrev;

	float rtV, dnV;
	int flag;
	for (int r = 1; r < sz.height - 1; r++)
	{
		mapPrev = *(pMap + 1);
		featPrev = *(pFeat + 1);

		for (int c = 1; c < sz.width - 1; c++)
		{
			rtV = (abs(featPrev - *pFeat)>thresh ? abs(featPrev - *pFeat):0.0f) + mapPrev;
			dnV = (abs(*pFeatdn - *pFeat)>thresh ? abs(*pFeatdn - *pFeat):0.0f) + *pMapdn;
			
			if (rtV < *pMap)
			{
				*pMap = rtV;
			}
			if (dnV < *pMap)
			{
				*pMap = dnV;
			}

			mapPrev = *pMap;
			pMap--; pMapdn--;
			featPrev = *pFeat;
			pFeat--; pFeatdn--;
		}


		pMapdn -= 2; pMap -= 2;
		pFeatdn -= 2; pFeat -= 2;
	}
}

cv::Mat fastGeodesic(const std::vector<cv::Mat> featureMaps)
{
	assert(featureMaps[0].type() == CV_8UC1);

	Size sz = featureMaps[0].size();
	Mat ret = Mat::zeros(sz, CV_32FC1);
	if (sz.width < 3 || sz.height < 3)
		return ret;


	for (int i = 0; i < featureMaps.size(); i++)
	{
		// determines the threshold for clipping
		float thresh = getThreshForGeo(featureMaps[i]);
		//cout << thresh << endl;
		Mat map = Mat::zeros(sz, CV_32FC1);
		Mat mapROI(map, Rect(1, 1, sz.width - 2, sz.height - 2));
		mapROI.setTo(Scalar(1000000000));

		rasterScanGeo(featureMaps[i], map, thresh);
		invRasterScanGeo(featureMaps[i], map, thresh);
		rasterScanGeo(featureMaps[i], map, thresh);

		ret += map;
	}

	return ret;

}

int findFrameMargin(const Mat& img, bool reverse)
{
	Mat edgeMap, edgeMapDil, edgeMask;
	Sobel(img, edgeMap, CV_16SC1, 0, 1);
	edgeMap = abs(edgeMap);
	edgeMap.convertTo(edgeMap, CV_8UC1);
	edgeMask = edgeMap < (SOBEL_THRESH * 255.0);
	dilate(edgeMap, edgeMapDil, Mat(), Point(-1, -1), 2);
	edgeMap = edgeMap == edgeMapDil;
	edgeMap.setTo(Scalar(0.0), edgeMask);


	if (!reverse)
	{
		for (int i = edgeMap.rows - 1; i >= 0; i--)
			if (mean(edgeMap.row(i))[0] > 0.6*255.0)
				return i + 1;
	}
	else
	{
		for (int i = 0; i < edgeMap.rows; i++)
			if (mean(edgeMap.row(i))[0] > 0.6*255.0)
				return edgeMap.rows - i;
	}

	return 0;
}

bool removeFrame(const cv::Mat& inImg, cv::Mat& outImg, cv::Rect &roi)
{
	if (inImg.rows < 2 * (FRAME_MAX + 3) || inImg.cols < 2 * (FRAME_MAX + 3))
	{
		roi = Rect(0, 0, inImg.cols, inImg.rows);
		outImg = inImg;
		return false;
	}

	Mat imgGray;
	cvtColor(inImg, imgGray, CV_RGB2GRAY);

	int up, dn, lf, rt;
	
	up = findFrameMargin(imgGray.rowRange(0, FRAME_MAX), false);
	dn = findFrameMargin(imgGray.rowRange(imgGray.rows - FRAME_MAX, imgGray.rows), true);
	lf = findFrameMargin(imgGray.colRange(0, FRAME_MAX).t(), false);
	rt = findFrameMargin(imgGray.colRange(imgGray.cols - FRAME_MAX, imgGray.cols).t(), true);

	int margin = MAX(up, MAX(dn, MAX(lf, rt)));
	if ( margin == 0 )
	{
		roi = Rect(0, 0, imgGray.cols, imgGray.rows);
		outImg = inImg;
		return false;
	}

	int count = 0;
	count = up == 0 ? count : count + 1;
	count = dn == 0 ? count : count + 1;
	count = lf == 0 ? count : count + 1;
	count = rt == 0 ? count : count + 1;

	// cut four border region if at least 2 border frames are detected
	if (count > 1)
	{
		margin += 2;
		roi = Rect(margin, margin, inImg.cols - 2*margin, inImg.rows - 2*margin);
		outImg = Mat(inImg, roi);

		return true;
	}

	// otherwise, cut only one border
	up = up == 0 ? up : up + 2;
	dn = dn == 0 ? dn : dn + 2;
	lf = lf == 0 ? lf : lf + 2;
	rt = rt == 0 ? rt : rt + 2;

	
	roi = Rect(lf, up, inImg.cols - lf - rt, inImg.rows - up - dn);
	outImg = Mat(inImg, roi);

	return true;
	
}

Mat doWork(
	const Mat& src,
	bool use_lab,
        bool remove_border,
	bool use_geodesic
	)
{
	Mat src_small;
	float w = (float)src.cols, h = (float)src.rows;
	float maxD = max(w,h);
	resize(src,src_small,Size((int)(MAX_IMG_DIM*w/maxD),(int)(MAX_IMG_DIM*h/maxD)),0.0,0.0,INTER_AREA);// standard: width: 300 pixel
	Mat srcRoi;
	Rect roi;
	// detect and remove the artifical frame of the image
        if (remove_border)
		removeFrame(src_small, srcRoi, roi);
	else 
	{
		srcRoi = src_small;
		roi = Rect(0, 0, src_small.cols, src_small.rows);
	}


	if (use_lab)        
		cvtColor(srcRoi, srcRoi, CV_RGB2Lab);
		
	/* Computing saliency */
	MBS mbs(srcRoi);
	mbs.computeSaliency(use_geodesic);
		
	Mat resultRoi=mbs.getSaliencyMap();
	Mat result = Mat::zeros(src_small.size(), CV_32FC1);

	normalize(resultRoi, Mat(result, roi), 0.0, 1.0, NORM_MINMAX);
		
	resize(result, result, src.size());
	return result;
}


/**
 * Main entry called from Matlab
 * @param nlhs number of left-hand-side arguments
 * @param plhs pointers to mxArrays in the left-hand-side
 * @param nrhs number of right-hand-side arguments
 * @param prhs pointers to mxArrays in the right-hand-side
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    // Check the number of arguments
    if (nrhs<1 || nrhs > 4 || nlhs>1)
        mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
    
    // Decide second arguments
    bool use_geodesic = false;
    bool use_lab = true;
    bool remove_border = true;
    if (nrhs >= 2)
	use_lab = mxIsLogicalScalarTrue(prhs[1]);
    if (nrhs >= 3)
	remove_border = mxIsLogicalScalarTrue(prhs[2]);
    if (nrhs >= 4)
	use_geodesic = mxIsLogicalScalarTrue(prhs[3]);
    
    // Apply
    Mat src(((MxArray)prhs[0]).toMat()), dst;
    dst = doWork(src, use_lab, remove_border, use_geodesic);
    plhs[0] = MxArray(dst);
}
