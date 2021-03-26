/* 
Author: Bin Fan (NLPR), bfan@nlpr.ia.ac.cn

This is a demo about how to use the 'LineMatcher' library, which conducts line matching by affine invariants. 
Please refer to the following paper for more details about the algorithm.

Bin Fan, Fuchao Wu and Zhanyi Hu. Line Matching Leveraged by Point Correspondences, In CVPR 2010, pp 390-391.

The line segmentation detector (LSD) used in the demo is from the website:
http://www.ipol.im/pub/algo/gjmr_line_segment_detector/

Rafael Grompone von Gioi, J¨¦r¨¦mie Jakubowicz, Jean-Michel Morel, Gregory Randall, LSD: A Fast Line Segment Detector with a False Detection Control, IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 32, no. 4, pp. 722-732, April 2010.

The matched points in the demo are obtained based on the opencv 2.2

If you have questions, please contact: bfan@nlpr.ia.ac.cn
*/

extern "C"
{
	#include "lsd.h"
};
#include <stdio.h>
#include <highgui.h>
#include <cv.h>
#include "LineMatcher.h"
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace cv;

#pragma comment(lib,"LineMatcher.lib")

#pragma comment(lib,"opencv_features2d220.lib")
#pragma comment(lib,"opencv_highgui220.lib")
#pragma comment(lib,"opencv_calib3d220.lib")
#pragma comment(lib,"opencv_imgproc220.lib")
#pragma comment(lib,"opencv_core220.lib")



bool GetImGray(IplImage* im, double x, double y, double &gray)
{
	if (x < 1 || x > (im->width-2) || y < 1 || y > (im->height-2))
		return false;
	int x1 = (int)x;
	int y1 = (int)y;
	int x2 = x1 + 1;
	int y2 = y1 + 1;
	int step = im->widthStep;
	int _x1 = x1 < 0 ? 0 : x1;
	int _y1 = y1 < 0 ? 0 : y1;
	int _x2 = x2 > (im->width-1) ? (im->width-1) : x2;
	int _y2 = y2 > (im->height-1) ? (im->height-1) : y2;
	gray = (x2 - x) * (y2 - y) * ((BYTE*)im->imageData)[_y1*step+_x1]/1.0 + 
		   (x - x1) * (y2 - y) * ((BYTE*)im->imageData)[_y1*step+_x2]/1.0 + 
		   (x2 - x) * (y - y1) * ((BYTE*)im->imageData)[_y2*step+_x1]/1.0 + 
		   (x - x1) * (y - y1) * ((BYTE*)im->imageData)[_y2*step+_x2]/1.0;
	return true;
}

Line* LineDetect(char* imfile, int &nLines)
{
	IplImage* im = cvLoadImage(imfile,CV_LOAD_IMAGE_GRAYSCALE);
	image_double image = new_image_double(im->width, im->height);
	unsigned char* im_src = (unsigned char*)im->imageData;
	int xsize = image->xsize;
	int ysize = image->ysize;
	int y,x;
	for (y = 0;y < ysize;y++)
	{
		for (x = 0;x < xsize;x++)
		{
			image->data[y * xsize + x] = im_src[y * im->widthStep + x];
		}
	}
	ntuple_list detected_lines = lsd(image);
	free_image_double(image);

	nLines = detected_lines->size;

	IplImage* smoothed_im = cvCreateImage(cvSize(im->width,im->height),IPL_DEPTH_8U,1);
	cvSmooth(im,smoothed_im);
	cvReleaseImage(&im);

	Line* pLines = new Line[nLines];
	int nCount = 0;
	int i,j;
	int dim = detected_lines->dim;
	for (i = 0;i < nLines;i++)
	{
		double x1 = detected_lines->values[i*dim+0];
		double y1 = detected_lines->values[i*dim+1];
		double x2 = detected_lines->values[i*dim+2];
		double y2 = detected_lines->values[i*dim+3];
		double len = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
		len = sqrt(len);
		if(len <= 20) continue;
		//if(len<0.04*(xsize+ysize)/2) continue;
		double a = y1 - y2;
		double b = -(x1 - x2);
		double c = x1 * y2 - x2 * y1;
		vector<CvPoint2D32f> Pts;
		int nPts = 0;
		CvPoint2D32f p;
		p.x = x1;
		p.y = y1;
		Pts.push_back(p);
		nPts++;
		if ( abs((a+0.0000000000001) / (b+0.0000000000001)) < 1 )
		{
			double iterval = (x2 - x1) / len;
			while (1)
			{
				p.x += iterval;
				if( (iterval > 0 && p.x > x2) || (iterval < 0 && p.x < x2) ) break;
				p.y = -(a * p.x + c) / b;
				Pts.push_back(p);
				nPts++;
			}
		}
		else
		{
			double iterval = (y2 - y1) / len;
			while (1)
			{
				p.y += iterval;
				if( (iterval > 0 && p.y > y2) || (iterval < 0 && p.y < y2) ) break;
				p.x = -(b * p.y + c) / a;
				Pts.push_back(p);
				nPts++;
			}
		}
		double line_grad_x, line_grad_y;
		int nValid = 0;
		line_grad_x = line_grad_y = 0;
		for (j = 0;j < nPts;j++)
		{
			double gray2;
			double gray1;
			if( !GetImGray(smoothed_im, Pts[j].x + 1, Pts[j].y, gray2) ) continue;
			if( !GetImGray(smoothed_im, Pts[j].x - 1, Pts[j].y, gray1) ) continue;
			line_grad_x += (gray2 - gray1) / 2;
			if( !GetImGray(smoothed_im, Pts[j].x, Pts[j].y + 1, gray2) ) continue;
			if( !GetImGray(smoothed_im, Pts[j].x, Pts[j].y - 1, gray1) ) continue;
			line_grad_y += (gray2 - gray1) / 2;
			nValid++;
		}
		if( nValid == 0) continue;
		line_grad_x /= nValid;
		line_grad_y /= nValid;
		double ExProd = (Pts[nPts-1].x - Pts[0].x) * line_grad_y - line_grad_x * (Pts[nPts-1].y - Pts[0].y);
		if (ExProd > 0)
		{
			pLines[nCount].StartPt = Pts[0];
			pLines[nCount].EndPt = Pts[nPts-1];
		}
		else
		{
			pLines[nCount].StartPt = Pts[nPts-1];
			pLines[nCount].EndPt = Pts[0];
		}
		pLines[nCount].length = len;
		pLines[nCount].Center.x = (x1 + x2) / 2;
		pLines[nCount].Center.y = (y1 + y2) / 2;
		pLines[nCount].para_a = a;
		pLines[nCount].para_b = b;
		pLines[nCount].para_c = c;
		nCount++;
	}
	cvReleaseImage(&smoothed_im);
	
	free_ntuple_list(detected_lines);
	
	nLines = nCount;
	return pLines;
}

void crossCheckMatching( Ptr<DescriptorMatcher>& descriptorMatcher,
                         const Mat& descriptors1, const Mat& descriptors2,
                         vector<DMatch>& filteredMatches12, int knn=1 )
{
    filteredMatches12.clear();
    vector<vector<DMatch> > matches12, matches21;
    descriptorMatcher->knnMatch( descriptors1, descriptors2, matches12, knn );
    descriptorMatcher->knnMatch( descriptors2, descriptors1, matches21, knn );
    for( size_t m = 0; m < matches12.size(); m++ )
    {
        bool findCrossCheck = false;
        for( size_t fk = 0; fk < matches12[m].size(); fk++ )
        {
            DMatch forward = matches12[m][fk];

            for( size_t bk = 0; bk < matches21[forward.trainIdx].size(); bk++ )
            {
                DMatch backward = matches21[forward.trainIdx][bk];
                if( backward.trainIdx == forward.queryIdx )
                {
                    filteredMatches12.push_back(forward);
                    findCrossCheck = true;
                    break;
                }
            }
            if( findCrossCheck ) break;
        }
    }
}


int PtMatch(char* imfile1, char* imfile2, vector<MatchPt> &match)
{
    Ptr<FeatureDetector> detector = FeatureDetector::create( "SIFT" );
    Ptr<DescriptorExtractor> descriptorExtractor = DescriptorExtractor::create( "SIFT" );
    Ptr<DescriptorMatcher> descriptorMatcher = DescriptorMatcher::create( "BruteForce" );
    if( detector.empty() || descriptorExtractor.empty() || descriptorMatcher.empty()  )
    {
        cout << "Can not create detector or descriptor exstractor or descriptor matcher of given types" << endl;
        return -1;
    }
	
    Mat img1 = imread( imfile1 );
    Mat img2 = imread( imfile2 );
    if( img1.empty() || img2.empty() )
    {
        cout << "Can not read images" << endl;
        return -1;
    }

    cout << endl << "Extracting keypoints from first image..." << endl;

    vector<KeyPoint> keypoints1;

    detector->detect( img1, keypoints1 );
    cout << keypoints1.size() << " points" << endl;

    cout << "Computing descriptors for keypoints from first image..." << endl;
    Mat descriptors1;
    descriptorExtractor->compute( img1, keypoints1, descriptors1 );

    cout << endl << "Extracting keypoints from second image..." << endl;
    vector<KeyPoint> keypoints2;
    detector->detect( img2, keypoints2 );
    cout << keypoints2.size() << " points" << endl;

    cout << "Computing descriptors for keypoints from second image..." << endl;
    Mat descriptors2;
    descriptorExtractor->compute( img2, keypoints2, descriptors2 );

    cout << "Matching descriptors..." << endl;
    vector<DMatch> filteredMatches;
	crossCheckMatching( descriptorMatcher, descriptors1, descriptors2, filteredMatches );
    vector<int> queryIdxs( filteredMatches.size() ), trainIdxs( filteredMatches.size() );
    for( size_t i = 0; i < filteredMatches.size(); i++ )
    {
        queryIdxs[i] = filteredMatches[i].queryIdx;
        trainIdxs[i] = filteredMatches[i].trainIdx;
    }

	vector<Point2f> points1; KeyPoint::convert(keypoints1, points1, queryIdxs);
    vector<Point2f> points2; KeyPoint::convert(keypoints2, points2, trainIdxs);

	int nPtMatches = points1.size();
	match.clear();
	for(int i = 0;i < nPtMatches;i++){
		MatchPt tmp;
		tmp.P1 = points1[i];
		tmp.P2 = points2[i];
		match.push_back(tmp);
	}

	return nPtMatches;
}

int main( )
{
    char* argv[3];
	argv[1]="E:\\1_A.jpg";
	argv[2]="E:\\1_B.jpg";

	//files of lines
	string dname1="E:\\1Aline.txt";
	string dname2="E:\\1Bline.txt";
	//file of matched points
	string fname="E:\\1ABpoint.txt"; 


	ofstream pout(fname,ios_base::trunc);
	ofstream lout1(dname1,ios_base::trunc);
	ofstream lout2(dname2,ios_base::trunc);


    int nLine1, nLine2;
	Line* pLine1 = LineDetect(argv[1], nLine1);

	Line* pLine2 = LineDetect(argv[2], nLine2);

	cout << "extracted lines : " << nLine1 << " , " << nLine2 << endl;

	CvFont font;
	double hScale=1;
	double vScale=1;
	cvInitFont(&font,CV_FONT_HERSHEY_PLAIN, hScale,vScale,0,1);
	IplImage* res_ima = cvLoadImage(argv[1]);
	for (int i = 0;i < nLine1;i++)
	{
		CvPoint start_pt = cvPoint((int)pLine1[i].StartPt.x,(int)pLine1[i].StartPt.y);
		CvPoint end_pt = cvPoint((int)pLine1[i].EndPt.x,(int)pLine1[i].EndPt.y);
		////////////////////////
	     lout1<<start_pt.x<<" "<<start_pt.y<<" "<<end_pt.x<<" "<<end_pt.y<<endl;
		/////////////////////
		cvLine(res_ima,start_pt,end_pt,CV_RGB(255*((i%3)==0),255*((i%3)==1),255*((i%3)==2)),1,CV_AA);
		CvPoint mid_pt;
		mid_pt.x = (start_pt.x + end_pt.x) / 2;
		mid_pt.y = (start_pt.y + end_pt.y) / 2;
		char text[100];
		sprintf(text,"%d", i+1);
		cvPutText(res_ima,text,mid_pt,&font,CV_RGB(255*((i%3)==0),255*((i%3)==1),255*((i%3)==2)));
	}
	lout1.close();/////////////////////

	IplImage* res_imb = cvLoadImage(argv[2]);
	for (int i = 0;i < nLine2;i++)
	{
		CvPoint start_pt = cvPoint((int)pLine2[i].StartPt.x,(int)pLine2[i].StartPt.y);
		CvPoint end_pt = cvPoint((int)pLine2[i].EndPt.x,(int)pLine2[i].EndPt.y);
		////////////////////////
		lout2<<start_pt.x<<" "<<start_pt.y<<" "<<end_pt.x<<" "<<end_pt.y<<endl;
		/////////////////////
		cvLine(res_imb,start_pt,end_pt,CV_RGB(255*((i%3)==0),255*((i%3)==1),255*((i%3)==2)),1,CV_AA);
		CvPoint mid_pt;
		mid_pt.x = (start_pt.x + end_pt.x) / 2;
		mid_pt.y = (start_pt.y + end_pt.y) / 2;
		char text[100];
		sprintf(text,"%d", i+1);
		cvPutText(res_imb,text,mid_pt,&font,CV_RGB(255*((i%3)==0),255*((i%3)==1),255*((i%3)==2)));
	}

    lout2.close();//////////////////
	cvNamedWindow("a");
	cvNamedWindow("b");
	cvShowImage("a",res_ima);
	cvShowImage("b",res_imb);

	cvWaitKey(0);
	cvReleaseImage(&res_ima);
	cvReleaseImage(&res_imb);
	cvDestroyWindow("a");
	cvDestroyWindow("b");

	vector<MatchPt> PtMatches;
	int nPtMatch = PtMatch(argv[1],argv[2],PtMatches );
	cout << nPtMatch << " matched points" << endl;

	for(int i=0;i<PtMatches.size();i++)
	{
		///////////////////////////
		pout<<PtMatches[i].P1.x<<" "<<PtMatches[i].P1.y<<" "<<PtMatches[i].P2.x<<" "<<PtMatches[i].P2.y<<endl;
		//////////////////////////
	}
	pout.close();

	
}