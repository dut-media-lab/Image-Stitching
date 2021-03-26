// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the LINEMATCHER_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// LINEMATCHER_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef LINEMATCHER_EXPORTS
#define LINEMATCHER_API __declspec(dllexport)
#else
#define LINEMATCHER_API __declspec(dllimport)
#endif

#include <vector>
using namespace std;
#include <cxcore.h>

struct Pt2LineInfo 
{
	double dis;          //点到直线的距离
	double dis2;         //点到直线垂线的距离
	double side;         //正、负号各表示一边，0在线上
};

struct Line 
{
	CvPoint2D32f Center;
	float length;
	float para_a;
	float para_b;
	float para_c;
	CvPoint2D32f StartPt;
	CvPoint2D32f EndPt;
};

struct MatchPt
{
	CvPoint2D32f P1;
	CvPoint2D32f P2;
};

struct MatchLine 
{
	int ID1;
	int ID2;
	double dis;
};

/* Line Matching by affine invariants computed from matched points. If 'bFast' = true, then 'angle' is the estimated 
global rotation between two images. If 'bFast' = false, then 'angle' is not used. */
LINEMATCHER_API int _stdcall MatchingLines(Line* pLine1,Line* pLine2,int nLine1,int nLine2,float angle,
											const vector<MatchPt> &match, vector<MatchLine> &outLM,bool bFast = false);