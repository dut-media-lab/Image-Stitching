#include "mex.h"
#include <math.h>

/* give a point and a quadrangle, find whether it is inside the quadrangle*/
bool JudgeInsideQuad( double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4, double px, double py)
{
    double Denominator1=0, Denominator2=0;
    double Numerator1=0, Numerator2=0, Numerator3=0;
    double Numerator4=0, Numerator5=0, Numerator6=0;
    bool InsideFlag1 = false;
    bool InsideFlag2 = false;

    Denominator1 = (px2-px1)*(py3-py1)-(px3-px1)*(py2-py1);
    Numerator1   = (px2-px)*(py3-py)-(px3-px)*(py2-py);
    Numerator2   = (px-px1)*(py3-py1)-(px3-px1)*(py-py1);
    Numerator3   = (px2-px1)*(py-py1)-(px-px1)*(py2-py1);
    
    if ((Denominator1>0 && Numerator1>0 && Numerator2>0 && Numerator3>0) || (Denominator1<0 && Numerator1<0 && Numerator2<0 && Numerator3<0))
    {
        InsideFlag1 = true; 
    }
    else
    {
        Denominator2 = (px3-px1)*(py4-py1)-(px4-px1)*(py3-py1);
        Numerator4   = (px3-px)*(py4-py)-(px4-px)*(py3-py);
        Numerator5   = (px-px1)*(py4-py1)-(px4-px1)*(py-py1);
        Numerator6   = (px3-px1)*(py-py1)-(px-px1)*(py3-py1);
        
        if ((Denominator2>0 && Numerator4>0 && Numerator5>0 && Numerator6>0) || (Denominator2<0 && Numerator4<0 && Numerator5<0 && Numerator6<0))
        {
            InsideFlag2 = true; 
        }
    }

    return InsideFlag1 || InsideFlag2; 

}

/* give a point and a quadrangle, find the bilinear interpolation coefficients alpha, beta*/
/*% Cal_AlphaInQuad 此函数计算四边形内部的一个点如何被四边形的四个顶点双线�?表示出来，假设此点在四边形中（内部包括边界）
% 具体公式：P=(1-alpha)*(1-beta)*P1+alpha*(1-beta)*P2+alpha*beta*P3+(1-alpha)*beta*P4
% （四边形顶点顺序为P1与P3相对，P2与P1纵坐标接近，处在相近的水平线上）*/
void Cal_AlphaBetaInQuad(double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4, double px, double py, double bi_coeff[2])
{
	/*version-1  by wangchao  */
    // double AlphaCoeff0=0, AlphaCoeff1=0, AlphaCoeff2=0;
    // double Alphadelta=0, Alpha1=0, Alpha2=0;
	
	// double BetaCoeff0=0, BetaCoeff1=0, BetaCoeff2=0;
    // double Betadelta=0, Beta1=0, Beta2=0;
	
	// bi_coeff[0]=0;  bi_coeff[1]=0;

    // AlphaCoeff2 = (px2-px1)*(py3-py4)-(py2-py1)*(px3-px4);
    // AlphaCoeff1 = (px2-px1)*(py4-py)-(py2-py1)*(px4-px)+(px1-px)*(py3-py4)-(py1-py)*(px3-px4);
    // AlphaCoeff0 = (px1-px)*(py4-py)-(py1-py)*(px4-px);

    // if (abs(AlphaCoeff2)<=1e-15) /*此时向量P1P2与向量P3P4平行*/
    // {
        // bi_coeff[0] = -AlphaCoeff0/AlphaCoeff1;
    // }
    // else
    // {
        // Alphadelta = sqrt(AlphaCoeff1*AlphaCoeff1 - 4*AlphaCoeff0*AlphaCoeff2);
        // Alpha1 = (-AlphaCoeff1+Alphadelta)/(2*AlphaCoeff2);
        // Alpha2 = (-AlphaCoeff1-Alphadelta)/(2*AlphaCoeff2);
        // if (Alpha1>=0 && Alpha1<=1)
        // {
            // bi_coeff[0] = Alpha1;
        // }
        // else
        // {
            // bi_coeff[0] = Alpha2;
        // }    
    // }
	
	// BetaCoeff2 = (px4-px1)*(py3-py2) - (py4-py1)*(px3-px2);
    // BetaCoeff1 = (px4-px1)*(py2-py) - (py4-py1)*(px2-px) + (px1-px)*(py3-py2) - (py1-py)*(px3-px2);
    // BetaCoeff0 = (px1-px)*(py2-py) - (py1-py)*(px2-px);
	
	// if (abs(BetaCoeff2)<=1e-15) /*此时向量P1P4与向量P2P3平行*/
    // {
        // bi_coeff[1] = -BetaCoeff0/BetaCoeff1;
    // }
    // else
    // {
        // Betadelta = sqrt(BetaCoeff1*BetaCoeff1-4*BetaCoeff0*BetaCoeff2);
        // Beta1 = (-BetaCoeff1+Betadelta)/(2*BetaCoeff2);
        // Beta2 = (-BetaCoeff1-Betadelta)/(2*BetaCoeff2);
        // if (Beta1>=0 && Beta1<=1)
        // {
           // bi_coeff[1] = Beta1;
        // }
        // else
        // {
           // bi_coeff[1] = Beta2;
        // }    
    // }
	
	/* version-2  by liaotianli using maple 18*/ 
	double Alpha=0, Beta=0;
	double BetaCoeff_a=0, BetaCoeff_b=0, BetaCoeff_c=0;
	double Betadelta = 0, Beta1=0, Beta2=0;
	
	BetaCoeff_a = (px1-px4)*(py2-py3)-(py1-py4)*(px2-px3);
	BetaCoeff_b = (-py1+py2-py3+py4)*px+(px1-px2+px3-px4)*py-2*px1*py2+px1*py3+2*px2*py1-px2*py4-px3*py1+px4*py2;
	BetaCoeff_c = (py1-py2)*px+(-px1+px2)*py+px1*py2-px2*py1;
		
	if (fabs(BetaCoeff_a)<=1e-15)
	{	
		if (fabs(BetaCoeff_b)>1e-15)  //  to avoid system crush
		{
			Beta = -BetaCoeff_c/BetaCoeff_b;
		}
	}
	else
	{
        Betadelta = BetaCoeff_b*BetaCoeff_b-4*BetaCoeff_a*BetaCoeff_c;
        if (Betadelta>=0)  //  to avoid system crush
        {
            Beta1 = (-BetaCoeff_b + sqrt(Betadelta))/(2*BetaCoeff_a);
            Beta2 = (-BetaCoeff_b - sqrt(Betadelta))/(2*BetaCoeff_a);
            if (Beta1>=0 && Beta1<=1)
            {
                Beta = Beta1;
            }
            else
            {
                Beta = Beta2;
            }
        }
	}
    
    if (fabs(Beta*px1-Beta*px2+Beta*px3-Beta*px4-px1+px2)>1e-15)  //  to avoid system crush
    {
        Alpha = (Beta*px1-Beta*px4+px-px1)/(Beta*px1-Beta*px2+Beta*px3-Beta*px4-px1+px2);
    }
	   
    if (Alpha>=0 && Alpha<=1 && Beta>=0 && Beta<=1)
    {
        bi_coeff[0] = Alpha;
        bi_coeff[1] = Beta;
    }
    else
    {
        bi_coeff[0]=0; bi_coeff[1]=0;  // to avoid system crush
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Input/output variables. */
    double *img; /* original input image. */
    double dch;  /*  height of warped image */
    double dcw;  /*  width of warped image */
    double dC1;  /* number of mesh rows */ 
    double dC2;  /* number of mesh columns*/
    double *X;      /* X positions of the original mesh grid. */
    double *Y;      /* Y positions of the original mesh grid. */ 
    double *wX;     /* X positions of the warped mesh grid. */
    double *wY;     /* Y positions of the warped mesh grid. */   
    double *off;   /* offset of warped image r.t. original image */ 
              
    double *warped_img;     /* Warped img (img after being warped with mesh transformation */
                              
    /* Intermediate variables.*/
    int ch, cw, C1, C2;  /* int-type ch, cw, C1, C2*/
    int imgm, imgn;      /* 'm' and 'n': size of the original image. */    
    int x, y;            /* x,y positions in each warped grid*/
    int mleft, mright, mup, mdown;   /* the external rectangle of each warped grid */

	int indij;          /* 1D-indice of X,Y (or wX, wY)*/
    int i, j;          

    bool InsideFlag;        /* judge whether a point is inside quadrangle */
	double alpha, beta;   /* alpha, beta coefficients of bilinear representation */
	double bi_coeff[2]={0};
    double posx, posy;
    int px, py, cidx, sidx;        /*contain the indexes of the warped and original images respectively */
    
    /* Since MATLAB works with 2D matrices and C works with 1D arrays we need the following variables for generating the 
     * final stitched RGB 2D matrix for MATLAB, this image/matrix contains the warped img2. */
    int ch1canv; /* Displacement for channel 1 in image 1 (R). */
    int ch2canv; /* Displacement for channel 2 in image 1 (G). */
    int ch3canv; /* Displacement for channel 3 in image 1 (B). */
    
    int ch1img; /* Displacement for channel 1 in image 2 (R). */
    int ch2img; /* Displacement for channel 2 in image 2 (G). */
    int ch3img; /* Displacement for channel 3 in image 2 (B). */
    
    
    /* Check for proper number of arguments. */    
    if (nrhs!=10)
    {
        mexErrMsgTxt("Wrong number of inputs.");
    }
    else if (nlhs!=1)
    {
        mexErrMsgTxt("Wrong number of output arguments.");
    }
	if (!mxIsDouble(prhs[0]))
	{
		mexErrMsgTxt("Input must be type of double.");
	}
	
    /* Assign scalars to inputs. */
    dch  = mxGetScalar(prhs[1]);
    dcw  = mxGetScalar(prhs[2]);
    dC1  = mxGetScalar(prhs[3]);
    dC2  = mxGetScalar(prhs[4]);
    
    ch = (int)dch;  /* height of warped image */
    cw = (int)dcw;  /* width of warped image */
    C1 = (int)dC1;  /* rows of mesh grid */
    C2 = (int)dC2;  /* columns of mesh grid */
    
    /* Assign pointers to inputs. */   
    img      = mxGetPr(prhs[0]);
    X        = mxGetPr(prhs[5]);
    Y        = mxGetPr(prhs[6]);
    wX       = mxGetPr(prhs[7]);
    wY       = mxGetPr(prhs[8]); 
    off      = mxGetPr(prhs[9]);  

    imgm = mxGetM(prhs[0]);
    imgn = mxGetN(prhs[0])/3; /* It is an RGB image, which means that for MATLAB it is an Mx(Nx3) image. So, we have  
                                  * to divide N by 3 in order to get the proper size of the image in C. */
    
    /* Create matrix for the return arguments. */
    plhs[0] = mxCreateDoubleMatrix( ch, cw*3, mxREAL); 

    /* Assign pointers to output canvas (warped image2). */
    warped_img  = mxGetPr(plhs[0]);       
    
    /* Initialize displacements. */
    ch1canv = 0;
    ch2canv = ch*cw;
    ch3canv = ch*cw*2;
    
    ch1img = 0;
    ch2img = imgn*imgm;
    ch3img = imgn*imgm*2;
    
    /*initialization*/
    for (j=0; j<cw;j++)
    {
        for (i=0;i<ch;i++)
        {
            warped_img[j*ch+i+ch1canv]=0;
            warped_img[j*ch+i+ch2canv]=0;
            warped_img[j*ch+i+ch3canv]=0;
        }
    }
	
    
    /* Start computations. */
    /* For each grid. */
    for(j=0;j<C2;j++)
    {
        for(i=0;i<C1;i++)
        {            
			indij = j*(C1+1)+i;
            mleft  = (int)floor( fmin(fmin(wX[indij], wX[indij+1]), fmin(wX[indij+C1+1], wX[indij+C1+2])) );
            mright = (int)ceil( fmax(fmax(wX[indij], wX[indij+1]), fmax(wX[indij+C1+1], wX[indij+C1+2])) );
            mup    = (int)floor( fmin(fmin(wY[indij], wY[indij+C1+1]), fmin(wY[indij+1], wY[indij+C1+2])) );
            mdown  = (int)ceil( fmax(fmax(wY[indij], wY[indij+C1+1]), fmax(wY[indij+1], wY[indij+C1+2])) );     
            /* Get grid pixel for current grid. */
            for(x=mleft; x <= mright; x++)
            {
                for(y=mup; y <= mdown ; y++)
                {
                    /* use bilinear method to interpolate  */
                    InsideFlag = JudgeInsideQuad( wX[indij], wY[indij], wX[indij+C1+1], wY[indij+C1+1], wX[indij+C1+2], wY[indij+C1+2], wX[indij+1], wY[indij+1], x, y );
                    if (InsideFlag)
                    {
						Cal_AlphaBetaInQuad(wX[indij], wY[indij], wX[indij+C1+1], wY[indij+C1+1], wX[indij+C1+2], wY[indij+C1+2], wX[indij+1], wY[indij+1], x, y, bi_coeff);
						alpha = bi_coeff[0];
						beta  = bi_coeff[1];
						posx = (1-alpha)*(1-beta)*X[indij] + alpha*(1-beta)*X[indij+C1+1] + alpha*beta*X[indij+C1+2] + (1-alpha)*beta*X[indij+1];
                        posy = (1-alpha)*(1-beta)*Y[indij] + alpha*(1-beta)*Y[indij+C1+1] + alpha*beta*Y[indij+C1+2] + (1-alpha)*beta*Y[indij+1];
                        px = (int)floor(posx);
                        py = (int)floor(posy);
                        cidx = (int)( (x+off[0]-1)*ch + y+off[1]-1 );
                        sidx = (int)(px-1)*imgm + py-1; 
                        /* texture mapping via bilinear interpolation */
                        warped_img[cidx+ch1canv] = (px+1-posx)*(py+1-posy)*img[sidx+ch1img] + (px+1-posx)*(posy-py)*img[sidx+1+ch1img] + (posx-px)*(py+1-posy)*img[sidx+imgm+ch1img] + (posx-px)*(posy-py)*img[sidx+imgm+1+ch1img];
                        warped_img[cidx+ch2canv] = (px+1-posx)*(py+1-posy)*img[sidx+ch2img] + (px+1-posx)*(posy-py)*img[sidx+1+ch2img] + (posx-px)*(py+1-posy)*img[sidx+imgm+ch2img] + (posx-px)*(posy-py)*img[sidx+imgm+1+ch2img];
                        warped_img[cidx+ch3canv] = (px+1-posx)*(py+1-posy)*img[sidx+ch3img] + (px+1-posx)*(posy-py)*img[sidx+1+ch3img] + (posx-px)*(py+1-posy)*img[sidx+imgm+ch3img] + (posx-px)*(posy-py)*img[sidx+imgm+1+ch3img];
                    }
                }
            }
        }
    }
            
    /* End.*/
    return;
}
