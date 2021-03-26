/*-------------------------------------------------------------------------
* Robust Multi-Structure Fitting and Model Selection via Preference
* Correlation
* -------------------------------------------------------------------------
* The demo code in this package implements the method for multi-structure
* robust fitting and model selection proposed in:
* T.-J. Chin, H. Wang and D. Suter
* Robust Multi-Structure Fitting and Model Selection via Preference
* Correlation
* Technical Report, ACVT, Adelaide, South Australia, 2011.
*
* Copyright (c) 2011 Tat-Jun Chin
* School of Computer Science, The University of Adelaide, South Australia
* http://www.cs.adelaide.edu.au/~tjchin
*
* An older version of the method is described in
* T.-J. Chin, H. Wang and D. Suter
* Robust Fitting of Multiple Structures: The Statistical Learning Approach
* In Proc. Int. Conf. on Computer Vision 2009, Kyoto, Japan
*
* The program is free for non-commercial academic use. Any commercial use
* is strictly prohibited without the authors' consent. Please acknowledge
* the authors by citing the above paper in any academic publications that
* have made use of this package or part of it.
*
* If you encounter any problems or questions please send email to
* tjchin@cs.adelaide.edu.au.*/
        
#include "mex.h"
#include "math.h"
#include "time.h"

double intersect(double *, double *, int, double);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Input/output variables.*/
    double *K;
    double *A1;
    double *A2;
    double *H;
  
    /* Intermediate variables.*/
    int m1,n1,m2,n2;
    int i,j;
  
    /* Check for proper number of arguments. */    
    if (nrhs != 3)
    {
        mexErrMsgTxt("Three inputs required.");
    }
    else if (nlhs > 1)
    {
        mexErrMsgTxt("Too many output arguments.");
    }
  
    /* Assign pointers to inputs.*/
    A1 = mxGetPr(prhs[0]);
    A2 = mxGetPr(prhs[1]);
    H = mxGetPr(prhs[2]);
  
    /* Get sizes of input matrices.*/
    m1 = mxGetM(prhs[0]);
    n1 = mxGetN(prhs[0]);
    m2 = mxGetM(prhs[1]);
    n2 = mxGetN(prhs[1]);
  
    /* Create matrix for the return argument. */
    plhs[0] = mxCreateDoubleMatrix(n1,n2,mxREAL); 
  
    /* Assign pointers to output.*/
    K = mxGetPr(plhs[0]);  
  
    /* Start computations.*/
    for(j=0;j<n1;j++)
    {
        for(i=0;i<n2;i++)
        {
            K[ j + i*n1 ] = intersect( A1+j*m1, A2+i*m2, m1, *H);         
        }
    }
    
    /* Bye bye.*/
    return;
}

double intersect(double *x, double *z, int pix, double param)
{   
    /* Intersection Kernel- new faster version.*/
    int *symtab;
    double accum;
    int n;
    
    /* Initialize symbol table.*/
    symtab = (int *)mxCalloc(pix,sizeof(int));

    /* Fill the tables.*/
    for(n=0;n<param;n++)
    {   
        symtab[(int)x[n]-1] = 1;        
    }

    /* Intersect.*/
    accum = 0;
    for(n=0;n<param;n++)
    {
        if (symtab[(int)z[n]-1]==1)
        {
            accum++;
        }
    }

    /* Housekeeping.*/
    mxFree(symtab);

    /* Normalized intersection.*/    
    return(accum/param); 
}
