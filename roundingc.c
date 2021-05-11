#include "mex.h"
#include <stdio.h>
#include <math.h>


void mexFunction (
  int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
  
{

	double *x, *f, *w;
	int rows_x, columns_x; 


	if (nrhs!=1) {
		mexErrMsgTxt("Sintax: rounding(x)\n");
	  }

	/******************************************************************/
	  /* The dimensions of vector x are checked */
	  x = mxGetPr(prhs[0]);
	  rows_x = mxGetM(prhs[0]);
	  columns_x = mxGetN(prhs[0]);

	  /* printf("Dimensions of x: %d x %d\n", rows_x, columns_x);
	  printf("Input: %f\n", *x); */

	  plhs[0] = mxCreateDoubleMatrix(rows_x, columns_x, mxREAL);
	  plhs[1] = mxCreateDoubleMatrix(rows_x, columns_x, mxREAL);
	  f = mxGetPr(plhs[0]);   /* x rounded*/
	  w = mxGetPr(plhs[1]);   /* x rounded "the wrong way"*/

	  if (*x==0) {
		  *f = 0;
		  *w = 1;
	  }
	  else if (*x>0) {
		if ((*x - floor(*x))<=0.5) {
			*f = floor(*x);
			*w = *f + 1;
		}
		else {
			*f = ceil(*x);
			*w = *f - 1;
		}
	  }
	  else {
		if ((ceil(*x) - *x)<=0.5) {
			*f = ceil(*x);
			*w = *f - 1;
		}
		else {
			*f = floor(*x);
			*w = *f + 1;
		}
	  }
}
