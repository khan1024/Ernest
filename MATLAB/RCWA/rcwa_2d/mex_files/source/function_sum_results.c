#include "mex.h"
#include "math.h"

/*
Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
Czech Technical University in Prage, Optical Physics Group, Czech Republic
*/

void DisplayMatrix(char *Name, double *Data, int M, int N)
{ /* Display matrix data */
  int m, n;
  mexPrintf("%s = \n", Name);
  for(m = 1; m < (M+1); m++, mexPrintf("\n"))
    for(n = 1; n < (N+1); n++)
      mexPrintf("%8.4f ", Data[(m-1) + M*(n-1)]);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int r,s,n_x;
    double *results_r, *results_i, *result_r, *result_i;
    
    double *fourier_1d_u_r, *fourier_1d_u_i, *fourier_1d_v_r, *fourier_1d_v_i;
    int m,n,M,N;
    long index_s, index, i_u,i_v;
    int r_result,s_result;
    
    r=mxGetM(prhs[1]);
    s=mxGetN(prhs[1]);

    M=(int)(mxGetScalar(prhs[0]));
    n_x=2*M+1;
    N=M;
    
    result_r = mxGetPr(prhs[1]);
    result_i = mxGetPi(prhs[1]);
    
    fourier_1d_u_r = mxGetPr(prhs[2]);
    fourier_1d_u_i = mxGetPi(prhs[2]);
    
    fourier_1d_v_r = mxGetPr(prhs[3]);
    fourier_1d_v_i = mxGetPi(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(2*M+1, 2*M+1, mxCOMPLEX);

    results_r=mxGetPr(plhs[0]);
    results_i=mxGetPi(plhs[0]);
    
    for(m=-M;m<(M+1);m++){
        for(n=-N;n<(N+1);n++){
                    
	    r_result=m+M+1;
            s_result=n+N+1;
                    
            index=((s_result-1)*(2*M+1))+(r_result-1);
		    
	    i_u=m+M;
	    i_v=n+N;
	    
            results_r[index]=result_r[index]+fourier_1d_u_r[i_u]*fourier_1d_v_r[i_v]-fourier_1d_u_i[i_u]*fourier_1d_v_i[i_v];
            results_i[index]=result_i[index]+fourier_1d_u_r[i_u]*fourier_1d_v_i[i_v]+fourier_1d_u_i[i_u]*fourier_1d_v_r[i_v];
          
        }
    }
    
}
