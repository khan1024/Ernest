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
    int n_x;
    double *alpha, *beta;
    double *alpha_m, *beta_n;
    int m,n,j,l;
    long index;
    int r_result,s_result;

    n_x=(int)(mxGetScalar(prhs[0]));
    
    alpha_m = mxGetPr(prhs[1]);
    beta_n = mxGetPr(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(n_x*n_x, n_x*n_x, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n_x*n_x, n_x*n_x, mxREAL);

    alpha=mxGetPr(plhs[0]);
    beta=mxGetPr(plhs[1]);    
    
    for(m=0;m<n_x;m++){        
      for(n=0;n<n_x;n++){            
	for(j=0;j<n_x;j++){                
	  for(l=0;l<n_x;l++){                    
		  
	    if(m==j && n==l){		    
		    
	      r_result=m*n_x+n+1;                    
	      s_result=j*n_x+l+1;                    
                    
	      index=((s_result-1)*n_x*n_x)+(r_result-1);                    
                    
	      alpha[index]=alpha_m[m];                    
	      beta[index]=beta_n[n];		  
	      
	    }
	  }
	}
      }
    }
    
}
