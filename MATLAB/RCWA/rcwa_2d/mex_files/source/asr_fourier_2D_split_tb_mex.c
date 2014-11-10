#include "mex.h"
#include "math.h"

/*
Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
Czech Technical University in Prage, Optical Physics Group, Czech Republic
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int r,s;
    double *result_r, *result_i;
    double *results_r, *results_i;
    int m,n,j,l;
    long index_result;
    
    int M,N;
    long MN;
    long i1, i2;
    
    int rr_inv_toep;
    
    double *inv_f_toep_r, *inv_f_toep_i, *fourier_1D_r, *fourier_1D_i;
    
    int r_inv_toep, s_inv_toep;
    long i_t,i_f;
    double real_part, imag_part;
    
    M=(int)(mxGetScalar(prhs[0]));
        
    result_r = mxGetPr(prhs[1]);
    result_i = mxGetPi(prhs[1]);
    r=mxGetM(prhs[1]);
    s=mxGetN(prhs[1]); 
    
    inv_f_toep_r = mxGetPr(prhs[2]);
    inv_f_toep_i = mxGetPi(prhs[2]);
    rr_inv_toep=mxGetM(prhs[2]);
    
    fourier_1D_r = mxGetPr(prhs[3]);
    fourier_1D_i = mxGetPi(prhs[3]);
    
    N=M;
    MN=(2*M+1)*(2*M+1);
    plhs[0] = mxCreateDoubleMatrix(MN,MN, mxCOMPLEX);

    results_r=mxGetPr(plhs[0]);
    results_i=mxGetPi(plhs[0]);
    
    for(m=-M;m<M+1;m++){
        for(n=-N;n<N+1;n++){
            i1=(n+N+1)+(m+M)*(2*N+1);
            for(j=-M;j<M+1;j++){
                for(l=-N;l<N+1;l++){
                    
                    i2=(l+N+1)+(j+M)*(2*N+1);        
                    
                    index_result=((i2-1)*MN)+(i1-1);
                    
                    r_inv_toep=n+N+1;
                    s_inv_toep=l+N+1;
                    
                    i_f=m-j+2*M;
                    i_t=((s_inv_toep-1)*rr_inv_toep)+(r_inv_toep-1);
                    
                    real_part=inv_f_toep_r[i_t]*fourier_1D_r[i_f]-inv_f_toep_i[i_t]*fourier_1D_i[i_f];
                    imag_part=inv_f_toep_r[i_t]*fourier_1D_i[i_f]+inv_f_toep_i[i_t]*fourier_1D_r[i_f];
                    
                    results_r[index_result]=result_r[index_result]+real_part;
                    results_i[index_result]=result_i[index_result]+imag_part;
                                              
                }
            }
        }
    }
    
    
}
