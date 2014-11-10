#include "mex.h"
#include "math.h"

/*
Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
Czech Technical University in Prage, Optical Physics Group, Czech Republic
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double x_input,y_input,x_1,y_1,x_2,y_2;
    int direction;
    double *delete;
    double a,b,y_line;

    x_input=(double)(mxGetScalar(prhs[0]));
    y_input=(double)(mxGetScalar(prhs[1]));
    x_1=(double)(mxGetScalar(prhs[2]));
    y_1=(double)(mxGetScalar(prhs[3]));
    x_2=(double)(mxGetScalar(prhs[4]));
    y_2=(double)(mxGetScalar(prhs[5]));
    
    direction=(int)(mxGetScalar(prhs[6]));
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    delete=mxGetPr(plhs[0]);
    
    delete[0]=1;
    
    if(x_1==x_2){
        if(direction==2){
            if(x_input<x_1){
                delete[0]=2; 
            }
        }
    
        if(direction==1){  
            if(x_input>=x_1){
                delete[0]=2;
            } 
        }
    }   
    else{      
        a=(y_2-y_1)/(x_2-x_1);      
        b=y_1-a*x_1;
      
        y_line=a*x_input+b;
      
        if(direction==1){	
            if(y_input<=y_line){
                delete[0]=2;
            }
        }
        else if(direction==2){
            if(y_input>=y_line){    
                delete[0]=2;
            }
        }
    }
}
