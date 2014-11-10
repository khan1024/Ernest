% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

function result=asr_fourier_2D_mex(n_x,results)

result=zeros(n_x*n_x,n_x*n_x);

for m=0:1:(n_x-1)
    for n=0:1:(n_x-1) 
        for j=0:1:(n_x-1)
            for l=0:1:(n_x-1)
                result(m*n_x+n+1,j*n_x+l+1)=results(n_x+m-j,n_x+n-l);
            end
        end
    end
end