% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

function result=asr_fourier_2D(disc_x,disc_y,disc_u,disc_v,f,M_X,N_X,G,mex_speedup,use_grid_asr_speed_up)

% based on Camfr (camfr.sourceforge.net) Fourier procedure

M=2*M_X;
N=2*N_X;

d_u=disc_u(end)-disc_u(1);
d_v=disc_v(end)-disc_v(1);

results=zeros(2*M+1,2*N+1)+1e-100*1i;

for i=1:1:length(disc_u)-1
    
    use_grid_asr_speed_up=1;
    if use_grid_asr_speed_up==1
        if i==1
            asr_function_only_v=asr_fourier_only_asr(disc_y,disc_v,N,d_v,G);
            fourier_1D_v=zeros(1,2*N+1);       
            fourier_1D_v=asr_fourier(f(:,i),disc_y,disc_v,N,d_v,G,use_grid_asr_speed_up,asr_function_only_v);
        else
            fourier_1D_v=zeros(1,2*N+1);      
            fourier_1D_v=asr_fourier(f(:,i),disc_y,disc_v,N,d_v,G,use_grid_asr_speed_up,asr_function_only_v);
        end
    else
        fourier_1D_v=zeros(1,2*N+1);
        fourier_1D_v=asr_fourier(f(:,i),disc_y,disc_v,N,d_v,G,use_grid_asr_speed_up,0);
    end
    
    disc_i_x=[disc_x(i),disc_x(i+1)];
    disc_i_u=[disc_u(i),disc_u(i+1)];
    f_i_u=1;
    
    use_grid_asr_speed_up=2;
    fourier_1D_u=zeros(1,2*M+1);
    fourier_1D_u=asr_fourier(f_i_u,disc_i_x,disc_i_u,M,d_u,G,use_grid_asr_speed_up,0);
    
    fourier_1D_u(1)=fourier_1D_u(1)+1e-100*1i;
    fourier_1D_v(1)=fourier_1D_v(1)+1e-100*1i;

    if mex_speedup==1
    results=function_sum_results(M,results,fourier_1D_u,fourier_1D_v);
    else

    for m=-M:1:M
        for n=-N:1:N
            results(m+M+1,n+N+1)=results(m+M+1,n+N+1)+fourier_1D_u(m+M+1)*fourier_1D_v(n+N+1);
        end
    end
    end

end

n_x=M+1;
result=zeros(n_x*n_x,n_x*n_x);

if mex_speedup==1
    results=results+1e-100*1i;
    result=asr_fourier_2D_mex(n_x,results);
else
    
    for m=0:1:(n_x-1)    
        for n=0:1:(n_x-1)        
            for j=0:1:(n_x-1)            
                for l=0:1:(n_x-1)                
                    result(m*n_x+n+1,j*n_x+l+1)=results(n_x+m-j,n_x+n-l);            
                end                
            end            
        end        
    end
    
end
