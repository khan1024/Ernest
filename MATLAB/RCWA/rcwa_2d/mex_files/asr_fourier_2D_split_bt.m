% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

function result=asr_fourier_2D_split_bt(disc_x,disc_y,disc_u,disc_v,ff,M,N,G,mex_speedup,use_grid_asr_speed_up)

% based on Camfr (camfr.sourceforge.net) Fourier procedure

f=1./ff; % vypocet 1/eps

d_u=disc_u(end)-disc_u(1);
d_v=disc_v(end)-disc_v(1);

MN=(2*M+1)*(2*N+1);

result=zeros(MN,MN);

for i=1:1:length(disc_v)-1
    
    use_grid_asr_speed_up=1;
    if use_grid_asr_speed_up==1
        if i==1
            asr_function_only_u=asr_fourier_only_asr(disc_x,disc_u,2*N,d_u,G);
            fourier_1D_u=zeros(1,4*N+1);        
            fourier_1D_u=asr_fourier(f(i,:),disc_x,disc_u,2*N,d_u,G,use_grid_asr_speed_up,asr_function_only_u);
        else
            fourier_1D_u=zeros(1,4*N+1);        
            fourier_1D_u=asr_fourier(f(i,:),disc_x,disc_u,2*N,d_u,G,use_grid_asr_speed_up,asr_function_only_u);
        end
    else
        fourier_1D_u=zeros(1,4*N+1);
        fourier_1D_u=asr_fourier(f(i,:),disc_x,disc_u,2*N,d_u,G,use_grid_asr_speed_up,0);
    end
    
    f_toep=zeros(2*N+1,2*N+1);
    
    for i1=-N:1:N
        for i2=-N:1:N
            f_toep(i1+N+1,i2+N+1)=fourier_1D_u(i1-i2+2*N+1);
        end
    end
    
    inv_f_toep=zeros(2*N+1,2*N+1);
    inv_f_toep=inv(f_toep);
    
    disc_i_y=[disc_y(i),disc_y(i+1)];
    disc_i_v=[disc_v(i),disc_v(i+1)];
    f_i_v=1;
    
    use_grid_asr_speed_up=2;
    fourier_1D_v=zeros(1,4*M+1);
    fourier_1D_v=asr_fourier(f_i_v,disc_i_y,disc_i_v,2*M,d_v,G,use_grid_asr_speed_up,0);
    
    if mex_speedup==1
        inv_f_toep=inv_f_toep+1e-100*1i;
        fourier_1D_v=fourier_1D_v+1e-100*1i;
        result=result+1e-100*1i;
        result=asr_fourier_2D_split_bt_mex(M,result,inv_f_toep,fourier_1D_v);
    else        

        for m=-M:1:M        
            for n=-N:1:N            
                i1=(n+N+1)+(m+M)*(2*N+1);            
                for j=-M:1:M                
                    for l=-N:1:N                    
                        i2=(l+N+1)+(j+M)*(2*N+1);                    
                        result(i1,i2)=result(i1,i2)+fourier_1D_v(n-l+2*M+1)*inv_f_toep(m+N+1,j+N+1);                
                    end                    
                end                
            end            
        end

    
    end
end