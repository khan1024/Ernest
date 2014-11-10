factorization_method=0; % default Li's factorization

if grating==0
    cd('structures');
    copyfile(open_grating_file,'temporary_f.m')
    temporary_f;
    delete('temporary_f.m');
    cd('..');
end

if grating==0 && optimalization_grating_0==1

else
    %{
    homogenous_layer=0.*thickness;
    for i=1:1:number_of_layers   
        temp_px=eval(['coordinate_x_' num2str(i)]);
        eval(['new_px_' num2str(i), '=temp_px;']);
    
        temp_py=eval(['coordinate_y_' num2str(i)]);
        eval(['new_py_' num2str(i), '=temp_py;']);
    
        temp_n_mat=eval(['r_index_' num2str(i)]);
        eval(['new_n_mat_' num2str(i), '=temp_n_mat;']);
        
        if abs((max(max(temp_n_mat)))-min(min(temp_n_mat)))<1E-12 % je homogenni        
            homogenous_layer(i)=1; % je homogenni        
            current_hom_n(i)=temp_n_mat(1,1);    
        else            
            homogenous_layer(i)=2; % neni homogenni    
        end
        
    end
    
    exist_nvm=exist('nvm_method','var');
    if exist_nvm==1
        if length(nvm_method)~=number_of_layers
            nvm_method_(1:number_of_layers)=0;
            nvm_method_(1:length(nvm_method))=nvm_method_(1:length(nvm_method))+nvm_method;    
        else
            nvm_method_=nvm_method;
        end
    else        
        nvm_method_(1:number_of_layers)=0;
    end
    %}
    
end

% zatím se budou počítat N_xx, ... pro všechny kroky
if factorization_method~=0
    if measurement>=2 && measurement<=6
        exist_loop=exist('loop_index','var');
        
        if exist_loop==0
            mex_speedup=1;        
            use_grid_asr_speed_up=1;
            clear N_xx N_yy N_xy
            if use_fftw==2            
                if factorization_method==1 || factorization_method==2                
            
                    N_xx=asr_fourier_2D(p_x,p_y,p_x,p_y,N_x_input.*N_x_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);            
                    N_yy=asr_fourier_2D(p_x,p_y,p_x,p_y,N_y_input.*N_y_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);            
                    N_xy=asr_fourier_2D(p_x,p_y,p_x,p_y,N_x_input.*N_y_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);            
            
                elseif factorization_method==3 || factorization_method==4                
                
                    N_xx=asr_fourier_2D(p_x,p_y,p_x,p_y,N_x_input.*conj(N_x_input),N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);                
                    N_yy=asr_fourier_2D(p_x,p_y,p_x,p_y,N_y_input.*conj(N_y_input),N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);                
                    N_xy=asr_fourier_2D(p_x,p_y,p_x,p_y,N_x_input.*conj(N_y_input),N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);                
                    N_yx=asr_fourier_2D(p_x,p_y,p_x,p_y,N_y_input.*conj(N_x_input),N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);                
            
                end
                
            elseif use_fftw==1            
                if factorization_method==1 || factorization_method==2            
                
                    N_xx=fft2D(N_x_input.*N_x_input,[N_X,N_Y,n_sq]);
                    N_yy=fft2D(N_y_input.*N_y_input,[N_X,N_Y,n_sq]);
                    N_xy=fft2D(N_x_input.*N_y_input,[N_X,N_Y,n_sq]);
            
                elseif factorization_method==3 || factorization_method==4
                
                    N_xx=fft2D(N_x_input.*conj(N_x_input),[N_X,N_Y,n_sq]);
                    N_yy=fft2D(N_y_input.*conj(N_y_input),[N_X,N_Y,n_sq]);
                    N_xy=fft2D(N_x_input.*conj(N_y_input),[N_X,N_Y,n_sq]);
                    N_yx=fft2D(N_y_input.*conj(N_x_input),[N_X,N_Y,n_sq]);
            
                end
                
            end
            loop_index=1; % zatím se budou počítat N_xx, ... pro všechny kroky
            
        elseif exist_loop==1
            loop_index=1;
        end
    else
        n_x=2*N_X+1;
        n_y=2*N_Y+1;
        
        n_sq=n_x*n_y;
        
        mex_speedup=1;    
        use_grid_asr_speed_up=1;    
        clear N_xx N_yy N_xy    
        if use_fftw==2
            if factorization_method==1 || factorization_method==2
                
            N_xx=asr_fourier_2D(p_x,p_y,p_x,p_y,N_x_input.*N_x_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
            N_yy=asr_fourier_2D(p_x,p_y,p_x,p_y,N_y_input.*N_y_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);    
            N_xy=asr_fourier_2D(p_x,p_y,p_x,p_y,N_x_input.*N_y_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
            
            elseif factorization_method==3 || factorization_method==4
                
                N_xx=asr_fourier_2D(p_x,p_y,p_x,p_y,N_x_input.*conj(N_x_input),N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
                N_yy=asr_fourier_2D(p_x,p_y,p_x,p_y,N_y_input.*conj(N_y_input),N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);    
                N_xy=asr_fourier_2D(p_x,p_y,p_x,p_y,N_x_input.*conj(N_y_input),N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
                N_yx=asr_fourier_2D(p_x,p_y,p_x,p_y,N_y_input.*conj(N_x_input),N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
                
            end
        elseif use_fftw==1
            if factorization_method==1 || factorization_method==2
            
                N_xx=fft2D(N_x_input.*N_x_input,[N_X,N_Y,n_sq]);            
                N_yy=fft2D(N_y_input.*N_y_input,[N_X,N_Y,n_sq]);            
                N_xy=fft2D(N_x_input.*N_y_input,[N_X,N_Y,n_sq]);
        
            elseif factorization_method==3 || factorization_method==4
            
                N_xx=fft2D(N_x_input.*conj(N_x_input),[N_X,N_Y,n_sq]);            
                N_yy=fft2D(N_y_input.*conj(N_y_input),[N_X,N_Y,n_sq]);            
                N_xy=fft2D(N_x_input.*conj(N_y_input),[N_X,N_Y,n_sq]);
                N_yx=fft2D(N_y_input.*conj(N_x_input),[N_X,N_Y,n_sq]);
            end
        end
    end
      
end