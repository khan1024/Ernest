% normal field method

% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

mex_speedup=1;
factorization_Li=1; % 1-Li (default), 2-Lalanne
parameter_alpha=0.5; % Lalanne's factorization parameter

if factorization_Li==1
    disc_x=eval(['new_px_' num2str(layer)]);
    disc_y=eval(['new_py_' num2str(layer)]);
    n_mat_temp=eval(['new_n_mat_' num2str(layer)]);
    eps_mat=n_mat_temp.^2;
    
    use_grid_asr_speed_up=1;

    if use_fftw==1    
        E=fft2D(eps_mat,[N_X,N_Y,n_sq]);    
        E_inv=fft2D(1./eps_mat,[N_X,N_Y,n_sq]);
    else        
        E=asr_fourier_2D(disc_x,disc_y,disc_x,disc_y,eps_mat,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up); % based on Camfr (camfr.sourceforge.net) Fourier procedure
        E_inv=asr_fourier_2D(disc_x,disc_y,disc_x,disc_y,1./eps_mat,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
    end
    
    %N_xx=asr_fourier_2D(disc_x,disc_y,disc_x,disc_y,N_x_input.*N_x_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
    %N_yy=asr_fourier_2D(disc_x,disc_y,disc_x,disc_y,N_y_input.*N_y_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
    %N_xy=asr_fourier_2D(disc_x,disc_y,disc_x,disc_y,N_y_input.*N_x_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
end

%inv_E=inv(E);

if factorization_Li==2 || zeta~=0
    disc_x=eval(['new_px_' num2str(layer)]);
    disc_y=eval(['new_py_' num2str(layer)]);
    n_mat_temp=eval(['new_n_mat_' num2str(layer)]);
    eps_mat=n_mat_temp.^2;
    
    inv_eps_mat=1./eps_mat;
    E_inv=asr_fourier_2D(disc_x,disc_y,disc_x,disc_y,inv_eps_mat,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
    inv_E_inv=inv(E_inv); % Lalanne's paper A=E_inv
end

inv_E=inv(E);

if zeta==0
    if factorization_Li==1
        % equation 34 [1]
        delta_epsilon=E-inv(E_inv);
        delta_N_xy=delta_epsilon*N_xy;
        
        %G=[-alpha*beta, alpha*alpha-(mu*k_0^2).*tb_E;
        %    (mu*k_0^2).*bt_E-beta*beta, alpha*beta];
        
        G=[-alpha*beta+(mu*k_0^2).*delta_N_xy, alpha*alpha+(mu*k_0^2).*(-E+delta_epsilon*N_yy);
            (mu*k_0^2).*(E-delta_epsilon*N_xx)-beta*beta, alpha*beta-(mu*k_0^2).*delta_N_xy];
        
        % equation 33 [1]
        F=[alpha*inv_E*beta, I_mat.*(mu*k_0^2)-alpha*inv_E*alpha;
            beta*inv_E*beta-I_mat.*(mu*k_0^2), -beta*inv_E*alpha];

    elseif factorization_Li==2 
        % equation 1b, Lalanne, JOSA A, vol. 14, p. 1592, 1997
        
        G_1=parameter_alpha*inv_E_inv+(1-parameter_alpha)*E;        
        G_2=parameter_alpha*E+(1-parameter_alpha)*inv_E_inv;
        

        F=[alpha*inv_E*beta, I_mat.*(mu*k_0^2)-alpha*inv_E*alpha;
            beta*inv_E*beta-I_mat.*(mu*k_0^2), -beta*inv_E*alpha];       
        

        G=[-alpha*beta, alpha*alpha-(mu*k_0^2).*G_2;
            (mu*k_0^2).*G_1-beta*beta, alpha*beta];
        
        clear G_1 G_2 inv_E_inv
    end
else % zeta ~ 0 does not work
    % equation 34 [1]
    G=[(mu*k_0^2*sin(zeta))*inv_E_inv-alpha*beta, alpha*alpha-(mu*k_0^2)*(cos(zeta)^2*tb_E+sin(zeta)^2*inv_E_inv);
        (mu*k_0^2)*(cos(zeta)^2*bt_E+sin(zeta)^2*inv_E_inv)-beta*beta, alpha*beta-(mu*k_0^2*sin(zeta))*inv_E_inv];
        
    % equation 33 [1]
    F=[alpha*inv_E*beta-I_mat.*(mu*k_0^2*sin(zeta)), I_mat.*(mu*k_0^2)-alpha*inv_E*alpha;
        beta*inv_E*beta-I_mat.*(mu*k_0^2), I_mat.*(mu*k_0^2*sin(zeta))-beta*inv_E*alpha];
end