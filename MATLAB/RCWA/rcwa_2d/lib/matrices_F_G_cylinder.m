% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

transformation_4D_to_2D;

if factorization_method==0
% equation 1b, Lalanne, JOSA A, vol. 14, p. 1592, 1997

inv_E=inv(E);
inv_A=inv(A);
        
G_1=parameter_alpha*inv_A+(1-parameter_alpha)*E;        
G_2=parameter_alpha*E+(1-parameter_alpha)*inv_A;

F=[alpha*inv_E*beta, I_mat.*(mu*k_0^2)-alpha*inv_E*alpha;
    beta*inv_E*beta-I_mat.*(mu*k_0^2), -I_mat.*beta*inv_E*alpha];       
        
G=[-alpha*beta, alpha*alpha-(mu*k_0^2).*G_2;
    (mu*k_0^2).*G_1-beta*beta, alpha*beta];

elseif factorization_method==1

mex_speedup=1;
dx=0.01;
px=(0:dx:Lambda_x)-Lambda_x/2;
py=(0:dx:Lambda_y)-Lambda_y/2;

disc_x=px;    
disc_y=py;

switch factorization_method
    case 0 % Li factorization
    case 1 % moje normal vector method
        Nf=length(dx/2:dx:Lambda_x-dx/2); Nx=(1:Nf)'*ones(1,Nf); Ny=ones(Nf,1)*(1:Nf);
        phi_n=angle(Nx-Nf/2+1i*(Ny-Nf/2));
        
        th=phi_n;
        a=cos(th);
        b=sin(th);
        
        % obraceně zadané 
        N_x_input=b; %a
        N_y_input=a; %b
        %[cone_x,cone_y]=meshgrid((dx/2:dx:Lambda_x-dx/2)-Lambda_x/2,(dx/2:dx:Lambda_y-dx/2)-Lambda_y/2);
        %cone_z=-sqrt((cone_x).^2+(cone_y).^2);
        %[N_x_input,N_y_input] = gradient(cone_z,dx,dx);
        %figure;contour(cone_x,cone_y,cone_z), hold on, quiver(cone_x,cone_y,N_x_input,N_y_input)
    case 2
    case 3
end

use_grid_asr_speed_up=1;
N_xx=asr_fourier_2D(disc_x,disc_y,disc_x,disc_y,N_x_input.*N_x_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);
N_yy=asr_fourier_2D(disc_x,disc_y,disc_x,disc_y,N_y_input.*N_y_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);    
N_xy=asr_fourier_2D(disc_x,disc_y,disc_x,disc_y,N_x_input.*N_y_input,N_X,N_Y,1,mex_speedup,use_grid_asr_speed_up);

% equation 34 [1]

E_inv=A;

delta_epsilon=E-inv(E_inv);        
delta_N_xy=delta_epsilon*N_xy;
        
        
%G=[-alpha*beta, alpha*alpha-(mu*k_0^2).*tb_E;
%    (mu*k_0^2).*bt_E-beta*beta, alpha*beta];        
        
G=[-alpha*beta+(mu*k_0^2).*delta_N_xy, alpha*alpha+(mu*k_0^2).*(-E+delta_epsilon*N_yy);
    (mu*k_0^2).*(E-delta_epsilon*N_xx)-beta*beta, alpha*beta-(mu*k_0^2).*delta_N_xy];        
        
% equation 33 [1]      
F=[alpha*inv_E*beta, I_mat.*(mu*k_0^2)-alpha*inv_E*alpha;
    beta*inv_E*beta-I_mat.*(mu*k_0^2), -beta*inv_E*alpha];
end