% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%[1] Li, JOSA A, Vol.14, No.10, 1997 - New formulation of the Fourier modal method for crossed surface-relief gratings 
%[2] Li, J. Opt. A: Pure Appl. Opt., Vol.5, 2003 - Fourier modal method for crossed
% anisotropic gratings with arbitrary permittivity and permeability tensors

exist_zeta0=exist('zeta0','var');
    
if exist_zeta0==0
    zeta0=0;
end

epsilon_1=n_1*n_1;
epsilon_3=n_3*n_3;

if theta0==0
    theta0=1E-10;
end
theta=theta0*pi/180;
phi=phi0*pi/180;
zeta=zeta0*pi/180;
psi=psi0*pi/180;
k_0=2*pi/lambda;

if grating==0
    if optimalization_grating_0==1
        optimalization_input;
    end
    Lambda_x=new_px_1(1,end)-new_px_1(1,1);
    Lambda_y=new_py_1(1,end,1)-new_py_1(1,1);
else
    homogenous_layer=zeros(number_of_layers,1);
end

K_1=2*pi/Lambda_x;
K_2=2*pi/Lambda_y;

K_x_m=K_1;
K_y_m=K_2;

% equation 11 [1]
k_p1=2*pi*sqrt(epsilon_1*mu)/lambda;
k_m1=2*pi*sqrt(epsilon_3*mu)/lambda;

alpha_0=k_p1*sin(theta)*cos(phi);
beta_0=k_p1*sin(theta)*sin(phi+zeta);
gamma_p00=k_p1*cos(theta);

m_K=(-N_X:1:N_X);
n_K=(-N_Y:1:N_Y);

n_x=2*N_X+1;
n_y=2*N_Y+1;
n_sq=n_x*n_y;
nPermCoeffs=2*n_x-1;
maxPermCoeff=(nPermCoeffs-1)/2;

I_mat=speye(n_sq);%eye(n_sq);
I_mat2=speye(2*n_sq);%eye(2*n_sq);
zero=sparse(zeros(n_sq));

exist_field_calculation=exist('calculate_field','var');
    
if exist_field_calculation==0
    calculate_field=2;
end
%------------------------------------------------------

alpha_m=zeros(2*N_X+1,1);
beta_m=zeros(2*N_X+1,1);

%alpha_1, alpha_3
%beta_1, beta_3

%alpha_1_m=zeros(2*N_X+1,1);
%alpha_3_m=zeros(2*N_X+1,1);
%beta_1_n=zeros(2*N_X+1,1);
%beta_3_n=zeros(2*N_X+1,1);

% equation 14 [1]
alpha_m=alpha_0+m_K*K_1;
beta_n=beta_0+n_K*K_2;

%alpha_1_m=2*pi*sqrt(epsilon_1*mu)/lambda*sin(theta)*cos(phi)+m_K*K_1;
%alpha_3_m=2*pi*sqrt(epsilon_3*mu)/lambda*sin(theta)*cos(phi)+m_K*K_1;

%beta_1_n=2*pi*sqrt(epsilon_1*mu)/lambda*sin(theta)*sin(phi+zeta)+n_K*K_2;
%beta_3_n=2*pi*sqrt(epsilon_3*mu)/lambda*sin(theta)*sin(phi+zeta)+n_K*K_2;

alpha=zeros(n_sq,n_sq);
beta=zeros(n_sq,n_sq);
%alpha_1=zeros(n_sq,n_sq);
%alpha_3=zeros(n_sq,n_sq);
%beta_1=zeros(n_sq,n_sq);
%beta_3=zeros(n_sq,n_sq);

% (alpha)_mn,jl=alpha_m*delta_mj*delta_nl
% (beta)_mn,jl=beta_n*delta_mj*delta_nl


% call mex-function or m-function
[alpha,beta]=alpha_beta_mex(n_x,alpha_m,beta_n);    
%[alpha_1,beta_1]=alpha_beta_mex(n_x,alpha_1_m,beta_1_n);
%[alpha_3,beta_3]=alpha_beta_mex(n_x,alpha_3_m,beta_3_n);

alpha=sparse(alpha);
beta=sparse(beta);

%alpha_1=sparse(alpha_1);
%beta_1=sparse(beta_1);

%alpha_3=sparse(alpha_3);
%beta_3=sparse(beta_3);

alpha_beta=alpha*beta;

phi_i=zeros(1,n_sq);

for i=1:1:n_sq
    if abs(alpha(i,i))>1E-8
        phi_i(i)=atan(beta(i,i)/alpha(i,i));
        if phi_i(i)<-1E-8
            phi_i(i)=phi_i(i)+pi;
        end
    else
        if abs(beta(i,i))<1E-8
            phi_i(i)=phi;
        else
            phi_i(i)=pi/2;
        end
    end
end

%-----------------------------------------------------
% nastavení z-ovych komponent difraktovanych radu --- reflexe
% equation 15 [1]

gamma_1=sparse(zeros(n_sq,n_sq));
gamma_3=sparse(zeros(n_sq,n_sq));

gamma_3=sqrt(I_mat*(2*pi/lambda*n_3)^2-sec(zeta)^2*(alpha^2+beta^2-2*alpha*beta*sin(zeta)));
gamma_1=sqrt(I_mat*(2*pi/lambda*n_1)^2-sec(zeta)^2*(alpha^2+beta^2-2*alpha*beta*sin(zeta)));

for i=1:1:n_sq
    %gamma_3(i,i)=sqrt((2*pi/lambda*n_3)^2-sec(zeta)^2*(alpha(i,i)^2+beta(i,i)^2-2*alpha(i,i)*beta(i,i)*sin(zeta)));
    %gamma_1(i,i)=sqrt((2*pi/lambda*n_1)^2-sec(zeta)^2*(alpha(i,i)^2+beta(i,i)^2-2*alpha(i,i)*beta(i,i)*sin(zeta)));
    
    % change of signs, equation 16 [1]
    
    % if imag(gamma)<0 (? nektere clanky maji obracene)
    if (imag(gamma_1(i,i))>0)
        gamma_1(i,i)=-1*gamma_1(i,i);
    end
    
    % vlastni uprava --- equation 16 [1]
    if abs(imag(gamma_1(i,i)))<abs(real(gamma_1(i,i)))
        if (real(gamma_1(i,i))<0)
            gamma_1(i,i)=-1*gamma_1(i,i);
        end
    end
    
    if (imag(gamma_3(i,i))>0)
        gamma_3(i,i)=-1*gamma_3(i,i);
    end
    
    if abs(imag(gamma_3(i,i)))<abs(real(gamma_3(i,i)))
        if (real(gamma_3(i,i))<0)
            gamma_3(i,i)=-1*gamma_3(i,i);
        end
    end
end

inv_E=zeros(n_sq,n_sq);
inv_E_inv=zeros(n_sq,n_sq);
G=zeros(2*n_sq,2*n_sq);
F=zeros(2*n_sq,2*n_sq);
FG=zeros(2*n_sq,2*n_sq);

%gamma=zeros(2*n_sq,2*n_sq,number_of_layers);

eigvec_E=zeros(2*n_sq,2*n_sq,number_of_layers);
eigvec_H=zeros(2*n_sq,2*n_sq,number_of_layers);

for layer=number_of_layers:-1:1
    if homogenous_layer(layer)==1
        rcwa_procedure_homogenous_layer;
        if calculate_field==1
            inv_E_hom=1/(current_hom_n(layer))^2.*I_mat;
            eval(['inv_E_' num2str(layer),'=inv_E_hom;']);
            clear inv_E_hom
        end
    else
        
        switch grating        
            case 1            
                epsilon_1=n_layer*n_layer;            
                epsilon_2=n_cylinder*n_cylinder;        
            case 2            
                epsilon_1=n_layer*n_layer;            
                epsilon_2=n_cylinder*n_cylinder;            
                epsilon_1d=0;            
                epsilon_1i=0;            
                epsilon_2d=-(abs(real(epsilon_1))+abs(imag(epsilon_1))*1i)+epsilon_2;            
                epsilon_2i=-(abs(real(1/epsilon_1))-imag(1/epsilon_1)*1i)+1/epsilon_2;    
        end        
        
        if grating==0
            if factorization_method==0 || nvm_method_(layer)==0
                matrices_F_G_rectangle;
            elseif (nvm_method_(layer)==1 && factorization_method==1) || (nvm_method_(layer)==1 && factorization_method==2)
                normal_vector_method;
            elseif (nvm_method_(layer)==1 && factorization_method==3) || (nvm_method_(layer)==1 && factorization_method==4)
                complex_polarization_basis;
            end
        elseif grating==1 || grating==2        
            matrices_F_G_cylinder;    
        end        
        
        FG=F*G;    
        clear F
    
        % equation 35 [1]    
        [eigvec_E(:,:,layer),gamma_temp]=eig(FG);    
        eval(['gamma_2_' num2str(layer), '=sparse(diag(sqrt(diag(gamma_temp))/(sqrt(mu)*k_0*cos(zeta))));']);    
    
        clear FG gamma_temp    
    
        %{
        % equation 37 [1]
        for i=1:1:2*n_sq
            if (imag(gamma(i,i,layer))>0)
                gamma(i,i,layer)=-1*gamma(i,i,layer);
            end            
            
            if abs(imag(gamma(i,i,layer)))<abs(real(gamma(i,i,layer)))            
                if (real(gamma(i,i,layer))<0)
                    gamma(i,i,layer)=-1*gamma(i,i,layer);
                end                
            end            
        end        
        %}    
    
        current_gamma=eval(['gamma_2_' num2str(layer)]);
    
        for k=1:1:2*n_sq        
            if imag(current_gamma(k,k))>0            
                current_gamma(k,k)=-current_gamma(k,k);        
            end            
            
            if (abs(imag(current_gamma(k,k)))<abs(real(current_gamma(k,k))))            
                if real(current_gamma(k,k))<0                
                    current_gamma(k,k)=-current_gamma(k,k);            
                end                
            end            
        end        
        
        eval(['gamma_2_' num2str(layer), '=current_gamma;']);    
    
        % equation 36 [1]    
        %eigvec_H(:,:,layer)=G*eigvec_E(:,:,layer)*sparse(I_mat2/(k_0*gamma(:,:,layer))); % mu*k_0    
        eigvec_H(:,:,layer)=G*eigvec_E(:,:,layer)*sparse(sec(zeta)*I_mat2/(k_0*eval(['gamma_2_' num2str(layer)]))); % mu*k_0    
        clear G
        
        if calculate_field==1
            eval(['inv_E_' num2str(layer),'=inv_E;']);
        end
        clear inv_E
    end
end

if grating==0
    if optimalization_grating_0==1
        optimalization_output;
    end
end

if calculate_field==1
    S_matrix5_1;
else
    S_matrix5;
end

diffraction_efficiency;