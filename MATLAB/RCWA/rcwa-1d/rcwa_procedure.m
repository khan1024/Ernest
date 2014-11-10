% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% [1]
% M. G. Moharam, Drew A. Pommet, Eric B. Grann, and T. K. Gaylord, 
% "Stable implementation of the rigorous coupled-wave analysis for surface-relief
% gratings: enhanced transmittance matrix approach," J. Opt. Soc. Am. A 12, 1077-1086 (1995) 

% [2]
% M. G. Moharam, Eric B. Grann, Drew A. Pommet, and T. K. Gaylord, 
% "Formulation for stable and efficient implementation of the rigorous coupled-wave
% analysis of binary gratings," J. Opt. Soc. Am. A 12, 1068-1076 (1995)

% [3]
% Lifeng Li, "Use of Fourier series in the analysis of discontinuous periodic
% structures," J. Opt. Soc. Am. A 13, 1870-1876 (1996)

% [4]
% Nikolay M. Lyndin, Olivier Parriaux, and Alexander V. Tishchenko, 
% "Modal analysis and suppression of the Fourier modal method instabilities 
% in highly conductive gratings," J. Opt. Soc. Am. A 24, 3781-3788 (2007) 
%break
initialization2_2;

order_max=(number_of_orders-1)/2;

n=(-order_max:order_max)';

% position of the zero diffraction order 
p=(number_of_orders+1)/2;

% incident angle theta
theta=theta0*pi/180;

% wave vector k_0
k_0=2*pi/(lambda*1E-6);

% psi_0 and phi_0 definition
if polarization==0
    psi=psi0*pi/180; % radian conversion
    if theta==0
        theta=1E-20; % due to some instabilities
    end
    phi=phi0*pi/180; % radian conversion
else
    psi=0;
    phi=0;
end

% grating vector K
K=2*pi/(Lambda*1E-6);

M=number_of_orders-1;

% identity matrix and matrix of zeros
I=eye(number_of_orders);
zero=zeros(number_of_orders);

% Floquet theorem, eg. (6) in [2]
k_x=k_0*n1*sin(theta)*cos(phi)+n*K;

% matrices K_x and K_x_2
k_x_2=k_x.^2;
K_x=diag(k_x/k_0);
K_x_2=K_x^2;

% k_y, ..., matrix k_y_2_I
k_y=k_0*n1*sin(theta)*sin(phi);
k_y_k_0=k_y/k_0;

k_y_2=k_y^2;
k_y_2_I=k_y_2*I;

% inclination angle phi_i
if k_y==0
    phi_i=zeros(number_of_orders,1);
else
    phi_i=atan(k_y./k_x);
    diff_phi_i=phi_i*180/pi;
end


% eq. (7) in [2], region 1
for m=1:number_of_orders
    if sqrt(k_x_2(m)+k_y_2)<k_0*n1
        k_1_z(m)=sqrt((k_0*n1)^2-k_x_2(m)-k_y_2);
    else
        k_1_z(m)=-1i*sqrt((k_x_2(m)+k_y_2)-(k_0*n1)^2);
    end
end

% eq. (7), region 3

for m=1:number_of_orders
    if sqrt(k_x_2(m)+k_y_2)<k_0*n3
        k_3_z(m)=sqrt((k_0*n3)^2-k_x_2(m)-k_y_2);
    else
        k_3_z(m)=-1i*sqrt((k_x_2(m)+k_y_2)-(k_0*n3)^2);
    end
end

% some changes in $k_{1,z}$ a $k_{3,z}$
    
k_1_z_n=zeros(1,number_of_orders);    
k_3_z_n=zeros(1,number_of_orders);

for i=1:1:number_of_orders        
    if imag(k_1_z(i))==0            
        k_1_z_n(i)=-k_1_z(i);     
    else        
        k_1_z_n(i)=k_1_z(i);        
    end
    
    if imag(k_3_z(i))==0            
        k_3_z_n(i)=-k_3_z(i);        
    else        
        k_3_z_n(i)=k_3_z(i);        
    end    
end

% x,y,z base matrix instabilities
if polarization==0
    if old_change_matrix_base==2 && change_matrix_base==1
        change_matrix_base=2;
    end
end
if change_matrix_base==2 && polarization==0
    %if theta0~=0 && phi0~=0
    if  length(find(k_1_z==0))>=2 || length(find(k_3_z==0))>=2
        disp('change_matrix_base=1')
        change_matrix_base=1;
    end
end

% x,y,z base in conical diffraction
if polarization==0 && matrix_algorithm==1 && change_matrix_base==2
    if set_calculation_field==1
        disp('change_matrix_base=1')
        change_matrix_base=1;
    end
end

% initialization of matrix algorithms

if matrix_algorithm==1
    S_matrix_initialization;
elseif matrix_algorithm==2
    T_matrix_initialization;
end


% planar diffraction angle=(-90,90)
% conical diffraction angle=(0,90)

D_R_angles_temp=zeros(number_of_orders,1);    
D_T_angles_temp=zeros(number_of_orders,1);    
for i=1:number_of_orders         
    D_R_angles_temp(i)=conj(-atan(sqrt(k_x(i)^2+k_y^2)/k_1_z_n(i))*180/pi)';        
    D_T_angles_temp(i)=conj(-atan(sqrt(k_x(i)^2+k_y^2)/k_3_z_n(i))*180/pi)';        
        
    if polarization~=0            
        if k_x(i)<0               
            D_R_angles_temp(i)=-conj(D_R_angles_temp(i));            
        else            
            D_T_angles_temp(i)=-conj(D_T_angles_temp(i));            
        end        
    end
    
end

D_R_angles=D_R_angles_temp(imag(D_R_angles_temp)>-1e-15);    
D_T_angles=D_T_angles_temp(imag(D_T_angles_temp)>-1e-15);
 
for l=number_of_layers:-1:1

    epsg=ng(1)^2;
    epsr=nr(1)^2;
    
    switch switch_permittivity
        case 0
            eps_fourier_series;  
        case 1 % grating from file        
            
            if grating_relative_dimensions==1            
                eps_fourier_series_grating_10_relative;        
            elseif grating_relative_dimensions==2            
                eps_fourier_series_grating_10_absolute;        
            end
    end

    
    % Toeplitz matrices
    E_l(:,:,l)=toeplitz(fliplr(v(1:number_of_orders)),v(number_of_orders:2*number_of_orders-1));
    A_l(:,:,l)=toeplitz(fliplr(i_v(1:number_of_orders)),i_v(number_of_orders:2*number_of_orders-1));
    %clear Sinc v_m v_0 v i_vm i_vo i_vp v_p i_v
    
    if polarization==0 % conical diffraction ------------------------------
        
        % eigenvalue equation for U_x, eq. (60) in [2]
        A=K_x_2-E_l(:,:,l);            
        eigenvalue_equation_1=k_y_2_I./k_0./k_0+K_x_2-E_l(:,:,l);             
        [W1(:,:,l),Q_pom1]=eig(eigenvalue_equation_1);            
        Q1(:,:,l)=diag(sqrt(diag(Q_pom1))); % =faster sqrt(Q_pom1)
        
        clear eigenvalue_equation_1 Q_pom1
        
        % eigenvalue equation for S_x, eq. (60) in [2]
        inv_E_l=inv(E_l(:,:,l));
        B=K_x*inv_E_l*K_x-I;            
        if faktorization==1                
            eigenvalue_equation_2=k_y_2_I./k_0./k_0+B/A_l(:,:,l);            
        elseif faktorization==2                
            eigenvalue_equation_2=k_y_2_I./k_0./k_0+B*E_l(:,:,l);            
        end        
        
        [W2(:,:,l),Q_pom2]=eig(eigenvalue_equation_2);            
        Q2(:,:,l)=diag(sqrt(diag(Q_pom2))); % =faster sqrt(Q_pom2)     
            
        clear eigenvalue_equation_2 Q_pom2        
            
        inv_A=inv(A);            
        inv_B=inv(B);        
        
        % eqs. (65) in [2]
        V11(:,:,l)=inv_A*W1(:,:,l)*Q1(:,:,l);            
        V12(:,:,l)=k_y_k_0*inv_A*K_x*W2(:,:,l);            
        V21(:,:,l)=k_y_k_0*inv_B*K_x*inv_E_l*W1(:,:,l);            
        V22(:,:,l)=inv_B*W2(:,:,l)*Q2(:,:,l);            
            
        if matrix_algorithm==1 && change_matrix_base==2                
            WW(:,:,l)=[zero, W2(:,:,l); V11(:,:,l), V12(:,:,l)];                
            VV(:,:,l)=-[W1(:,:,l), zero; V21(:,:,l), V22(:,:,l)];            
        end        
        
        if change_matrix_base==1
            % eqs. (72) in [2]
            Vss(:,:,l)=F_c*V11(:,:,l);
            Wss(:,:,l)=F_c*W1(:,:,l)+F_s*V21(:,:,l);
            Vsp(:,:,l)=F_c*V12(:,:,l)-F_s*W2(:,:,l);
            Wsp(:,:,l)=F_s*V22(:,:,l);
            Wpp(:,:,l)=F_c*V22(:,:,l);
            Vpp(:,:,l)=F_c*W2(:,:,l)+F_s*V12(:,:,l);
            Wps(:,:,l)=F_c*V21(:,:,l)-F_s*W1(:,:,l);
            Vps(:,:,l)=F_s*V11(:,:,l);
        end
            
    elseif polarization==1 % TM polarization ------------------------------         
        %%{
        % eq. 35 in [2]
        if faktorization==1 % correct factorization
            eigenvalue_equation=A_l(:,:,l)\(K_x*(E_l(:,:,l)\K_x)-I);
        else % incorrect faktorization      
            eigenvalue_equation=E_l(:,:,l)*(K_x*A_l(:,:,l)*K_x-I);  
        end

        [W(:,:,l),Q_2]=eig(eigenvalue_equation);
        Q(:,:,l)=diag(sqrt(diag(Q_2))); % eigenvalue matrix, =faster sqrt(Q_2)

	
        %{
        for k=1:1:number_of_orders
            if imag(Q(k,k,l))<0
                Q(k,k,l)=-Q(k,k,l);
            end            
            
            if (abs(imag(Q(k,k,l)))<abs(real(Q(k,k,l))))
                if real(Q(k,k,l))<0
                    Q(k,k,l)=-Q(k,k,l);
                end                
            end            
        end
        %}
        
        V(:,:,l)=A_l(:,:,l)*W(:,:,l)*Q(:,:,l);
        
        %%{
        % spurious modes [4]
        for k=1:1:number_of_orders            
            if imag(Q(k,k,l))<-20                
                Q(k,k,l)=Q(k,k,l)+1e40;
            end
        end
    
        %%}
        %%}
        %RCWA_NVM;
    elseif polarization==2 % TE polarization ------------------------------
        
        % eq. (15) in [2]
        eigenvalue_equation=K_x_2-E_l(:,:,l);
        [W(:,:,l),Q_2]=eig(eigenvalue_equation);

        clear eigenvalue_equation
        Q(:,:,l)=diag(sqrt(diag(Q_2))); % eigenvalue matrix, =faster sqrt(Q_2)
        V(:,:,l)=W(:,:,l)*Q(:,:,l); % eigenvector matrix - field H
    end

end

if exist('get_local_absoption')==0
    get_local_absoption=2;
end

if get_local_absoption==1 || get_electric_intensity==1
    set_calculation_field=1;
    plot_field_region_2=1;
end

if get_poynting_vector==1
    set_calculation_field=1;
    plot_field_region_2=1;
end

if polarization==0        
    if matrix_algorithm==1            
        if change_matrix_base==2
            if set_calculation_field==1
                error('set change_matrix_base=1')
            else
                S_matrix_W_S_grating_conical; %
            end
        elseif change_matrix_base==1                
            %if set_calculation_field==2
                %S_matrix_W_S_grating_conical;
            %else
            S_Li_conical_final;
            %end
        end        
    elseif matrix_algorithm==2
        T_matrix_algorithm;
    end    
else    
    if matrix_algorithm==1
        if set_calculation_field==2
            S_matrix_W_S_grating; % fast algorithm, but
        else            
            S_Li_final; % we need to know c coefficients           
        end        
    else        
        T_matrix_algorithm;        
    end
end

if polarization==1 %TM polarization        
    D_R=(R.*conj(R).*real(k_1_z_k_0/(n1*cos(theta))));
    D_T=(T.*conj(T).*(real(k_3_z_k_0/((n3)^2)))/(cos(theta)/n1));
    
elseif polarization==2 %TE polarization
    D_R=real((R.*conj(R).*real(k_1_z_k_0/(n1*cos(theta)))));
    D_T=real((T.*conj(T).*real(k_3_z_k_0/(n1*cos(theta)))));
    
else % conical diffraction
    if change_matrix_base==1
        R_x=(k_1_z_n'.*cos(phi_i).*R_p./(-1i*k_0*n1*n1)-sin(phi_i).*R_s);
        R_y=cos(phi_i).*R_s+k_1_z_n'.*sin(phi_i).*R_p./(-1i*k_0*n1*n1);
            
        if imag(n3)==0 % ? n1=imag
            T_x=(-1i*n1)*k_3_z_n'.*cos(phi_i).*T_p./(k_0*n3*n3)-sin(phi_i).*T_s;
            T_y=cos(phi_i).*T_s+(-1i*n1)*k_3_z_n'.*sin(phi_i).*T_p./(k_0*n3*n3);
            
        else
            T_x=(1i*n1)*conj(k_3_z_n)'.*cos(phi_i).*T_p./(k_0*n3*n3)-sin(phi_i).*T_s;
            T_y=cos(phi_i).*T_s+(1i*n1)*conj(k_3_z_n)'.*sin(phi_i).*T_p./(k_0*n3*n3);
            
        end
        R_z=(R_x.*k_x + R_y.*k_y)./k_1_z;
        T_z=-(T_x.*k_x + T_y.*k_y)./k_3_z;
        D_R=(R_s.*conj(R_s).*real(k_1_z_k_0/(n1*cos(theta))))+(R_p.*conj(R_p).*real(k_1_z_k_0/(n1^2)/(n1*cos(theta))));            
        D_T=(T_s.*conj(T_s).*real(k_3_z_k_0/(n1*cos(theta))))+(T_p.*conj(T_p).*real(k_3_z_k_0/(n3^2)/(n1*cos(theta))));
        
    elseif change_matrix_base==2
        D_R=(abs(R_x).^2 + abs(R_y).^2 + abs(R_z).^2).* real(k_1_z./(k_0*n1*cos(theta)));
        D_T=(abs(T_x).^2 + abs(T_y).^2 + abs(T_z).^2).* real(k_3_z./(k_0*n1*cos(theta)));
        %D_R=D_R';            
        %D_T=D_T';
        
    end
    
end

%-----------------------------------------------------------
% space harmonics and fields

if get_poynting_vector==1
    perioda=Lambda*1E-6;
    poynting_vector;
    plot_field_region_2=2;
end

if get_local_absoption==1 || get_electric_intensity==1
    perioda=Lambda*1E-6;
    local_absorption;
    plot_field_region_2=2;
end


if plot_field_region_1==1 || plot_field_region_1==2
    perioda=Lambda*1E-6;
    if plot_field_region_1==1
        %try
            field1_7;
        %catch
        %    fprintf('switch to older version');
        %    field1_6;
        %end
    end
    if plot_field_region_2==1
        %try
            field2_7;
        %catch
        %    fprintf('switch to older version');
        %    field2_5;
        %end 
    end
    if plot_field_region_3==1
        %try
            field3_6;
        %catch
        %    fprintf('switch to older version');
        %    field3_5;
        %end
    end
    if plot_field_region_1==1 || plot_field_region_2==1 || plot_field_region_3==1
        new_figure=1;
        if polarization==0
            field_graphs_conical;
        elseif polarization==1 || polarization==2
            if get_poynting_vector==1
            else
                field_graphs;
            end
        end
    end
end