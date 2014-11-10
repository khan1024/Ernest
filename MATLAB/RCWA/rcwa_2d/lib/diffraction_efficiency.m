% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic


%u_x=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
%u_y=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);

psi=pi/2;

u_x=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
u_y=cos(psi)*cos(theta)*sin(zeta+phi)+sin(psi)*cos(zeta+phi);

input_1_2=u_x;
input_2_2=u_y;

normalize_2=A_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_2_2^2+B_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_1_2^2+C_mn_3((n_sq+1)/2,(n_sq+1)/2)*(input_1_2*conj(input_2_2)+input_2_2*conj(input_1_2));

incid_2=zeros(2*n_sq,1);
incid_2((n_x-1)/2+(n_x-1)*n_x/2+1)=input_1_2/sqrt(normalize_2);
incid_2((n_x-1)/2+(n_x-1)*n_x/2+1+n_sq)=input_2_2/sqrt(normalize_2);
%---

if calculate_field==1
    R_ud_temp=R_ud(:,:,number_of_layers+2);
    clear R_ud
    R_ud=R_ud_temp;
    clear R_ud_temp
end

Li2_u_n1=R_ud*incid_2;
Li2_d_0=T_dd*incid_2;

R_x_Li2=Li2_u_n1(1:n_sq);
R_y_Li2=Li2_u_n1(n_sq+1:2*n_sq);
T_x_Li2=Li2_d_0(1:n_sq);
T_y_Li2=Li2_d_0(n_sq+1:2*n_sq);

% equation 34a, 34b [2] 
D_R_3=A_mn_3*(R_y_Li2.*conj(R_y_Li2))+B_mn_3*(R_x_Li2.*conj(R_x_Li2))+C_mn_3*(R_x_Li2.*conj(R_y_Li2)+R_y_Li2.*conj(R_x_Li2));
D_T_3=A_mn_1*(T_y_Li2.*conj(T_y_Li2))+B_mn_1*(T_x_Li2.*conj(T_x_Li2))+C_mn_1*(T_x_Li2.*conj(T_y_Li2)+T_y_Li2.*conj(T_x_Li2));

D_R_90=real(D_R_3);
D_T_90=real(D_T_3);

%----------------------------------------------------------
psi=0;

u_x=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
u_y=cos(psi)*cos(theta)*sin(zeta+phi)+sin(psi)*cos(zeta+phi);

input_1_2=u_x;
input_2_2=u_y;

normalize_2=A_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_2_2^2+B_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_1_2^2+C_mn_3((n_sq+1)/2,(n_sq+1)/2)*(input_1_2*conj(input_2_2)+input_2_2*conj(input_1_2));

incid_2=zeros(2*n_sq,1);
incid_2((n_x-1)/2+(n_x-1)*n_x/2+1)=input_1_2/sqrt(normalize_2);
incid_2((n_x-1)/2+(n_x-1)*n_x/2+1+n_sq)=input_2_2/sqrt(normalize_2);
%---

Li2_u_n1=R_ud*incid_2;
Li2_d_0=T_dd*incid_2;

R_x_Li2=Li2_u_n1(1:n_sq);
R_y_Li2=Li2_u_n1(n_sq+1:2*n_sq);
T_x_Li2=Li2_d_0(1:n_sq);
T_y_Li2=Li2_d_0(n_sq+1:2*n_sq);

% equation 34a, 34b [2] 
D_R_3=A_mn_3*(R_y_Li2.*conj(R_y_Li2))+B_mn_3*(R_x_Li2.*conj(R_x_Li2))+C_mn_3*(R_x_Li2.*conj(R_y_Li2)+R_y_Li2.*conj(R_x_Li2));
D_T_3=A_mn_1*(T_y_Li2.*conj(T_y_Li2))+B_mn_1*(T_x_Li2.*conj(T_x_Li2))+C_mn_1*(T_x_Li2.*conj(T_y_Li2)+T_y_Li2.*conj(T_x_Li2));

D_R_0=real(D_R_3);
D_T_0=real(D_T_3);

%----------------------------------------------------------
psi=psi0*pi/180;

u_x=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
u_y=cos(psi)*cos(theta)*sin(zeta+phi)+sin(psi)*cos(zeta+phi);

input_1_2=u_x;
input_2_2=u_y;

normalize_2=A_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_2_2^2+B_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_1_2^2+C_mn_3((n_sq+1)/2,(n_sq+1)/2)*(input_1_2*conj(input_2_2)+input_2_2*conj(input_1_2));

incid_2=zeros(2*n_sq,1);
incid_2((n_x-1)/2+(n_x-1)*n_x/2+1)=input_1_2/sqrt(normalize_2);
incid_2((n_x-1)/2+(n_x-1)*n_x/2+1+n_sq)=input_2_2/sqrt(normalize_2);
%---

Li2_u_n1=R_ud*incid_2;
Li2_d_0=T_dd*incid_2;

R_x_Li2=Li2_u_n1(1:n_sq);
R_y_Li2=Li2_u_n1(n_sq+1:2*n_sq);
T_x_Li2=Li2_d_0(1:n_sq);
T_y_Li2=Li2_d_0(n_sq+1:2*n_sq);

% equation 34a, 34b [2] 
D_R_3=A_mn_3*(R_y_Li2.*conj(R_y_Li2))+B_mn_3*(R_x_Li2.*conj(R_x_Li2))+C_mn_3*(R_x_Li2.*conj(R_y_Li2)+R_y_Li2.*conj(R_x_Li2));
D_T_3=A_mn_1*(T_y_Li2.*conj(T_y_Li2))+B_mn_1*(T_x_Li2.*conj(T_x_Li2))+C_mn_1*(T_x_Li2.*conj(T_y_Li2)+T_y_Li2.*conj(T_x_Li2));

D_R=real(D_R_3);
D_T=real(D_T_3);


% quick export of real diffraction orders

if measurement==0
    
    m=-N_X:1:N_X;
    n=m;
    diff_orders=zeros((2*N_X+1)^2,3);
    for i=1:1:length(m)    
        for j=1:1:length(n)        
            r_index=1+(m(j)+N_X)+(2*N_X+1)*(n(i)+N_X);        
            diff_orders(r_index,:)=[r_index,n(i),m(j)];    
        end        
    end
    
    real_gamma_1=[];
    real_gamma_3=[];
    
    % real reflected + % real transmitted
    [real_gamma_1]=find(imag(diag(gamma_1))==0);
    [real_gamma_3]=find(imag(diag(gamma_3))==0);
    
    total_D_R_0=[];
    total_D_R_90=[];
    total_D_R=[];
    
    for i=1:1:length(real_gamma_1)
        total_D_R_0=[total_D_R_0; diff_orders(real_gamma_1(i),2), diff_orders(real_gamma_1(i),3), D_R_0(real_gamma_1(i))];
        total_D_R_90=[total_D_R_90; diff_orders(real_gamma_1(i),2), diff_orders(real_gamma_1(i),3), D_R_90(real_gamma_1(i))];
        total_D_R=[total_D_R; diff_orders(real_gamma_1(i),2), diff_orders(real_gamma_1(i),3), D_R(real_gamma_1(i))];
    end
    
    total_D_T_0=[];
    total_D_T_90=[];
    total_D_T=[];
    
    for i=1:1:length(real_gamma_3)
        total_D_T_0=[total_D_T_0; diff_orders(real_gamma_3(i),2), diff_orders(real_gamma_3(i),3), D_T_0(real_gamma_3(i))];
        total_D_T_90=[total_D_T_90; diff_orders(real_gamma_3(i),2), diff_orders(real_gamma_3(i),3), D_T_90(real_gamma_3(i))];
        total_D_T=[total_D_T; diff_orders(real_gamma_3(i),2), diff_orders(real_gamma_3(i),3), D_T(real_gamma_3(i))];
    end   
    
end