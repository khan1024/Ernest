% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% S-maticový algoritmus upravený podle Liova článku

% počáteční inicializace
W_2_11(:,:,number_of_layers)=zeros(2*n_sq,2*n_sq);
W_2_12(:,:,number_of_layers)=zeros(2*n_sq,2*n_sq);
W_2_21(:,:,number_of_layers)=zeros(2*n_sq,2*n_sq);
W_2_22(:,:,number_of_layers)=zeros(2*n_sq,2*n_sq);
X_1(:,:,number_of_layers)=zeros(n_sq,n_sq);
X_2(:,:,number_of_layers)=zeros(n_sq,n_sq);
b_uu(:,:,number_of_layers)=zeros(2*n_sq,2*n_sq);
b_ud(:,:,number_of_layers)=zeros(2*n_sq,2*n_sq);
b_du(:,:,number_of_layers)=zeros(2*n_sq,2*n_sq);
b_dd(:,:,number_of_layers)=zeros(2*n_sq,2*n_sq);
a_uu(:,:,number_of_layers+1)=zeros(2*n_sq,2*n_sq);
a_ud(:,:,number_of_layers+1)=zeros(2*n_sq,2*n_sq);
a_du(:,:,number_of_layers+1)=zeros(2*n_sq,2*n_sq);
a_dd(:,:,number_of_layers+1)=zeros(2*n_sq,2*n_sq);

if number_of_layers>2
s_1(:,:,number_of_layers-1)=zeros(4*n_sq,4*n_sq);
s_1_l(:,:,number_of_layers-1)=zeros(4*n_sq,4*n_sq);
end

%----------------------------------------------------------
thickness=tloustka;

k_x00=k_0*n_1*sin(theta)*cos(phi); % vymazano k_0 díky vydělení
k_y00=k_0*n_1*sin(theta)*sin(phi+zeta); % vymazano k_0 díky vydělení
k_z00=k_0*n_1*cos(theta); % vymazano k_0 díky vydělení

A_mn_3=(I_mat*(k_0*n_1)^2-alpha*alpha)/(mu*k_0*cos(zeta)*gamma_1);
B_mn_3=(I_mat*(k_0*n_1)^2-beta*beta)/(mu*k_0*cos(zeta)*gamma_1);
C_mn_3=(alpha*beta-I_mat*(k_0*n_1)^2*sin(zeta))/(mu*k_0*cos(zeta)*gamma_1);

W_3_11=[I_mat zero; zero I_mat];
W_3_12=[I_mat zero ; zero I_mat];
W_3_21=[C_mn_3, A_mn_3; -B_mn_3, -C_mn_3];
W_3_22=[-C_mn_3, -A_mn_3; B_mn_3, C_mn_3];

for i=1:1:number_of_layers
    W_2_11(:,:,i)=eigvec_E(:,:,i);
    W_2_12(:,:,i)=eigvec_E(:,:,i);
    W_2_21(:,:,i)=-eigvec_H(:,:,i);
    W_2_22(:,:,i)=eigvec_H(:,:,i);
end

A_mn_1=(I_mat*(k_0*n_3)^2-alpha*alpha)/(mu*k_0*cos(zeta)*gamma_3);
B_mn_1=(I_mat*(k_0*n_3)^2-beta*beta)/(mu*k_0*cos(zeta)*gamma_3);
C_mn_1=(alpha*beta-I_mat*(k_0*n_3)^2*sin(zeta))/(mu*k_0*cos(zeta)*gamma_3);

W_1_11=[zero zero ; zero zero];
W_1_12=[I_mat zero ;zero, I_mat];
W_1_21=[zero zero ; zero zero];
W_1_22=[-C_mn_1, -A_mn_1; B_mn_1, C_mn_1];

%clear W_1 W_2 V_1 V_2
%---------------------------------------------------------

s_2=inv([W_3_11 -W_2_12(:,:,1) ; W_3_21 -W_2_22(:,:,1)])*[W_2_11(:,:,1) -W_3_12 ; W_2_21(:,:,1) -W_3_22];
for i=1:1:number_of_layers-1
    s_1(:,:,number_of_layers-i)=inv([W_2_11(:,:,i) -W_2_12(:,:,i+1) ; W_2_21(:,:,i) -W_2_22(:,:,i+1)])*[W_2_11(:,:,i+1) -W_2_12(:,:,i) ; W_2_21(:,:,i+1) -W_2_22(:,:,i)];
end
s_0=inv([W_2_11(:,:,number_of_layers) -W_1_12 ; W_2_21(:,:,number_of_layers) -W_1_22])*[W_1_11 -W_2_12(:,:,number_of_layers) ; W_1_21 -W_2_22(:,:,number_of_layers)];

for i=1:1:number_of_layers
    X_1(:,:,number_of_layers+1-i)=diag(exp(-1i*diag(gamma(1:n_sq,1:n_sq,i))*thickness(i))); %tloustka(i)
    X_2(:,:,number_of_layers+1-i)=diag(exp(-1i*diag(gamma(n_sq+1:2*n_sq,n_sq+1:2*n_sq,i))*thickness(i))); %tloustka(i)
end

for i=1:1:number_of_layers
    if i==number_of_layers
        prvni_m=[I_mat zero zero zero ; zero I_mat zero zero; zero zero X_1(:,:,number_of_layers) zero; zero zero zero X_2(:,:,number_of_layers)];
        druha_m=[X_1(:,:,number_of_layers) zero zero zero ; zero X_2(:,:,number_of_layers) zero zero; zero zero I_mat zero; zero zero zero I_mat];
        
        s_2_l=prvni_m*s_2*druha_m;
    else
        prvni_m=[I_mat zero zero zero ; zero I_mat zero zero; zero zero X_1(:,:,i) zero; zero zero zero X_2(:,:,i)];
        druha_m=[X_1(:,:,i) zero zero zero ; zero X_2(:,:,i) zero zero; zero zero I_mat zero; zero zero zero I_mat];
        
        s_1_l(:,:,i)=prvni_m*s_1(:,:,i)*druha_m;
    end
end

%---------------------------------------------------------

for i=1:1:number_of_layers-1
    b_uu(:,:,i)=s_1_l(1:2*n_sq,1:2*n_sq,i);
    b_ud(:,:,i)=s_1_l(1:2*n_sq,1+2*n_sq:4*n_sq,i);
    b_du(:,:,i)=s_1_l(1+2*n_sq:4*n_sq,1:2*n_sq,i);
    b_dd(:,:,i)=s_1_l(2*n_sq+1:4*n_sq,2*n_sq+1:4*n_sq,i);
end

b_uu(:,:,number_of_layers)=s_2_l(1:2*n_sq,1:2*n_sq);
b_ud(:,:,number_of_layers)=s_2_l(1:2*n_sq,1+2*n_sq:4*n_sq);
b_du(:,:,number_of_layers)=s_2_l(1+2*n_sq:4*n_sq,1:2*n_sq);
b_dd(:,:,number_of_layers)=s_2_l(2*n_sq+1:4*n_sq,2*n_sq+1:4*n_sq);

a_uu(:,:,1)=s_0(1:2*n_sq,1:2*n_sq);
a_ud(:,:,1)=s_0(1:2*n_sq,1+2*n_sq:4*n_sq);
a_du(:,:,1)=s_0(1+2*n_sq:4*n_sq,1:2*n_sq);
a_dd(:,:,1)=s_0(2*n_sq+1:4*n_sq,2*n_sq+1:4*n_sq);

nove_I_mat=[I_mat zero; zero I_mat];
for i=1:1:number_of_layers
    %[dd_U,dd_S,dd_V]=svd(nove_I_mat-a_ud(:,:,i)*b_du(:,:,i));
    %dd=dd_V*inv(dd_S)*dd_U';
    
    %[dd1_U,dd1_S,dd1_V]=svd(nove_I_mat-b_du(:,:,i)*a_ud(:,:,i));
    %dd1=dd1_V*inv(dd1_S)*dd1_U';
    
    a_uu(:,:,i+1)=b_uu(:,:,i)*inv(nove_I_mat-a_ud(:,:,i)*b_du(:,:,i))*a_uu(:,:,i);
    a_ud(:,:,i+1)=b_ud(:,:,i)+b_uu(:,:,i)*a_ud(:,:,i)*inv(nove_I_mat-b_du(:,:,i)*a_ud(:,:,i))*b_dd(:,:,i);
    a_du(:,:,i+1)=a_du(:,:,i)+a_dd(:,:,i)*b_du(:,:,i)*inv(nove_I_mat-a_ud(:,:,i)*b_du(:,:,i))*a_uu(:,:,i);
    a_dd(:,:,i+1)=a_dd(:,:,i)*inv(nove_I_mat-b_du(:,:,i)*a_ud(:,:,i))*b_dd(:,:,i);  
end

T_dd=a_dd(:,:,number_of_layers+1);
R_ud=a_ud(:,:,number_of_layers+1);

psi=pi/2;
% nova_uprava
u_x=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
u_y=cos(psi)*cos(theta)*sin(zeta+phi)+sin(psi)*cos(zeta+phi);

input_1_2=u_x;
input_2_2=u_y;

normalize_2=A_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_2_2^2+B_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_1_2^2+C_mn_3((n_sq+1)/2,(n_sq+1)/2)*(input_1_2*conj(input_2_2)+input_2_2*conj(input_1_2));

incid_2=zeros(2*n_sq,1);
incid_2((n_x-1)/2+(n_x-1)*n_x/2+1)=input_1_2/sqrt(normalize_2);
incid_2((n_x-1)/2+(n_x-1)*n_x/2+1+n_sq)=input_2_2/sqrt(normalize_2);
%---------------------------------------------------------

Li2_u_n1=R_ud*incid_2;
Li2_d_0=T_dd*incid_2;

R_x_Li2=Li2_u_n1(1:n_sq);
R_y_Li2=Li2_u_n1(n_sq+1:2*n_sq);
T_x_Li2=Li2_d_0(1:n_sq);
T_y_Li2=Li2_d_0(n_sq+1:2*n_sq);

% vzorec 34a, 34b [2] 
D_R_3=A_mn_3*(R_y_Li2.*conj(R_y_Li2))+B_mn_3*(R_x_Li2.*conj(R_x_Li2))+C_mn_3*(R_x_Li2.*conj(R_y_Li2)+R_y_Li2.*conj(R_x_Li2));
D_T_3=A_mn_1*(T_y_Li2.*conj(T_y_Li2))+B_mn_1*(T_x_Li2.*conj(T_x_Li2))+C_mn_1*(T_x_Li2.*conj(T_y_Li2)+T_y_Li2.*conj(T_x_Li2));

D_R_90=real(D_R_3);
D_T_90=real(D_T_3);

%----------------------------------------------------------
psi=0;
% nova_uprava
u_x=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
u_y=cos(psi)*cos(theta)*sin(zeta+phi)+sin(psi)*cos(zeta+phi);

input_1_2=u_x;
input_2_2=u_y;

normalize_2=A_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_2_2^2+B_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_1_2^2+C_mn_3((n_sq+1)/2,(n_sq+1)/2)*(input_1_2*conj(input_2_2)+input_2_2*conj(input_1_2));

incid_2=zeros(2*n_sq,1);
incid_2((n_x-1)/2+(n_x-1)*n_x/2+1)=input_1_2/sqrt(normalize_2);
incid_2((n_x-1)/2+(n_x-1)*n_x/2+1+n_sq)=input_2_2/sqrt(normalize_2);
%---------------------------------------------------------

Li2_u_n1=R_ud*incid_2;
Li2_d_0=T_dd*incid_2;

R_x_Li2=Li2_u_n1(1:n_sq);
R_y_Li2=Li2_u_n1(n_sq+1:2*n_sq);
T_x_Li2=Li2_d_0(1:n_sq);
T_y_Li2=Li2_d_0(n_sq+1:2*n_sq);

% vzorec 34a, 34b [2] 
D_R_3=A_mn_3*(R_y_Li2.*conj(R_y_Li2))+B_mn_3*(R_x_Li2.*conj(R_x_Li2))+C_mn_3*(R_x_Li2.*conj(R_y_Li2)+R_y_Li2.*conj(R_x_Li2));
D_T_3=A_mn_1*(T_y_Li2.*conj(T_y_Li2))+B_mn_1*(T_x_Li2.*conj(T_x_Li2))+C_mn_1*(T_x_Li2.*conj(T_y_Li2)+T_y_Li2.*conj(T_x_Li2));

D_R_0=real(D_R_3);
D_T_0=real(D_T_3);