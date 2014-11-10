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

W_3_11=[I_mat zero; -1i*K_1_z zero];
W_3_12=[sin(psi)*I_mat zero ; 1i*n_1*cos(theta)*sin(psi)*I_mat zero];
W_3_21=[zero I_mat ; zero -1i*K_1_z./(n_1^2)];
W_3_22=[-1i*n_1*cos(psi)*I_mat zero ; cos(psi)*cos(theta)*I_mat zero];

for i=1:1:number_of_layers
    W_2_11(:,:,i)=[F_c*W_1(:,:,i)-F_s*W_2(:,:,i); F_c*V_2(:,:,i)+F_s*V_1(:,:,i)];
    W_2_12(:,:,i)=[F_c*W_1(:,:,i)-F_s*W_2(:,:,i); -F_c*V_2(:,:,i)-F_s*V_1(:,:,i)];
    W_2_21(:,:,i)=[F_c*V_1(:,:,i)-F_s*V_2(:,:,i); F_c*W_2(:,:,i)+F_s*W_1(:,:,i)];
    W_2_22(:,:,i)=[-F_c*V_1(:,:,i)+F_s*V_2(:,:,i); F_c*W_2(:,:,i)+F_s*W_1(:,:,i)];
end

W_1_11=[zero zero ; zero zero];
W_1_12=[I_mat zero ;1i*K_3_z zero];
W_1_21=[zero zero ; zero zero];
W_1_22=[zero I_mat; zero 1i*K_3_z/(n_3^2)];

%clear W_1 W_2 V_1 V_2
%---------------------------------------------------------

s_2=inv([W_3_11 -W_2_12(:,:,1) ; W_3_21 -W_2_22(:,:,1)])*[W_2_11(:,:,1) -W_3_12 ; W_2_21(:,:,1) -W_3_22];
for i=1:1:number_of_layers-1
    s_1(:,:,number_of_layers-i)=inv([W_2_11(:,:,i) -W_2_12(:,:,i+1) ; W_2_21(:,:,i) -W_2_22(:,:,i+1)])*[W_2_11(:,:,i+1) -W_2_12(:,:,i) ; W_2_21(:,:,i+1) -W_2_22(:,:,i)];
end
s_0=inv([W_2_11(:,:,number_of_layers) -W_1_12 ; W_2_21(:,:,number_of_layers) -W_1_22])*[W_1_11 -W_2_12(:,:,number_of_layers) ; W_1_21 -W_2_22(:,:,number_of_layers)];

for i=1:1:number_of_layers
    X_1(:,:,number_of_layers+1-i)=diag(exp(-k_0*diag(Q(1:n_sq,1:n_sq,i))*tloustka(i)));; %tloustka(i)
    X_2(:,:,number_of_layers+1-i)=diag(exp(-k_0*diag(Q(n_sq+1:2*n_sq,n_sq+1:2*n_sq,i))*tloustka(i))); %tloustka(i)
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

%---------------------------------------------------------
d_n1=zeros(2*n_sq,1);
d_n1((n_sq+1)/2,1)=1;
%---------------------------------------------------------
u_n1=R_ud*d_n1;
d_0=T_dd*d_n1;

R_s=u_n1(1:n_sq);
R_p=u_n1(n_sq+1:2*n_sq);
T_s=d_0(1:n_sq);
T_p=d_0(n_sq+1:2*n_sq);

%if vypocet_pole_2==1
%    S_pole_n_vrstev_konicka;
%end