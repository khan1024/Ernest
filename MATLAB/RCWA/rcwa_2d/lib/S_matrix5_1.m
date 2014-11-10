% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% S-maticový algoritmus upravený podle Liova článku
%----------------------------------------------------------

W_1=zeros(2*n_sq,2*n_sq,number_of_layers+1); % number_of_layers+2 is I 
W_2=zeros(2*n_sq,2*n_sq,number_of_layers+2);
%phi_minus=zeros(2*n_sq,2*n_sq,number_of_layers+1);

new_I_mat=sparse([I_mat zero; zero I_mat]);
new_Z_mat=sparse([zero zero; zero zero]);

%-------------------------------------------------
% matrix alocation
phi_minus_1=new_I_mat;

% phi_minus, W_1, W_2 v mrizce
for i=1:1:number_of_layers
    W_1(:,:,(number_of_layers+2)-i)=eigvec_E(:,:,i);
    W_2(:,:,(number_of_layers+2)-i)=eigvec_H(:,:,i);
    %phi_minus(:,:,(number_of_layers+2)-i)=diag(exp(-1i*diag(gamma(:,:,i))*thickness(i)));
    current_phi_minus=sparse(diag(exp(-1i*diag(eval(['gamma_2_' num2str(i)]))*thickness(i))));
    eval(['phi_minus_' num2str((number_of_layers+2)-i), '=current_phi_minus;']);
end

S_express=1;
%---------------------------------------------------------
calculate_field=1;clear phi_minus
if calculate_field==1
    s_11=zeros(2*n_sq,2*n_sq,number_of_layers);
    s_12=zeros(2*n_sq,2*n_sq,number_of_layers);
    s_21=zeros(2*n_sq,2*n_sq,number_of_layers);
    s_22=zeros(2*n_sq,2*n_sq,number_of_layers);
    for_inv_s_22=zeros(2*n_sq,2*n_sq);
    %phi_minus=zeros(2*n_sq,2*n_sq,number_of_layers);
end

% matrix alocation
% počáteční inicializace

%R_du=zeros(2*n_sq,2*n_sq);
%T_uu=zeros(2*n_sq,2*n_sq);
T_dd=zeros(2*n_sq,2*n_sq);

R_ud_tilde=zeros(2*n_sq,2*n_sq);
%R_du_tilde=zeros(2*n_sq,2*n_sq);
%T_uu_tilde=zeros(2*n_sq,2*n_sq);
T_dd_tilde=zeros(2*n_sq,2*n_sq);

% pre-inicialization

if calculate_field==1
    R_ud(:,:,1)=zeros(2*n_sq,2*n_sq);
else
    R_ud=new_Z_mat;
end

%R_du=new_Z_mat;
%T_uu=new_I_mat;
T_dd=new_I_mat;
%-------------------------------------------------

W_1(:,:,1)=sparse([I_mat zero ; zero I_mat]);
%W_1(:,:,2)=eigvec_E;
%W_1(:,:,number_of_layers+2)=sparse([I_mat zero ; zero I_mat]);

A_mn_3=sparse((I_mat*(k_0*n_1)^2-alpha*alpha)/(mu*k_0*cos(zeta)*gamma_1));
B_mn_3=sparse((I_mat*(k_0*n_1)^2-beta*beta)/(mu*k_0*cos(zeta)*gamma_1));
C_mn_3=sparse((alpha*beta-I_mat*(k_0*n_1)^2*sin(zeta))/(mu*k_0*cos(zeta)*gamma_1));

A_mn_1=sparse((I_mat*(k_0*n_3)^2-alpha*alpha)/(mu*k_0*cos(zeta)*gamma_3));
B_mn_1=sparse((I_mat*(k_0*n_3)^2-beta*beta)/(mu*k_0*cos(zeta)*gamma_3));
C_mn_1=sparse((alpha*beta-I_mat*(k_0*n_3)^2*sin(zeta))/(mu*k_0*cos(zeta)*gamma_3));

W_2(:,:,1)=[-C_mn_1, -A_mn_1; B_mn_1, C_mn_1]; % output
%W_2(:,:,2)=eigvec_H;
W_2(:,:,number_of_layers+2)=[-C_mn_3, -A_mn_3; B_mn_3, C_mn_3]; % input

for i=1:1:number_of_layers+1
    
    %R_ud_tilde=phi_minus(:,:,i)*R_ud*phi_minus(:,:,i);
    %T_dd_tilde=T_dd*phi_minus(:,:,i);
    
    if calculate_field==1
        %R_ud_vlnka=phi_minus(:,:,i)*R_ud(:,:,i)*phi_minus(:,:,i);
        %T_dd_vlnka=T_dd*phi_minus(:,:,i);
        R_ud_tilde=eval(['phi_minus_' num2str(i)])*R_ud(:,:,i)*eval(['phi_minus_' num2str(i)]);
        T_dd_tilde=T_dd*eval(['phi_minus_' num2str(i)]);
    else
    T_dd_tilde=T_dd*eval(['phi_minus_' num2str(i)]);
    R_ud_tilde=eval(['phi_minus_' num2str(i)])*R_ud*eval(['phi_minus_' num2str(i)]);
    end
    
    %if S_express==1
    %else
    %    T_uu_tilde=phi_minus*T_uu;
    %end
    
    if i==number_of_layers+1 % W_1(:,:,number_of_layers+2)=sparse([I_mat zero ; zero I_mat])
        Q_01=W_1(:,:,i);
    else
        Q_01=(W_1(:,:,i+1))\W_1(:,:,i);
    end
    
    Q_02=(W_2(:,:,i+1))\W_2(:,:,i);
    
    %Q_01=inv(W_1(:,:,i+1))*W_1(:,:,i);
    %Q_02=inv(W_2(:,:,i+1))*W_2(:,:,i);

    F=Q_01*(new_I_mat+R_ud_tilde);
    G=Q_02*(new_I_mat-R_ud_tilde);
    
    tau=inv(F+G);
    
    if calculate_field==1
        R_ud(:,:,i+1)=new_I_mat-2*G*tau;
    else
        R_ud=new_I_mat-2*G*tau;
    end
    T_dd=2*T_dd_tilde*tau;
    
    %if S_express==1
    %else
    %    T_uu(:,:,i+1)=(F*tau*Q_02+G*tau*Q_01)*T_uu_tilde(:,:,i);
    %    R_du(:,:,i+1)=R_du(:,:,i)+T_dd_tilde(:,:,i)*tau*(Q_02-Q_01)*T_uu_tilde(:,:,i);
    %end
    if calculate_field==1
        t_11=1/2*(Q_01+Q_02);    
        t_12=1/2*(Q_01-Q_02);    
    
        if i==1        
            for_inv_s_22=t_11;    
        end        
        
        inv_t_22=inv(t_11);    
        inv_t_22_t_12=inv_t_22*t_12;    
    
        s_11(:,:,i)=t_11-t_12*inv_t_22_t_12; % t_21=t_12    
        s_12(:,:,i)=t_12*inv_t_22;    
        s_21(:,:,i)=-inv_t_22_t_12;  % t_21=t_12    
        s_22(:,:,i)=inv_t_22;    
    end
end

if calculate_field==1
    
    psi=psi0*pi/180;
    u_x=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
    u_y=cos(psi)*cos(theta)*sin(zeta+phi)+sin(psi)*cos(zeta+phi);
    u_z=-cos(psi)*sin(theta);

    input_1_2=u_x;
    input_2_2=u_y;

    normalize_2=A_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_2_2^2+B_mn_3((n_sq+1)/2,(n_sq+1)/2)*input_1_2^2+C_mn_3((n_sq+1)/2,(n_sq+1)/2)*(input_1_2*conj(input_2_2)+input_2_2*conj(input_1_2));

    incid_2=zeros(2*n_sq,1);
    incid_2((n_x-1)/2+(n_x-1)*n_x/2+1)=input_1_2/sqrt(normalize_2);
    incid_2((n_x-1)/2+(n_x-1)*n_x/2+1+n_sq)=input_2_2/sqrt(normalize_2);

    d_n1=incid_2;


    R=R_ud(:,:,number_of_layers+2)*d_n1;
    T=T_dd*d_n1;
    
    T_x_Li2=T(1:n_sq);
    T_y_Li2=T(n_sq+1:2*n_sq);
    
    c_m(:,:,number_of_layers)=s_12(:,:,1)*for_inv_s_22*T; % u_1_0
    c_p(:,:,1)=s_21(:,:,end)*(s_11(:,:,end)\(R-s_12(:,:,end)*d_n1))+s_22(:,:,end)*d_n1;
    %d_n_plus_1_n(:,:,number_of_layers)=sparse(phi_minus(:,:,end))*c_p(:,:,1); % výpočet do dalších částí
    d_n_plus_1_n(:,:,number_of_layers)=eval(['phi_minus_' num2str(number_of_layers+1)])*c_p(:,:,1); % výpočet do dalších částí

    if number_of_layers>=2
        for i=number_of_layers-1:-1:1
            c_m(:,:,number_of_layers-i)=R_ud(:,:,i+2)*d_n_plus_1_n(:,:,i+1);
            c_p(:,:,number_of_layers+1-i)=s_21(:,:,i+1)*(s_11(:,:,i+1)\(c_m(:,:,number_of_layers-i)-s_12(:,:,i+1)*d_n_plus_1_n(:,:,i+1)))+s_22(:,:,i+1)*d_n_plus_1_n(:,:,i+1);
            %d_n_plus_1_n(:,:,i)=sparse(phi_minus(:,:,i+1))*c_p(:,:,number_of_layers+1-i); % výpočet do dalších částí
            d_n_plus_1_n(:,:,i)=eval(['phi_minus_' num2str(i+1)])*c_p(:,:,number_of_layers+1-i); % výpočet do dalších částí 
        end
    end
    
end

%clear R_ud T_dd matrix_G T_dd_vlnka R_ud_vlnka matrix_F Q_01 Q_02 R_ud_vlnka phi_minus tau
clear s_11 s_12 s_21 s_22

set_waveguide_framework=2; % 1-yes, 2-no

xy_slice_view=2;
xz_slice_view=2;
yz_slice_view=2;

symmetry=2;
if calculate_field==1
    if view_field_slice==1
        xy_slice;
    elseif view_field_slice==2
        xz_slice;
    elseif view_field_slice==3
        yz_slice;
    end
end
