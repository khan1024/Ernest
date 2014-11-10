% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% S-maticový algoritmus upravený podle Liova článku
%----------------------------------------------------------

W_1=zeros(2*n_sq,2*n_sq,number_of_layers+1); % number_of_layers+2 = I 
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
% matrix alocation
R_ud=zeros(2*n_sq,2*n_sq);
%R_du=zeros(2*n_sq,2*n_sq);
%T_uu=zeros(2*n_sq,2*n_sq);
T_dd=zeros(2*n_sq,2*n_sq);

R_ud_tilde=zeros(2*n_sq,2*n_sq);
%R_du_tilde=zeros(2*n_sq,2*n_sq);
%T_uu_tilde=zeros(2*n_sq,2*n_sq);
T_dd_tilde=zeros(2*n_sq,2*n_sq);

% pre-inicialization
R_ud=new_Z_mat;
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
    
    T_dd_tilde=T_dd*eval(['phi_minus_' num2str(i)]);
    R_ud_tilde=eval(['phi_minus_' num2str(i)])*R_ud*eval(['phi_minus_' num2str(i)]);
   
    
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
    
    R_ud=new_I_mat-2*G*tau;
    T_dd=2*T_dd_tilde*tau;
    
    %if S_express==1
    %else
    %    T_uu(:,:,i+1)=(F*tau*Q_02+G*tau*Q_01)*T_uu_tilde(:,:,i);
    %    R_du(:,:,i+1)=R_du(:,:,i)+T_dd_tilde(:,:,i)*tau*(Q_02-Q_01)*T_uu_tilde(:,:,i);
    %end
end
