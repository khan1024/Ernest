% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

X=diag(exp(-k_0*diag(Q(:,:,layer))*tloustka(layer)));

% !! úprava na pozici 2,2, puvodne F_c*V_2-F_s*V_1
H=[F_c*W_1(:,:,layer)-F_s*W_2(:,:,layer),    F_c*W_1(:,:,layer)-F_s*W_2(:,:,layer) ; ...
   -(F_c*V_2(:,:,layer)+F_s*V_1(:,:,layer)), F_c*V_2(:,:,layer)+F_s*V_1(:,:,layer) ; ...
   -F_c*V_1(:,:,layer)+F_s*V_2(:,:,layer),   -(-F_c*V_1(:,:,layer)+F_s*V_2(:,:,layer)) ; ...
    F_c*W_2(:,:,layer)+F_s*W_1(:,:,layer),    F_c*W_2(:,:,layer)+F_s*W_1(:,:,layer)];

%[H_U,H_S,H_V]=svd(H);
%inv_H=H_V*inv(H_S)*H_U';
inv_H=inv(H);

a_b=inv_H*f_g;

a=a_b(1:2*n_sq, : );
b=a_b((2*n_sq+1):4*n_sq, : );

%[a_U,a_S,a_V]=svd(a);
%a_1=a_V*inv(a_S)*a_U';
a_1=inv(a);

a_1_X=a_1*X;
X_b_a_1_X=X*b*a_1_X; % rodis davá X na druhou stranu

f_g=H*[I_mat zero ; zero I_mat; X_b_a_1_X];

nove_T=nove_T*a_1_X;