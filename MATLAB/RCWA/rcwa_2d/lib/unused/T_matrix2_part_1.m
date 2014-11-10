% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

X=diag(exp(-k_0*diag(Q(:,:,layer))*tloustka(layer)));

K_x_T2=K_x*k_0;
K_y_T2=K_y*k_0;
K_1_z_T2=K_1_z*k_0;
K_3_z_T2=K_3_z*k_0;


f_g=[I_mat, zero; zero, I_mat; -K_x_T2*K_y_T2/(k_0*K_3_z_T2), -(K_y_T2.^2+K_3_z_T2.^2)/(k_0*K_3_z_T2); (K_x_T2.^2+K_3_z_T2.^2)/(k_0*K_3_z_T2), K_x_T2*K_y_T2/(k_0*K_3_z_T2)];
nove_T=[I_mat, zero; zero, I_mat];

% !! úprava na pozici 2,2, puvodne F_c*V_2-F_s*V_1
F=[W_2,  W_2 ; W_1, W_1 ; 1i*V_2, -1i*V_2; 1i*V_1, -1i*V_1];

%[H_U,H_S,H_V]=svd(H);
%inv_H=H_V*inv(H_S)*H_U';
inv_F=inv(F);

a_b=inv_F*f_g;

a=a_b(1:2*n_sq, : );
b=a_b((2*n_sq+1):4*n_sq, : );

%[a_U,a_S,a_V]=svd(a);
%a_1=a_V*inv(a_S)*a_U';
a_1=inv(a);

a_1_X=a_1*X;
X_b_a_1_X=X*b*a_1_X; % rodis davá X na druhou stranu

f_g=F*[I_mat zero ; zero I_mat; X_b_a_1_X];

nove_T=nove_T*a_1_X;