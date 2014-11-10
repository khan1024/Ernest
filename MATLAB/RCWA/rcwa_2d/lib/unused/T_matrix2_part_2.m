% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

u_x=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
u_y=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);
u_z=-cos(psi)*sin(theta);
k_x00=n_1*sin(theta)*cos(phi); % vymazano k_0 díky vydělení
k_y00=n_1*sin(theta)*sin(phi); % vymazano k_0 díky vydělení
k_z00=n_1*cos(theta); % vymazano k_0 díky vydělení

incident=zeros(1,4*n_sq);

incident(N_inc)=u_x;
incident(n_sq+N_inc)=u_y;
incident(2*n_sq+N_inc)=k_y00*u_z-k_z00*u_y;
incident(3*n_sq+N_inc)=k_z00*u_x-k_x00*u_z;

transfer=[-I_mat, zero; zero, -I_mat; -K_x_T2*K_y_T2/(k_0*K_1_z_T2), -(K_y_T2.^2+K_1_z_T2.^2)/(k_0*K_1_z_T2); (K_x_T2.^2+K_1_z_T2.^2)/(k_0*K_1_z_T2), K_x_T2*K_y_T2/(k_0*K_1_z_T2)]; %oproti papíru opačná znaménka  
transfer( : , 2*n_sq+1:4*n_sq)=f_g;
    
R_T1=transfer\incident.';
R_x=R_T1(1:n_sq);
R_y=R_T1(n_sq+1:2*n_sq);
T1=R_T1(2*n_sq+1:size(R_T1));
    
T=nove_T*T1;
T_x=T(1:n_sq);    
T_y=T(n_sq+1:2*n_sq);

%R_z=(R_x.*k_x_mn' + R_y.*k_y_mn')./diag(k_1_z);
%T_z=-(T_x.*k_x_mn' + T_y.*k_y_mn')./diag(k_3_z);

R_z=(R_x.*diag(k_x_mn) + R_y.*diag(k_y_mn))./diag(k_1_z);
T_z=-(T_x.*diag(k_x_mn) + T_y.*diag(k_y_mn))./diag(k_3_z);