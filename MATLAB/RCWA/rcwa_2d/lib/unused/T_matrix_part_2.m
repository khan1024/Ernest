% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

incident=zeros(1,4*n_sq);

incident(N_inc)=sin(psi);
incident(n_sq+N_inc)=1i*n_1*cos(theta)*sin(psi);
incident(2*n_sq+N_inc)=-1i*n_1*cos(psi);
incident(3*n_sq+N_inc)=cos(psi)*cos(theta);

transfer=[-I_mat zero; 1i*K_1_z zero; zero -I_mat; zero 1i*K_1_z/(n_1^2)];   
transfer( : , 2*n_sq+1:4*n_sq)=f_g;
    
R_T1=transfer\incident.';
R_s=R_T1(1:n_sq);
R_p=R_T1(n_sq+1:2*n_sq);
T1=R_T1(2*n_sq+1:size(R_T1));
    
T=nove_T*T1;
T_s=T(1:n_sq);    
T_p=T(n_sq+1:size(T));