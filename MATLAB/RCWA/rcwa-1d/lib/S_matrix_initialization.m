% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

if change_matrix_base==1
    if polarization==0
        Vss=zeros(number_of_orders,number_of_orders,number_of_layers);
        Wss=zeros(number_of_orders,number_of_orders,number_of_layers);
        Vsp=zeros(number_of_orders,number_of_orders,number_of_layers);
        Wsp=zeros(number_of_orders,number_of_orders,number_of_layers);
        Wpp=zeros(number_of_orders,number_of_orders,number_of_layers);
        Vpp=zeros(number_of_orders,number_of_orders,number_of_layers);
        Wps=zeros(number_of_orders,number_of_orders,number_of_layers);
        Vps=zeros(number_of_orders,number_of_orders,number_of_layers); 
    else
    end
end


F_c=diag(cos(phi_i));    
F_s=diag(sin(phi_i));    
    
k_1_z_k_0=k_1_z/k_0;    
Y1=diag(k_1_z_k_0);    
Z1=diag(k_1_z_k_0/((n1)^2));
    
k_3_z_k_0=k_3_z/k_0;    
Y3=diag(k_3_z_k_0);    
Z3=diag(k_3_z_k_0/((n3)^2));

K_y_n=k_y_k_0*I;
inv_Y1=diag(1./k_1_z_k_0);
inv_Y3=diag(1./k_3_z_k_0);

if matrix_algorithm==1 && change_matrix_base==2
    WW=zeros(2*number_of_orders,2*number_of_orders,number_of_layers);
    VV=zeros(2*number_of_orders,2*number_of_orders,number_of_layers);
end