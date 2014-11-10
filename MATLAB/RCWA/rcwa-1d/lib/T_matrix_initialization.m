% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% M. G. Moharam, Drew A. Pommet, Eric B. Grann, and T. K. Gaylord, 
% "Stable implementation of the rigorous coupled-wave analysis for surface-relief
% gratings: enhanced transmittance matrix approach," J. Opt. Soc. Am. A 12, 1077-1086 (1995) 

if polarization==0        
    X=zeros(2*number_of_orders,2*number_of_orders,number_of_layers);        
    a_1_X=zeros(2*number_of_orders,2*number_of_orders,number_of_layers);        
    b_a_1_X=zeros(2*number_of_orders,2*number_of_orders,number_of_layers); 
    Vss=zeros(number_of_orders,number_of_orders,number_of_layers);        
    Wss=zeros(number_of_orders,number_of_orders,number_of_layers);        
    Vsp=zeros(number_of_orders,number_of_orders,number_of_layers);        
    Wsp=zeros(number_of_orders,number_of_orders,number_of_layers);        
    Wpp=zeros(number_of_orders,number_of_orders,number_of_layers);        
    Vpp=zeros(number_of_orders,number_of_orders,number_of_layers);        
    Wps=zeros(number_of_orders,number_of_orders,number_of_layers);        
    Vps=zeros(number_of_orders,number_of_orders,number_of_layers);        
    X1=zeros(number_of_orders,number_of_orders,number_of_layers);        
    X2=zeros(number_of_orders,number_of_orders,number_of_layers);
    if change_matrix_base==2
        new_W=zeros(2*number_of_orders,2*number_of_orders,number_of_layers);  
        new_V=zeros(2*number_of_orders,2*number_of_orders,number_of_layers);
    end
else
    %if S_matice~1            
        X=zeros(number_of_orders,number_of_orders,number_of_layers);            
        a_1_X=zeros(number_of_orders,number_of_orders,number_of_layers);            
        b_a_1_X=zeros(number_of_orders,number_of_orders,number_of_layers);        
    %end
end


if polarization==0        
    a=zeros(2*number_of_orders,2*number_of_orders);
    a_1=zeros(2*number_of_orders,2*number_of_orders);
    b=zeros(2*number_of_orders,2*number_of_orders);
    R_T1=zeros(4*number_of_orders,1);
    T1=zeros(2*number_of_orders,1);
    a_b=zeros(4*number_of_orders,2*number_of_orders);
    f_g=zeros(4*number_of_orders,2*number_of_orders);
    X_b_a_1_x=zeros(2*number_of_orders,2*number_of_orders);
else
    a=zeros(number_of_orders,number_of_orders);        
    a_1=zeros(number_of_orders,number_of_orders);        
    b=zeros(number_of_orders,number_of_orders);        
    R=zeros(1,number_of_orders);        
    R_T1=zeros(1,number_of_orders);        
    T1=zeros(1,number_of_orders);        
    T_1=zeros(1,number_of_orders);        
    a_b=zeros(2*number_of_orders,number_of_orders);        
    f_g=zeros(2*number_of_orders,number_of_orders);        
    X_b_a_1_x=zeros(number_of_orders,number_of_orders);
    if change_matrix_base==2
        f_g_1=zeros(number_of_orders,number_of_orders);
        f_g_2=zeros(number_of_orders,number_of_orders);
        f_g_3=zeros(number_of_orders,number_of_orders);
        f_g_4=zeros(number_of_orders,number_of_orders);
    end
end

    
%------------------------------------------------------------------------
zero=zeros(number_of_orders);
if polarization==0
    new_T=[I zero; zero I];
else
    new_T=I;
end

F_c=diag(cos(phi_i));
F_s=diag(sin(phi_i));

% definition of $f_{L+1}$ and $g_{L+1}$

k_1_z_k_0=k_1_z/k_0;
Y1=diag(k_1_z_k_0);
Z1=diag(k_1_z_k_0/((n1)^2));

k_3_z_k_0=k_3_z/k_0;
Y3=diag(k_3_z_k_0);
Z3=diag(k_3_z_k_0/((n3)^2));

if polarization==0
    incident=zeros(1,4*number_of_orders);
else
    incident = zeros(1,2*number_of_orders);
end

if polarization==1 %TM
    incident(order_max+1)=1;
    incident(number_of_orders+order_max+1)=1i*cos(theta)/n1;
    f_g=[I;1i*Z3];
elseif polarization==2 %TE
    incident(order_max+1)=1;
    incident(number_of_orders+order_max+1)=1i*n1*cos(theta);
    f_g=[I;1i*Y3];
elseif polarization==0
    if change_matrix_base==1    
        incident(order_max+1)=sin(psi);    
        incident(number_of_orders+order_max+1)=1i*n1*cos(theta)*sin(psi);    
        incident(2*number_of_orders+order_max+1)=-1i*n1*cos(psi);    
        incident(3*number_of_orders+order_max+1)=cos(psi)*cos(theta);    
        f_g=[I zero; 1i*Y3 zero; zero I ; zero 1i*Z3];
    else
        inv_Y1=diag(1./k_1_z_k_0);
        inv_Y3=diag(1./k_3_z_k_0);
        
        incident(order_max+1)=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);    
        incident(number_of_orders+order_max+1)=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);    
        incident(2*number_of_orders+order_max+1)=-1i*n1*(cos(psi)*sin(phi)+sin(psi)*cos(theta)*cos(phi)); % změna znaménka    
        incident(3*number_of_orders+order_max+1)=1i*n1*(cos(psi)*cos(phi)-sin(psi)*cos(theta)*sin(phi));
        
        K_y_n=k_y_k_0*I;
        f_g_1=-1i*K_y_n*inv_Y3*K_x;    
        f_g_2=-1i*(Y3+K_y_n*inv_Y3*K_y_n);    
        f_g_3=1i*(Y3+K_x*inv_Y3*K_x);    
        f_g_4=1i*K_x*inv_Y3*K_y_n;    
        f_g=[I zero; zero I; f_g_1 f_g_2; f_g_3 f_g_4];    
        f_new=[I zero; zero I];    
        g_new=[f_g_1 f_g_2; f_g_3 f_g_4];
    end
end