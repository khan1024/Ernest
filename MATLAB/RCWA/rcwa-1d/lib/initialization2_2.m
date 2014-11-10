% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

    
A_l=zeros(number_of_orders,number_of_orders,number_of_layers);    
E_l=zeros(number_of_orders,number_of_orders,number_of_layers);       
   
    
if polarization==0 
    A=zeros(number_of_orders,number_of_orders,number_of_layers);        
    W1=zeros(number_of_orders,number_of_orders,number_of_layers);        
    Q1=zeros(number_of_orders,number_of_orders,number_of_layers);        
    W2=zeros(number_of_orders,number_of_orders,number_of_layers);        
    Q2=zeros(number_of_orders,number_of_orders,number_of_layers);        
    V11=zeros(number_of_orders,number_of_orders,number_of_layers);        
    V12=zeros(number_of_orders,number_of_orders,number_of_layers);        
    V21=zeros(number_of_orders,number_of_orders,number_of_layers);        
    V22=zeros(number_of_orders,number_of_orders,number_of_layers);        
else
    Q=zeros(number_of_orders,number_of_orders,number_of_layers);
    V=zeros(number_of_orders,number_of_orders,number_of_layers);    
    W=zeros(number_of_orders,number_of_orders,number_of_layers);
end

K_x=zeros(number_of_orders,number_of_orders);
K_x_2=zeros(number_of_orders,number_of_orders);
Q_2=zeros(number_of_orders,number_of_orders);
Y1=zeros(number_of_orders,number_of_orders);
Y3=zeros(number_of_orders,number_of_orders);
Z1=zeros(number_of_orders,number_of_orders);
Z3=zeros(number_of_orders,number_of_orders);

D_R=zeros(number_of_orders,1);
D_T=zeros(number_of_orders,1);
k_1_z=zeros(number_of_orders,1);
k_1_z_k_0=zeros(number_of_orders,1);
k_3_z=zeros(number_of_orders,1);
k_3_z_k_0=zeros(number_of_orders,1);
k_x=zeros(number_of_orders,1);
k_x_2=zeros(number_of_orders,1);

k_y_k_0=zeros(number_of_orders,1);
k_y_2_I=zeros(number_of_orders,1);
phi_i=zeros(number_of_orders,1);

n=zeros(number_of_orders,1);

i_v=zeros(1,2*number_of_orders-1);
v=zeros(1,2*number_of_orders-1);

K=0;
M=0;
epsG=0;
epsr=0;
epsg=0;
k_0=0;
k_c=0;
m=0;
v_0=0;