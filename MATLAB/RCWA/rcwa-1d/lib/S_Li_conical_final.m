% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% Lifeng Li, "Formulation and comparison of two recursive matrix algorithms 
% for modeling layered diffraction gratings," J. Opt. Soc. Am. A 13, 1024-1035 (1996) 

% počáteční inicializace
W_2_11(:,:,number_of_layers)=zeros(2*number_of_orders,2*number_of_orders);
W_2_12(:,:,number_of_layers)=zeros(2*number_of_orders,2*number_of_orders);
W_2_21(:,:,number_of_layers)=zeros(2*number_of_orders,2*number_of_orders);
W_2_22(:,:,number_of_layers)=zeros(2*number_of_orders,2*number_of_orders);
X_1(:,:,number_of_layers)=zeros(number_of_orders,number_of_orders);
X_2(:,:,number_of_layers)=zeros(number_of_orders,number_of_orders);
b_uu=zeros(2*number_of_orders,2*number_of_orders);
b_ud=zeros(2*number_of_orders,2*number_of_orders);
b_du=zeros(2*number_of_orders,2*number_of_orders);
b_dd=zeros(2*number_of_orders,2*number_of_orders);
%a_uu(:,:,number_of_layers+1)=zeros(2*number_of_orders,2*number_of_orders);
a_ud(:,:,number_of_layers+1)=zeros(2*number_of_orders,2*number_of_orders);
%a_du(:,:,number_of_layers+1)=zeros(2*number_of_orders,2*number_of_orders);
a_dd(:,:,number_of_layers+1)=zeros(2*number_of_orders,2*number_of_orders);

s_1_l=zeros(4*number_of_orders,4*number_of_orders);

if number_of_layers>2
s_1(:,:,number_of_layers-1)=zeros(4*number_of_orders,4*number_of_orders);
end

%----------------------------------------------------------

W_3_11=[I zero; -1i*Y1 zero];
W_3_12=[sin(psi)*I zero ; 1i*n1*cos(theta)*sin(psi)*I zero];
W_3_21=[zero I ; zero -1i*Z1];
W_3_22=[-1i*n1*cos(psi)*I zero ; cos(psi)*cos(theta)*I zero];

for i=1:1:number_of_layers
    W_2_11(:,:,i)=[Vss(:,:,i) Vsp(:,:,i);-Wss(:,:,i) -Wsp(:,:,i)];
    W_2_12(:,:,i)=[Vss(:,:,i) Vsp(:,:,i);Wss(:,:,i) Wsp(:,:,i)];
    W_2_21(:,:,i)=[-Wps(:,:,i) -Wpp(:,:,i);Vps(:,:,i) Vpp(:,:,i)];
    W_2_22(:,:,i)=[Wps(:,:,i) Wpp(:,:,i);Vps(:,:,i) Vpp(:,:,i)];
end

W_1_11=[zero zero ; zero zero];
W_1_12=[I zero ;1i*Y3 zero];
W_1_21=[zero zero ; zero zero];
W_1_22=[zero I; zero 1i*Z3];

%---------------------------------------------------------

s_2=([W_3_11 -W_2_12(:,:,1) ; W_3_21 -W_2_22(:,:,1)])\[W_2_11(:,:,1) -W_3_12 ; W_2_21(:,:,1) -W_3_22];
for i=1:1:number_of_layers-1
    s_1(:,:,number_of_layers-i)=([W_2_11(:,:,i) -W_2_12(:,:,i+1) ; W_2_21(:,:,i) -W_2_22(:,:,i+1)])\[W_2_11(:,:,i+1) -W_2_12(:,:,i) ; W_2_21(:,:,i+1) -W_2_22(:,:,i)];
end
s_0=([W_2_11(:,:,number_of_layers) -W_1_12 ; W_2_21(:,:,number_of_layers) -W_1_22])\[W_1_11 -W_2_12(:,:,number_of_layers) ; W_1_21 -W_2_22(:,:,number_of_layers)];

for i=1:1:number_of_layers
    X_1(:,:,number_of_layers+1-i)=diag(exp(-k_0*diag(Q1(:,:,i))*layer_thickness(i)));
    X_2(:,:,number_of_layers+1-i)=diag(exp(-k_0*diag(Q2(:,:,i))*layer_thickness(i)));
end

%a_uu(:,:,1)=s_0(1:2*number_of_orders,1:2*number_of_orders);
a_ud(:,:,1)=s_0(1:2*number_of_orders,1+2*number_of_orders:4*number_of_orders);
%a_du(:,:,1)=s_0(1+2*number_of_orders:4*number_of_orders,1:2*number_of_orders);
a_dd(:,:,1)=s_0(2*number_of_orders+1:4*number_of_orders,2*number_of_orders+1:4*number_of_orders);

new_I=[I zero; zero I];

for i=1:1:number_of_layers
    if i==number_of_layers
        first_m=[I zero zero zero ; zero I zero zero; zero zero X_1(:,:,number_of_layers) zero; zero zero zero X_2(:,:,number_of_layers)];
        second_m=[X_1(:,:,number_of_layers) zero zero zero ; zero X_2(:,:,number_of_layers) zero zero; zero zero I zero; zero zero zero I];
        
        s_2_l=sparse(first_m)*s_2*sparse(second_m);
        
        b_uu=s_2_l(1:2*number_of_orders,1:2*number_of_orders);
        b_ud=s_2_l(1:2*number_of_orders,1+2*number_of_orders:4*number_of_orders);
        b_du=s_2_l(1+2*number_of_orders:4*number_of_orders,1:2*number_of_orders);
        b_dd=s_2_l(2*number_of_orders+1:4*number_of_orders,2*number_of_orders+1:4*number_of_orders);
        
    else
        first_m=[I zero zero zero ; zero I zero zero; zero zero X_1(:,:,i) zero; zero zero zero X_2(:,:,i)];
        second_m=[X_1(:,:,i) zero zero zero ; zero X_2(:,:,i) zero zero; zero zero I zero; zero zero zero I];
        
        s_1_l=diag(first_m)*conj((diag(second_m))').*s_1(:,:,i);
        %s_1_l=sparse(first_m)*s_1(:,:,i)*sparse(second_m);
        
        b_uu=s_1_l(1:2*number_of_orders,1:2*number_of_orders);    
        b_ud=s_1_l(1:2*number_of_orders,1+2*number_of_orders:4*number_of_orders);    
        b_du=s_1_l(1+2*number_of_orders:4*number_of_orders,1:2*number_of_orders);    
        b_dd=s_1_l(2*number_of_orders+1:4*number_of_orders,2*number_of_orders+1:4*number_of_orders);
    end
    
    %a_uu(:,:,i+1)=(b_uu/(new_I-a_ud(:,:,i)*b_du))*a_uu(:,:,i);
    a_ud(:,:,i+1)=b_ud+b_uu*(a_ud(:,:,i)/(new_I-b_du*a_ud(:,:,i)))*b_dd;
    %a_du(:,:,i+1)=a_du(:,:,i)+a_dd(:,:,i)*(b_du/(new_I-a_ud(:,:,i)*b_du))*a_uu(:,:,i);
    a_dd(:,:,i+1)=(a_dd(:,:,i)/(new_I-b_du*a_ud(:,:,i)))*b_dd;
    
end

%---------------------------------------------------------

T_dd=a_dd(:,:,number_of_layers+1);
R_ud=a_ud(:,:,number_of_layers+1);

%---------------------------------------------------------
d_n1=zeros(2*number_of_orders,1);
d_n1((number_of_orders+1)/2,1)=1;
%---------------------------------------------------------
R=R_ud*d_n1;
T=T_dd*d_n1;

R_s=R(1:number_of_orders);
R_p=R(number_of_orders+1:2*number_of_orders);
T_s=T(1:number_of_orders);
T_p=T(number_of_orders+1:2*number_of_orders);

%---------------------------------------------------------------------

if plot_field_region_2==1
    S_field_n_layers_conical;
end

clear W_2_11 W_2_12 W_2_21 W_2_22 X_1 X_2 a_ud a_dd %a_du a_uu