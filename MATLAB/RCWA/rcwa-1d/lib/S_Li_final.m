% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% Lifeng Li, "Formulation and comparison of two recursive matrix algorithms
% for modeling layered diffraction gratings," J. Opt. Soc. Am. A 13, 1024-1035 (1996) 

% počáteční inicializace
W_2_11(:,:,number_of_layers)=zeros(number_of_orders,number_of_orders);
W_2_12(:,:,number_of_layers)=zeros(number_of_orders,number_of_orders);
W_2_21(:,:,number_of_layers)=zeros(number_of_orders,number_of_orders);
W_2_22(:,:,number_of_layers)=zeros(number_of_orders,number_of_orders);
X1(:,:,number_of_layers)=zeros(number_of_orders,number_of_orders);
b_uu=zeros(number_of_orders,number_of_orders);
b_ud=zeros(number_of_orders,number_of_orders);
b_du=zeros(number_of_orders,number_of_orders);
b_dd=zeros(number_of_orders,number_of_orders);
%a_uu(:,:,number_of_layers+1)=zeros(number_of_orders,number_of_orders);
a_ud(:,:,number_of_layers+1)=zeros(number_of_orders,number_of_orders);
%a_du(:,:,number_of_layers+1)=zeros(number_of_orders,number_of_orders);
a_dd(:,:,number_of_layers+1)=zeros(number_of_orders,number_of_orders);

s_1_l=zeros(2*number_of_orders,2*number_of_orders);

if number_of_layers>2
s_1(:,:,number_of_layers-1)=zeros(2*number_of_orders,2*number_of_orders);
end

%----------------------------------------------------------

W_3_11=I;
W_3_12=I;

if polarization==1
    W_3_21=-1i*Z1;
    W_3_22=(1i*cos(theta)/n1)*I;
    W_1_22=1i*Z3;
elseif polarization==2
    W_3_21=-1i*Y1;
    W_3_22=(1i*n1*cos(theta))*I;
    W_1_22=1i*Y3;
end

for i=1:1:number_of_layers
    W_2_11(:,:,i)=W(:,:,i);
    W_2_12(:,:,i)=W(:,:,i);
    W_2_21(:,:,i)=-V(:,:,i);
    W_2_22(:,:,i)=V(:,:,i);
end

W_1_11=zero;% na nic
W_1_12=I;
W_1_21=I;% na nic

%---------------------------------------------------------

%s_2=inv([W_3_11 -W_2_12(:,:,1) ; W_3_21 -W_2_22(:,:,1)])*[W_2_11(:,:,1) -W_3_12 ; W_2_21(:,:,1) -W_3_22];
%for i=1:1:number_of_layers-1
%    s_1(:,:,number_of_layers-i)=inv([W_2_11(:,:,i) -W_2_12(:,:,i+1) ; W_2_21(:,:,i) -W_2_22(:,:,i+1)])*[W_2_11(:,:,i+1) -W_2_12(:,:,i) ; W_2_21(:,:,i+1) -W_2_22(:,:,i)];
%end
%s_0=inv([W_2_11(:,:,number_of_layers) -W_1_12 ; W_2_21(:,:,number_of_layers) -W_1_22])*[W_1_11 -W_2_12(:,:,number_of_layers) ; W_1_21 -W_2_22(:,:,number_of_layers)];

s_2=([W_3_11 -W_2_12(:,:,1) ; W_3_21 -W_2_22(:,:,1)])\[W_2_11(:,:,1) -W_3_12 ; W_2_21(:,:,1) -W_3_22];
for i=1:1:number_of_layers-1
    s_1(:,:,number_of_layers-i)=([W_2_11(:,:,i) -W_2_12(:,:,i+1) ; W_2_21(:,:,i) -W_2_22(:,:,i+1)])\[W_2_11(:,:,i+1) -W_2_12(:,:,i) ; W_2_21(:,:,i+1) -W_2_22(:,:,i)];
end
s_0=([W_2_11(:,:,number_of_layers) -W_1_12 ; W_2_21(:,:,number_of_layers) -W_1_22])\[W_1_11 -W_2_12(:,:,number_of_layers) ; W_1_21 -W_2_22(:,:,number_of_layers)];

for i=1:1:number_of_layers
    X1(:,:,number_of_layers+1-i)=diag(exp(-k_0*diag(Q(:,:,i))*layer_thickness(i)));
end

%a_uu(:,:,1)=s_0(1:number_of_orders,1:number_of_orders);
a_ud(:,:,1)=s_0(1:number_of_orders,1+number_of_orders:2*number_of_orders);
%a_du(:,:,1)=s_0(1+number_of_orders:2*number_of_orders,1:number_of_orders);
a_dd(:,:,1)=s_0(number_of_orders+1:2*number_of_orders,number_of_orders+1:2*number_of_orders);


for i=1:1:number_of_layers
    if i==number_of_layers
        s_2_l=[I zero;zero X1(:,:,number_of_layers)]*s_2*[X1(:,:,number_of_layers) zero;zero I];
        
        b_uu=s_2_l(1:number_of_orders,1:number_of_orders);
        b_ud=s_2_l(1:number_of_orders,1+number_of_orders:2*number_of_orders);
        b_du=s_2_l(1+number_of_orders:2*number_of_orders,1:number_of_orders);
        b_dd=s_2_l(number_of_orders+1:2*number_of_orders,number_of_orders+1:2*number_of_orders);
    else
        %s_1_l=[I zero;zero X1(:,:,i)]*s_1(:,:,i)*[X1(:,:,i) zero;zero I];
        s_1_l=diag([I zero;zero X1(:,:,i)])*conj((diag([X1(:,:,i) zero;zero I]))').*s_1(:,:,i);
        
        b_uu=s_1_l(1:number_of_orders,1:number_of_orders);    
        b_ud=s_1_l(1:number_of_orders,1+number_of_orders:2*number_of_orders);    
        b_du=s_1_l(1+number_of_orders:2*number_of_orders,1:number_of_orders);    
        b_dd=s_1_l(number_of_orders+1:2*number_of_orders,number_of_orders+1:2*number_of_orders);
    end
    
    %a_uu(:,:,i+1)=(b_uu/(I-a_ud(:,:,i)*b_du))*a_uu(:,:,i);
    a_ud(:,:,i+1)=b_ud+b_uu*(a_ud(:,:,i)/(I-b_du*a_ud(:,:,i)))*b_dd;
    %a_du(:,:,i+1)=a_du(:,:,i)+a_dd(:,:,i)*(b_du/(I-a_ud(:,:,i)*b_du))*a_uu(:,:,i);
    a_dd(:,:,i+1)=(a_dd(:,:,i)/(I-b_du*a_ud(:,:,i)))*b_dd;
    
end

%---------------------------------------------------------

T_dd=a_dd(:,:,number_of_layers+1);
R_ud=a_ud(:,:,number_of_layers+1);

%---------------------------------------------------------
d_n1=zeros(number_of_orders,1);
d_n1((number_of_orders+1)/2,1)=1;
%---------------------------------------------------------

R=R_ud*d_n1;
T=T_dd*d_n1;

if plot_field_region_2==1
    S_field_n_layers;
end