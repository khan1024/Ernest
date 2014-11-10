% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% Na začátku vyřeším některé jednoduché výpočty
% Výpočet c_m(:,:,number_of_layers) --- z=d
c_m(:,:,number_of_layers)=s_0(1:2*number_of_orders,1+2*number_of_orders:4*number_of_orders)*(s_0(2*number_of_orders+1:4*number_of_orders,2*number_of_orders+1:4*number_of_orders)\T);


% Výpočet c_p(:,:,1) --- z=0

inv_s_n_11=inv(s_2(1:2*number_of_orders,1:2*number_of_orders));

u_n_n=inv_s_n_11*R-inv_s_n_11*s_2(1:2*number_of_orders,1+2*number_of_orders:4*number_of_orders)*d_n1;
c_p(:,:,1)=s_2(1+2*number_of_orders:4*number_of_orders,1:2*number_of_orders)*u_n_n+s_2(2*number_of_orders+1:4*number_of_orders,2*number_of_orders+1:4*number_of_orders)*d_n1;

d_n_plus_1_n(:,:,number_of_layers)=[X_1(:,:,number_of_layers) zero; zero X_2(:,:,number_of_layers)]*c_p(:,:,1); % výpočet do dalších částí 

%nepřehledný for-cyklus, nejlépe kouknou do poznámek
if number_of_layers>=2    
    for i=number_of_layers-1:-1:1
        c_m(:,:,number_of_layers-i)=a_ud(:,:,i+1)*d_n_plus_1_n(:,:,i+1); %využívám S-matice --- koeficient R_ud  ;    
    
        inv_s_2_11=inv(s_1(1:2*number_of_orders,1:2*number_of_orders,i));
    
        u_n_n=inv_s_2_11*c_m(:,:,number_of_layers-i)-inv_s_2_11*s_1(1:2*number_of_orders,1+2*number_of_orders:4*number_of_orders,i)*d_n_plus_1_n(:,:,i+1);    
        c_p(:,:,number_of_layers+1-i)=s_1(1+2*number_of_orders:4*number_of_orders,1:2*number_of_orders,i)*u_n_n+s_1(2*number_of_orders+1:4*number_of_orders,2*number_of_orders+1:4*number_of_orders,i)*d_n_plus_1_n(:,:,i+1);
       
        d_n_plus_1_n(:,:,i)=[X_1(:,:,i) zero; zero X_2(:,:,i)]*c_p(:,:,number_of_layers+1-i); % výpočet do dalších částí 
    end
end
