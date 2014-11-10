% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% Lifeng Li, "Note on the S-matrix propagation algorithm," J. Opt. Soc. Am. A 20, 655-660 (2003)

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

S_express=1;

% inicialization
W_1_11=zero;% na nic
W_1_12=I;
W_1_21=I;% na nic

phi_minus=zeros(number_of_orders,number_of_orders,number_of_layers+1);
W_2_1=zeros(number_of_orders,number_of_orders,number_of_layers);
W_2_2=zeros(number_of_orders,number_of_orders,number_of_layers);

R_du_tilde=zeros(number_of_orders,number_of_orders);
R_ud_tilde=zeros(number_of_orders,number_of_orders);
T_dd_tilde=zeros(number_of_orders,number_of_orders);
T_uu_tilde=zeros(number_of_orders,number_of_orders);

if plot_field_region_2==1
    s_11=zeros(number_of_orders,number_of_orders,number_of_layers+1);
    s_12=zeros(number_of_orders,number_of_orders,number_of_layers+1);
    s_21=zeros(number_of_orders,number_of_orders,number_of_layers+1);
    s_22=zeros(number_of_orders,number_of_orders,number_of_layers+1);
    for_inv_s_22=zeros(number_of_orders,number_of_orders);
end

R_ud=zeros(number_of_orders,number_of_orders,number_of_layers+1);
R_du=zeros(number_of_orders,number_of_orders,number_of_layers+1);
T_uu=zeros(number_of_orders,number_of_orders,number_of_layers+1);
T_dd=zeros(number_of_orders,number_of_orders,number_of_layers+1);

R_ud(:,:,1)=zero;
R_du(:,:,1)=zero;
T_uu(:,:,1)=I;
T_dd(:,:,1)=I;

for i=number_of_layers+1:-1:1
    if i==1
        phi_minus(:,:,i)=I;
    else
        index_phi=number_of_layers-i+2;
        index_W_V=number_of_layers-i+2;
        phi_minus(:,:,i)=diag(exp(-k_0*diag(Q(:,:,index_phi))*layer_thickness(index_phi)));
        W_2_1(:,:,i)=W(:,:,index_W_V);
        W_2_2(:,:,i)=V(:,:,index_W_V);
    end
end

for i=1:1:number_of_layers+1
    R_ud_tilde=phi_minus(:,:,i)*R_ud(:,:,i)*phi_minus(:,:,i);
    T_dd_tilde=T_dd(:,:,i)*phi_minus(:,:,i);
        
    if S_express==1
        % skip
    else
        T_uu_tilde=phi_minus(:,:,i)*T_uu(:,:,i);    
    end
    
    if i==1 || i==number_of_layers+1
        if i==1
            if S_express==1
                T_uu_tilde=phi_minus(:,:,i)*T_uu(:,:,i);
            end
            
            X_1=[W_1_11*T_uu_tilde; W_1_21*T_uu_tilde];
            X_2=[-W_2_1(:,:,i+1); -W_2_2(:,:,i+1)];
            
            Z=[W_2_1(:,:,i+1), -W_1_11*R_ud_tilde-W_1_12; -W_2_2(:,:,i+1), -W_1_21*R_ud_tilde-W_1_22];
        elseif i==number_of_layers+1
            if S_express==1
                T_uu_tilde=phi_minus(:,:,i)*T_uu(:,:,i);
            end
            
            X_1=[W_2_1(:,:,i)*T_uu_tilde; -W_2_2(:,:,i)*T_uu_tilde];
            X_2=[-W_3_12; -W_3_22];
        
            Z=[W_3_11, -W_2_1(:,:,i)*R_ud_tilde-W_2_1(:,:,i); W_3_21, W_2_2(:,:,i)*R_ud_tilde-W_2_2(:,:,i)];
        end
        
        inv_Z=inv(Z);
        matrix_1=inv_Z*X_1;
        matrix_2=inv_Z*X_2;
        
        R_ud(:,:,i+1)=matrix_2(1:number_of_orders,:);
        T_dd(:,:,i+1)=T_dd_tilde*matrix_2((number_of_orders+1):2*number_of_orders,:);
        
        if S_express==1
        else
            T_uu(:,:,i+1)=matrix_1(1:number_of_orders,:);
            R_du(:,:,i+1)=R_du(:,:,i)+T_dd_tilde*matrix_1((number_of_orders+1):2*number_of_orders,:);
        end
        
    else

        %Q_01=inv(W_2_1(:,:,i+1))*W_2_1(:,:,i);    
        %Q_02=inv(W_2_2(:,:,i+1))*W_2_2(:,:,i);
        
        Q_01=W_2_1(:,:,i+1)\W_2_1(:,:,i);    
        Q_02=W_2_2(:,:,i+1)\W_2_2(:,:,i);

        F=Q_01*(I+R_ud_tilde);
        G=Q_02*(I-R_ud_tilde);
                
        tau=inv(F+G);
    
        R_ud(:,:,i+1)=I-2*G*tau;
        T_dd(:,:,i+1)=2*T_dd_tilde*tau;
        
        if plot_field_region_2==1        
            t_11=1/2*(Q_01+Q_02);  
            t_12=1/2*(Q_01-Q_02);
        
            if i==1
                for_inv_s_22=t_11;
            end
            
            inv_t_22=inv(t_11);
            inv_t_22_t_12=inv_t_22*t_12;
        
            s_11(:,:,i)=t_11-t_12*inv_t_22_t_12; % t_21=t_12
            s_12(:,:,i)=t_12*inv_t_22;
            s_21(:,:,i)=-inv_t_22_t_12;  % t_21=t_12
            s_22(:,:,i)=inv_t_22;    
    
        end
end

clear t_11 t_12 inv_t_22 inv_t_22_t_12
    
        if S_express==1
            if i==number_of_layers
                T_uu(:,:,i+1)=(F*tau*Q_02+G*tau*Q_01)*T_uu_tilde;
                R_du(:,:,i+1)=R_du(:,:,i)+T_dd_tilde*tau*(Q_02-Q_01)*T_uu_tilde;    
            end
        else
            T_uu(:,:,i+1)=(F*tau*Q_02+G*tau*Q_01)*T_uu_tilde;
            R_du(:,:,i+1)=R_du(:,:,i)+T_dd_tilde*tau*(Q_02-Q_01)*T_uu_tilde;    
        end

    end
end

%---------------------------------------------------------
d_n1=zeros(number_of_orders,1);
d_n1((number_of_orders+1)/2,1)=1;
%---------------------------------------------------------

R=R_ud(:,:,number_of_layers+2)*d_n1;
T=T_dd(:,:,number_of_layers+2)*d_n1;

if plot_field_region_2==1
    c_m(:,:,number_of_layers)=s_12(:,:,1)*for_inv_s_22*T; % u_1_0
    c_p(:,:,1)=s_21(:,:,end)*(s_11(:,:,end)\(R-s_12(:,:,end)*d_n1))+s_22(:,:,end)*d_n1;
    d_n_plus_1_n(:,:,number_of_layers)=phi_minus(:,:,end)*c_p(:,:,1); % výpočet do dalších částí 

    if number_of_layers>=2
        for i=number_of_layers-1:-1:1
            c_m(:,:,number_of_layers-i)=R_ud(:,:,i+2)*d_n_plus_1_n(:,:,i+1);
            c_p(:,:,number_of_layers+1-i)=s_21(:,:,i+1)*(s_11(:,:,i+1)\(c_m(:,:,number_of_layers-i)-s_12(:,:,i+1)*d_n_plus_1_n(:,:,i+1)))+s_22(:,:,i+1)*d_n_plus_1_n(:,:,i+1);
            d_n_plus_1_n(:,:,i)=phi_minus(:,:,i+1)*c_p(:,:,number_of_layers+1-i); % výpočet do dalších částí 
        end
    end
end