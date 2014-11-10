% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% Lifeng Li, "Note on the S-matrix propagation algorithm," J. Opt. Soc. Am. A 20, 655-660 (2003)

phi_minus=zeros(2*number_of_orders,2*number_of_orders,number_of_layers+1);
W_2_1=zeros(2*number_of_orders,2*number_of_orders,number_of_layers);
W_2_2=zeros(2*number_of_orders,2*number_of_orders,number_of_layers);

%R_ud_tilde=zeros(2*number_of_orders,2*number_of_orders,number_of_layers+1);
%T_dd_tilde=zeros(2*number_of_orders,2*number_of_orders,number_of_layers+1);
%T_uu_tilde=zeros(2*number_of_orders,2*number_of_orders,number_of_layers+1);

%R_ud=zeros(2*number_of_orders,2*number_of_orders,number_of_layers+1);
%R_du=zeros(2*number_of_orders,2*number_of_orders,number_of_layers+1);
%T_uu=zeros(2*number_of_orders,2*number_of_orders,number_of_layers+1);
%T_dd=zeros(2*number_of_orders,2*number_of_orders,number_of_layers+1);

f_g_1=zeros(number_of_orders,number_of_orders);
f_g_2=zeros(number_of_orders,number_of_orders);
f_g_3=zeros(number_of_orders,number_of_orders);
f_g_4=zeros(number_of_orders,number_of_orders);

big_M_1=zeros(number_of_orders,number_of_orders);
big_M_2=zeros(number_of_orders,number_of_orders);
big_M_3=zeros(number_of_orders,number_of_orders);
big_M_4=zeros(number_of_orders,number_of_orders);

W_3_11=zeros(2*number_of_orders,2*number_of_orders);
W_3_12=zeros(2*number_of_orders,2*number_of_orders);
W_3_21=zeros(2*number_of_orders,2*number_of_orders);
W_3_22=zeros(2*number_of_orders,2*number_of_orders);

W_1_11=zeros(2*number_of_orders,2*number_of_orders);
W_1_12=zeros(2*number_of_orders,2*number_of_orders);
W_1_21=zeros(2*number_of_orders,2*number_of_orders);
W_1_22=zeros(2*number_of_orders,2*number_of_orders);
%------------------------------------------------------

f_g_1=-1i*K_y_n*inv_Y3*K_x;
f_g_2=-1i*(Y3+K_y_n*inv_Y3*K_y_n);
f_g_3=1i*(Y3+K_x*inv_Y3*K_x);
f_g_4=1i*K_x*inv_Y3*K_y_n;

big_M_1=1i*K_y_n*inv_Y1*K_x;
big_M_2=1i*(Y1+K_y_n*inv_Y1*K_y_n);
big_M_3=-1i*(Y1+K_x*inv_Y1*K_x);
big_M_4=-1i*K_x*inv_Y1*K_y_n;

incident_1=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
incident_2=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);
incident_3=-1i*n1*(cos(psi)*sin(phi)+sin(psi)*cos(theta)*cos(phi));
incident_4=1i*n1*(cos(psi)*cos(phi)-sin(psi)*cos(theta)*sin(phi));

W_3_11=[I, zero; zero, I];
W_3_12=[incident_1*I, zero; incident_2*I, zero];
W_3_21=[big_M_1 big_M_2; big_M_3 big_M_4];
W_3_22=[incident_3*I, zero; incident_4*I, zero];

S_express=1;

W_1_11=[zero, zero; zero, zero];% na nic
W_1_12=[I, zero; zero, I];
W_1_21=[zero, zero; zero, zero];% na nic
W_1_22=[f_g_1 f_g_2; f_g_3 f_g_4];

R_ud=[zero, zero; zero, zero];
R_du=[zero, zero; zero, zero];
T_uu=[I, zero; zero, I];
T_dd=[I, zero; zero, I];

for i=number_of_layers+1:-1:1
    if i==1
        phi_minus(:,:,i)=[I zero; zero I];
    else
        index_phi=number_of_layers-i+2;
        index_W_V=number_of_layers-i+2;
        phi_minus(:,:,i)=[diag(exp(-k_0*diag(Q1(:,:,index_phi))*layer_thickness(index_phi))), zero; zero, diag(exp(-k_0*diag(Q2(:,:,index_phi))*layer_thickness(index_phi)))];
        W_2_1(:,:,i)=WW(:,:,index_W_V);
        W_2_2(:,:,i)=VV(:,:,index_W_V);
    end
end

for i=1:1:number_of_layers+1
    %R_ud_tilde=phi_minus(:,:,i)*R_ud*phi_minus(:,:,i);
    %T_dd_tilde=T_dd*phi_minus(:,:,i);
    
    R_ud_tilde=diag(phi_minus(:,:,i))*conj((diag(phi_minus(:,:,i)))').*R_ud;
    T_dd_tilde=T_dd*sparse(phi_minus(:,:,i));
    
    if S_express==1
        % nic
    else
        T_uu_tilde=phi_minus(:,:,i)*T_uu;    
    end
    
    if i==1 || i==number_of_layers+1
        if i==1
            if S_express==1
                T_uu_tilde=phi_minus(:,:,i)*T_uu;
            end
            
            X_1=[W_1_11*T_uu_tilde; W_1_21*T_uu_tilde];
            X_2=[-W_2_1(:,:,i+1); -W_2_2(:,:,i+1)];
            
            Z=[W_2_1(:,:,i+1), -W_1_11*R_ud_tilde-W_1_12; -W_2_2(:,:,i+1), -W_1_21*R_ud_tilde-W_1_22];
        elseif i==number_of_layers+1
            if S_express==1
                T_uu_tilde=phi_minus(:,:,i)*T_uu;
            end
            
            X_1=[W_2_1(:,:,i)*T_uu_tilde; -W_2_2(:,:,i)*T_uu_tilde];
            X_2=[-W_3_12; -W_3_22];
        
            Z=[W_3_11, -W_2_1(:,:,i)*R_ud_tilde-W_2_1(:,:,i); W_3_21, W_2_2(:,:,i)*R_ud_tilde-W_2_2(:,:,i)];
        end
        
        %matice1=inv(Z)*X_1;
        %matice2=inv(Z)*X_2;
        
        matice1=Z\X_1;
        matice2=Z\X_2;
        
        R_ud=matice2(1:2*number_of_orders,:);
        T_dd=T_dd_tilde*matice2((2*number_of_orders+1):4*number_of_orders,:);
        
        if S_express==1
        else
            T_uu=matice1(1:2*number_of_orders,:);
            R_du=R_du+T_dd_tilde*matice1((2*number_of_orders+1):4*number_of_orders,:);
        end
        
    else

        %Q_01=inv(W_2_1(:,:,i+1))*W_2_1(:,:,i);    
        %Q_02=inv(W_2_2(:,:,i+1))*W_2_2(:,:,i);
        
        Q_01=W_2_1(:,:,i+1)\W_2_1(:,:,i);    
        Q_02=W_2_2(:,:,i+1)\W_2_2(:,:,i);

        F=Q_01*([I zero; zero I]+R_ud_tilde);
        G=Q_02*([I zero; zero I]-R_ud_tilde);
                
        tau=inv(F+G);
    
        R_ud=[I zero; zero I]-2*G*tau;
        T_dd=2*T_dd_tilde*tau;
    
        if S_express==1
            if i==number_of_layers
                T_uu=(F*tau*Q_02+G*tau*Q_01)*T_uu_tilde;
                R_du=R_du+T_dd_tilde*tau*(Q_02-Q_01)*T_uu_tilde;    
            end
        else
            T_uu=(F*tau*Q_02+G*tau*Q_01)*T_uu_tilde;
            R_du=R_du+T_dd_tilde*tau*(Q_02-Q_01)*T_uu_tilde;    
        end

    end
end

%---------------------------------------------------------
d_n1=zeros(2*number_of_orders,1);
d_n1((number_of_orders+1)/2,1)=1;
%---------------------------------------------------------

u_n1=R_ud*d_n1;
d_0=T_dd*d_n1;

R_x=u_n1(1:number_of_orders);
R_y=u_n1(number_of_orders+1:2*number_of_orders);
T_x=d_0(1:number_of_orders);
T_y=d_0(number_of_orders+1:2*number_of_orders);
R_z=(R_x.*k_x + R_y.*k_y)./k_1_z;
T_z=-(T_x.*k_x + T_y.*k_y)./k_3_z;