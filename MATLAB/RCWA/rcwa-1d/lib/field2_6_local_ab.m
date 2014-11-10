% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% P. Lalanne and M.P. Jurek, J. Mod. Opt. 45, 1357-1374 (1998)

mu=1.2566370614E-6;
epsilon=8.854187818E-12;

thickness=thickness_total*1e-6; %patch
if plot_field_region_2==1
    if grating==10
        x_min=data(1,2)*1e-6;
        x_max=data(1,2)*1e-6+number_of_plotted_period*perioda;
    else
        x_min=0;
        x_max=number_of_plotted_period*perioda;% kolik period mřížky si necháme vykreslit
    end
    x_res_=perioda/resolution_x;% resolution in x-direction
    z_res_=thickness/resolution_z;% rozlišení z
    z_22=thickness/resolution_z:thickness/resolution_z:thickness;
    %D Nadefinujeme si počáteční hodnoty
    if matrix_algorithm==2
        T_l=T1;
    end
    D_l_pred=0;%na počátku je tloušťka předešlé vrstvičky nulová
    
    x=x_min:x_res_:x_max;
    exp_xx=zeros(length(k_x),length(x));
    
    for i=1:length(x)
        exp_xx(:,i)=exp(-1i*k_x*x(i));
    end
    
    
    %D Při výpočtu postupujeme od \quotation{zhora} směrem %\quotation{dolů}.
    for l=1:1:number_of_layers 
        D_l=D_l_pred+layer_thickness(l);%celková tloušťka do (včetně) l-té vrstvičky  
        z=(D_l_pred+z_res_):z_res_:D_l; %pohybuju se uprostřed l-té vrstvičky
        absorption_in_layer_thickness(l)=length(z);
        %x=x_min:x_res_:x_max;
        z_plus=z-D_l_pred; % exponent do soustavy vázaných rovnic
        z_minus=z-D_l; % exponent do soustavy vázaných rovnic
        
        %epsilon_matrix_local=zeros(length(x),length(z));
        %D Nadefinujeme si potřebné matice
        if polarization==0 %conical
            if conical_field_component==1
                H_x=zeros(length(x),length(z));
            elseif conical_field_component==2
                H_y=zeros(length(x),length(z));
            elseif conical_field_component==3
                H_z=zeros(length(x),length(z));
            elseif conical_field_component==4
                D_x=zeros(length(x),length(z));
                E_x=zeros(length(x),length(z));
            elseif conical_field_component==5
                E_y=zeros(length(x),length(z));
            elseif conical_field_component==6
                E_z=zeros(length(x),length(z));
            end
        elseif polarization==1 %TM
            if planar_field_component==1
                D_x=zeros(length(x),length(z));
                E_x=zeros(length(x),length(z));
            elseif planar_field_component==2
                H_y=zeros(length(x),length(z));
            elseif planar_field_component==3
                E_z=zeros(length(x),length(z));
            end
        elseif polarization==2 %TE
            if planar_field_component==1
                H_x=zeros(length(x),length(z));
            elseif planar_field_component==2
                E_y=zeros(length(x),length(z));
            elseif planar_field_component==3
                H_z=zeros(length(x),length(z));
            end
        end
        
        %D V textu jsme odvodili následující vzoreček pro výpočet konstant
        %D \type{c_plus} --- $c^+_{l,m}$ a \type{c_minus} --- $c^-_{l,m}$
        %D \placeformula[-]
        %D \startformula
        %D \pmatrix{\mat{c}^{+}_{l} \cr \mat{c}^{-}_{l}}=
        %D \pmatrix{\mat{I} \cr \mat{b}_l\mat{a}^{-1}_l\mat{X}_l}\mat{T}_l.
        %D \stopformula
        
        if matrix_algorithm==1 % S-matrix
            if polarization==0
                c_1_plus=c_p(1:number_of_orders,:,l);            
                c_2_plus=c_p(number_of_orders+1:2*number_of_orders,:,l);            
                c_1_minus=c_m(1:number_of_orders,:,l);            
                c_2_minus=c_m(number_of_orders+1:2*number_of_orders,:,l);
            else
                c_plus=c_p(:,:,l);
                c_minus=c_m(:,:,l);
            end
        elseif matrix_algorithm==2 % T-matrix
            if polarization==0
                pok=[I zero ; zero I ; b_a_1_X(:,:,l)];
                c_plus_minus=pok*T_l;
                c_1_plus=c_plus_minus(1:number_of_orders,:);
                c_2_plus=c_plus_minus(number_of_orders+1:2*number_of_orders ,:);
                c_1_minus=c_plus_minus(2*number_of_orders+1:3*number_of_orders,:);
                c_2_minus=c_plus_minus(3*number_of_orders+1:size(c_plus_minus),:);
            else
                c_plus_minus=[I ; b_a_1_X(:,:,l)]*T_l;
                c_plus=c_plus_minus(1:number_of_orders,:);
                c_minus=c_plus_minus(number_of_orders+1:size(c_plus_minus),:);
            end
        end
        
        if polarization==0
            for z_l=1:1:length(z)
                
                c_1_plus_exp=c_1_plus.*exp(-k_0*diag(Q1(:,:,l))*z_plus(z_l));
                c_2_plus_exp=c_2_plus.*exp(-k_0*diag(Q2(:,:,l))*z_plus(z_l));
                c_1_minus_exp=c_1_minus.*exp(k_0*diag(Q1(:,:,l))*z_minus(z_l));
                c_2_minus_exp=c_2_minus.*exp(k_0*diag(Q2(:,:,l))*z_minus(z_l));
                
                if conical_field_component==1
                    U_x_harmonics=W1(:,:,l)*(-c_1_plus_exp+c_1_minus_exp);
                elseif conical_field_component==2
                    U_y_harmonics=V21(:,:,l)*(-c_1_plus_exp+c_1_minus_exp)+V22(:,:,l)*(-c_2_plus_exp+c_2_minus_exp);
                elseif conical_field_component==3
                    S_x_harmonics=W2(:,:,l)*(c_2_plus_exp+c_2_minus_exp);
                    S_y_harmonics=V11(:,:,l)*(c_1_plus_exp+c_1_minus_exp)+V12(:,:,l)*(c_2_plus_exp+c_2_minus_exp);
                    U_z_harmonics=1i*(-K_x*S_y_harmonics+k_y/k_0*S_x_harmonics);
                elseif conical_field_component==4
                    S_x_harmonics=W2(:,:,l)*(c_2_plus_exp+c_2_minus_exp);
                    G_x_harmonics=(1/377)*sqrt(epsilon*mu)*(V11(:,:,l)*Q1(:,:,l)*(+c_1_plus_exp+c_1_minus_exp)+V22(:,:,l)*Q2(:,:,l)*(c_2_plus_exp+c_2_minus_exp)); % ?V11 nebo V12
                elseif conical_field_component==5
                    S_y_harmonics=V11(:,:,l)*(c_1_plus_exp+c_1_minus_exp)+V12(:,:,l)*(c_2_plus_exp+c_2_minus_exp);
                elseif conical_field_component==6
                    U_x_harmonics=W1(:,:,l)*(-c_1_plus_exp+c_1_minus_exp);
                    U_y_harmonics=V21(:,:,l)*(-c_1_plus_exp+c_1_minus_exp)+V22(:,:,l)*(-c_2_plus_exp+c_2_minus_exp);
                    S_z_harmonics=1i*(E_l(:,:,l)\(K_x*U_y_harmonics-k_y/k_0*U_x_harmonics));    
                end
                    
                if conical_field_component==1                        
                    H_x(:,z_l)=-1i*sqrt(epsilon/mu)*U_x_harmonics'*exp_xx;                    
                elseif conical_field_component==2                        
                    H_y(:,z_l)=-1i*sqrt(epsilon/mu)*U_y_harmonics'*exp_xx; % jiné normování než u planární difrakce                    
                elseif conical_field_component==3                        
                    H_z(:,z_l)=-1i*sqrt(epsilon/mu)*U_z_harmonics'*exp_xx;                    
                elseif conical_field_component==4                        
                    if phi0==0
                        D_x(:,z_l)=G_x_harmonics'*exp_xx;
                        for x_l=1:1:length(x)
                            field2_correct; % správně vypočítané E_x
                            field2_correct_local_ab;
                        end
                    else
                        E_x(:,z_l)=S_x_harmonics'*exp_xx;
                        for x_l=1:1:length(x)
                            field2_correct_local_ab;
                        end
                    end
                elseif conical_field_component==5                        
                    E_y(:,z_l)=S_y_harmonics'*exp_xx;                    
                elseif conical_field_component==6                        
                    E_z(:,z_l)=S_z_harmonics'*exp_xx;                    
                end  
            end
            
            if conical_field_component==1
                field_H=H_x;
            elseif conical_field_component==2
                field_H=H_y;
            elseif conical_field_component==3
                field_H=H_z;
            elseif conical_field_component==4
                field_E=E_x;
            elseif conical_field_component==5
                field_E=E_y;
            elseif conical_field_component==6
                field_E=E_z;
            end

            if number_of_layers==1 || l==1
                if conical_field_component==1 || conical_field_component==2 || conical_field_component==3 
                    field_H_2=field_H;
                elseif conical_field_component==4 || conical_field_component==5 || conical_field_component==6
                    field_E_2=field_E;
                end   
            elseif number_of_layers>=2
                if conical_field_component==1 || conical_field_component==2 || conical_field_component==3 
                    field_H_2=[field_H_2 field_H];
                elseif conical_field_component==4 || conical_field_component==5 || conical_field_component==6
                    field_E_2=[field_E_2 field_E];
                end    
            end
            if conical_field_component==4
                epsilon_matrix_local=[epsilon_matrix_local epsilon_matrix_local_];
                clear epsilon_matrix_local_
            end
        elseif polarization==1
            for z_l=1:1:length(z)
                %epsilon_matrix_local_=zeros(length(x),length(z));
%D Vypočítané konstanty dosadíme do soustavy vázaných rovnic.                
%D \placeformula[-]
%D \startformula
%D U_{l,y,i}=\sum^n_{m=1}w_{l,i,m}\left\lbrace c^+_{l,m}\e^{-k_0\lambda_{l,m}(z-D_l+d_l)}
%D +c^-_{l,m}\e^{k_0\lambda_{l,m}(z-D_l)} \right\rbrace,
%D \stopformula

%D \placeformula[-]
%D \startformula
%D S_{l,x,i}=\sum^n_{m=1}v_{l,i,m}\left\lbrace -c^+_{l,m}\e^{-k_0\lambda_{l,m}(z-D_l+d_l)}
%D +c^-_{l,m}\e^{k_0\lambda_{l,m}(z-D_l)} \right\rbrace,
%D \stopformula                
                c_plus_exp=c_plus.*exp(-k_0*diag(Q(:,:,l))*z_plus(z_l));
                c_minus_exp=c_minus.*exp(k_0*diag(Q(:,:,l))*z_minus(z_l));
%D Vypočítáme si harmoniky

                if planar_field_component==1
                    %if ARCWA==1
                    %    S_harmonics=V(:,:,l)*(-c_plus_exp+c_minus_exp); % --- dostanu zákmity, ARCWA --- problemy se zobrazenim
                    %else
                        G_harmonics=1i*sqrt(epsilon*mu)*W(:,:,l)*Q(:,:,l)*(-c_plus_exp+c_minus_exp);
                    %end
                elseif planar_field_component==2
                    U_harmonics=W(:,:,l)*(c_plus_exp+c_minus_exp);
                elseif planar_field_component==3
                    U_harmonics=W(:,:,l)*(c_plus_exp+c_minus_exp);
                    f_harmonics=-sqrt(mu/epsilon)*(E_l(:,:,l)\K_x)*U_harmonics;
                end

                if planar_field_component==1
                    %E_x(:,z_l)=1i*sqrt(mu/epsilon)*S_harmonics'*exp_xx; % --- %dostanu zákmity
                    D_x(:,z_l)=G_harmonics'*exp_xx;
                    for x_l=1:1:length(x)
                        field2_correct; % správně vypočítané E_x
                        field2_correct_local_ab;
                    end
                elseif planar_field_component==2
                    H_y(:,z_l)=U_harmonics'*exp_xx;
                elseif planar_field_component==3
                    E_z(:,z_l)=f_harmonics'*exp_xx;
                end


            end
            if planar_field_component==1
                field_E=E_x; % zapamatujeme si field v každé vrstvičce
                if number_of_layers==1 || l==1
                    field_E_2=field_E;
                elseif number_of_layers>=2
                    field_E_2=[field_E_2 field_E];
                end
                epsilon_matrix_local=[epsilon_matrix_local epsilon_matrix_local_];
                clear epsilon_matrix_local_
            elseif planar_field_component==2
                field_H=H_y; % zapamatujeme si field v každé vrstvičce
                if number_of_layers==1 || l==1
                    field_H_2=field_H;
                elseif number_of_layers>=2
                    field_H_2=[field_H_2 field_H];
                end
            elseif planar_field_component==3
                field_E=E_z; % zapamatujeme si field v každé vrstvičce
                if number_of_layers==1 || l==1
                    field_E_2=field_E;
                elseif number_of_layers>=2
                    field_E_2=[field_E_2 field_E];
                end
            end
            %field_D_x(:,:,l)=D_x;
        elseif polarization==2 %TE polarization --- stejný postup odvození
            for z_l=1:1:length(z)
                
                c_plus_exp=c_plus.*exp(-k_0*diag(Q(:,:,l))*z_plus(z_l));
                c_minus_exp=c_minus.*exp(k_0*diag(Q(:,:,l))*z_minus(z_l));
                
                if planar_field_component==1
                    U_harmonics=V(:,:,l)*(-c_plus_exp+c_minus_exp);
                elseif planar_field_component==2
                    S_harmonics=W(:,:,l)*(c_plus_exp+c_minus_exp);
                elseif planar_field_component==3
                    U_z_harmonics=1i*K_x*W(:,:,l)*(c_plus_exp+c_minus_exp);
                end
                
                if planar_field_component==1                            
                    H_x(:,z_l)=1i*sqrt(epsilon/mu)*U_harmonics'*exp_xx;                        
                elseif planar_field_component==2                            
                    E_y(:,z_l)=S_harmonics'*exp_xx; 
                    for x_l=1:1:length(x)   
                        field2_correct_local_ab;
                    end
                elseif planar_field_component==3                            
                    H_z(:,z_l)=-1i*sqrt(epsilon/mu)*U_z_harmonics'*exp_xx;                        
                end
                
            end
            if planar_field_component==1
                field_H=H_x; % zapamatujeme si field v každé vrstvičce
                if number_of_layers==1 || l==1
                    field_H_2=field_H;
                elseif number_of_layers>=2
                    field_H_2=[field_H_2 field_H];
                end
            elseif planar_field_component==2
                field_E=E_y; % zapamatujeme si field v každé vrstvičce
                if number_of_layers==1 || l==1
                    field_E_2=field_E;
                elseif number_of_layers>=2
                    field_E_2=[field_E_2 field_E];
                end
                epsilon_matrix_local=[epsilon_matrix_local epsilon_matrix_local_];
            elseif planar_field_component==3
                field_H=H_z; % zapamatujeme si field v každé vrstvičce
                if number_of_layers==1 || l==1
                    field_H_2=field_H;
                elseif number_of_layers>=2
                    field_H_2=[field_H_2 field_H];
                end
            end
        end
        if matrix_algorithm==2
            T_l=a_1_X(:,:,l)*T_l; % postoupíme do další vrstvičky, potřebujeme nové $T_l$
        end
        D_l_pred=D_l; % výpočet tloušťky nad budoucí vrstvičkou
    end
    
    %D Pokud je navolena pouze jedna vrstvička je vše již jasné, pokud je
    %D navoleno více vrstviček, vytváříme matici postupným přidáváním
    %D polí z jednotlivých vrstviček.
   
end