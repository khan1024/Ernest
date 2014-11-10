% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% M. G. Moharam, Drew A. Pommet, Eric B. Grann, and T. K. Gaylord, 
% "Stable implementation of the rigorous coupled-wave analysis for surface-relief
% gratings: enhanced transmittance matrix approach," J. Opt. Soc. Am. A 12, 1077-1086 (1995) 

% M. G. Moharam, Eric B. Grann, Drew A. Pommet, and T. K. Gaylord, 
% "Formulation for stable and efficient implementation of the rigorous coupled-wave
% analysis of binary gratings," J. Opt. Soc. Am. A 12, 1068-1076 (1995)

for l=number_of_layers:-1:1
    if polarization==0 
        X1(:,:,l)=diag(exp(-k_0*diag(Q1(:,:,l))*layer_thickness(l)));    
        X2(:,:,l)=diag(exp(-k_0*diag(Q2(:,:,l))*layer_thickness(l)));    
        X(:,:,l)=[X1(:,:,l) zero; zero X2(:,:,l)];    
        
        if change_matrix_base==1
            auxiliary_matrix=[Vss(:,:,l) Vsp(:,:,l) Vss(:,:,l) Vsp(:,:,l) ; ...
                            Wss(:,:,l) Wsp(:,:,l) -Wss(:,:,l) -Wsp(:,:,l) ; ...
                            Wps(:,:,l) Wpp(:,:,l) -Wps(:,:,l) -Wpp(:,:,l) ; ...
                            Vps(:,:,l) Vpp(:,:,l) Vps(:,:,l) Vpp(:,:,l)];
        elseif change_matrix_base==2
            new_W(:,:,l)=[zero W2(:,:,l); V11(:,:,l) V12(:,:,l)];
            new_V(:,:,l)=-[W1(:,:,l) zero; V21(:,:,l) V22(:,:,l)];
        end
    else
        auxiliary_matrix=[W(:,:,l) W(:,:,l); V(:,:,l) -V(:,:,l)];   
        X(:,:,l)=diag(exp(-k_0*diag(Q(:,:,l))*layer_thickness(l)));
    end
    %--------------------------------------------------------------------
    %D \subsubject{Okrajové podmínky}

    %D Nejprve si vypočítám
    %D \placeformula[]
    %D \startformula
    %D \pmatrix{\mat{a}_L \cr \mat{b}_L}=\pmatrix{\mat{W}_L & \mat{W}_L\cr
    %D \mat{V}_L & -\mat{V}_L}^{-1}\pmatrix{\mat{f}_{L+1} \cr \mat{g}_{L+1}}
    %D \stopformula

    %D Pro výpočet inverzní matice se doporučuje dle článku použít single-value
    %D decomposition, mělo by se tím zabránit ztrátě přesnosti, spíše to zatím
    %D je jenom na optíž.

    %D Při některých výpočtech se objeví, že matice \type{a_b}
    %D nebo \type{a_1} jsou singulární, výpočet inverzní matice
    %D je potom nepřesný. Můžeme tomu zabránit pokud použijeme
    %D výpočet inverzní matice pomocí SVD. !!!!Pozor!!!! --- pro velké 
    %D matice je SVD velmi pomalé
    
    if change_matrix_base==2 && polarization==0
        
        inv_new_W_f_new=new_W(:,:,l)\f_new;
        inv_new_V_g_new=new_V(:,:,l)\g_new;
        
        a=0.5*(inv_new_W_f_new+inv_new_V_g_new);
        b=0.5*(inv_new_W_f_new-inv_new_V_g_new);
        
        %a=0.5*(inv(new_W(:,:,l))*f_new+inv(new_V(:,:,l))*g_new);
        %b=0.5*(inv(new_W(:,:,l))*f_new-inv(new_V(:,:,l))*g_new);
    else
        %a_b=inv(auxiliary_matrix)*f_g;
        a_b=auxiliary_matrix\f_g;
        %clear pomoc_U pomoc_S pomoc_V
        %D Vyjádřím si \type{a} a \type{b}
        if polarization==0
            a=a_b(1:2*number_of_orders, : );
            b=a_b(2*number_of_orders+1:size(a_b), : );    
        else
            a=a_b(1:number_of_orders, : );        
            b=a_b(number_of_orders+1:size(a_b), : );    
        end
    end

    %D Vyjádřím si důležitý člen, který je přítomný v některých maticích
    %D $\mat{X}_L\mat{b}_L\mat{a}^{-1}_L\mat{X}_L$

    %D Pro špatně podmíněné matice.
    %if rcond(a)<0.00005    
    %    [a_U,a_S,a_V]=svd(a);    
    %    a_1=a_V*inv(a_S)*a_U';
    %else
    %    a_1=inv(a);
    %end

    %clear a_U a_S a_V

    %a_1_X(:,:,l)=a_1*X(:,:,l);
    
    a_1_X(:,:,l)=a\X(:,:,l);
    b_a_1_X(:,:,l)=b*a_1_X(:,:,l);
    X_b_a_1_X=X(:,:,l)*b*a_1_X(:,:,l);

    %D Při postupu mřížkou nahoru musíme vypočítat nové \type{f_g} podle vzorce
    %D \placeformula[]
    %D \startformula
    %D \pmatrix{\mat{f}_{L} \cr \mat{g}_{L}}=\pmatrix{\mat{W}_L\left(\mat{I}+
    %D \mat{X}_L\mat{b}_L\mat{a}^{-1}_L\mat{X}_L \right)
    %D \cr \mat{V}_L\left(\mat{I}-\mat{X}_L\mat{b}_L\mat{a}^{-1}_L\mat{X}_L \right)}
    %D \stopformula

    if polarization==0
        if change_matrix_base==1
            f_g=auxiliary_matrix*[I zero ; zero I; X_b_a_1_X];
            clear auxiliary_matrix
        elseif change_matrix_base==2
            new_I=[I zero; zero I];
            f_g=[new_W(:,:,l)*(new_I+X_b_a_1_X) ; new_V(:,:,l)*(new_I-X_b_a_1_X)];
            f_new=f_g(1:2*number_of_orders,1:2*number_of_orders);
            g_new=f_g(2*number_of_orders+1:4*number_of_orders,1:2*number_of_orders);
        end
    else % TE, TM polarization    
        f_g=[W(:,:,l)*(I+X_b_a_1_X) ; V(:,:,l)*(I-X_b_a_1_X)];
    end
    
    %D Při průchodu budu vypočítávat \type{new_T}
    new_T=new_T*a_1_X(:,:,l);
end

%D \subsubject{Výpočet reflexí a transmisí}

%D Pomocí vzorce
%D \placeformula[-]
%D \startformula
%D \pmatrix{\mat{R} \cr \mat{T}_1}=\pmatrix{-\mat{I} & \mat{W}_1
%D \left(\mat{I}+\mat{X}_1\mat{b}_1\mat{a}^{-1}_1\mat{X}_1 \right) 
%D \cr \i\mat{Z}_1 & \mat{V}_1\left(\mat{I}-\mat{X}_1\mat{b}_1\mat{a}^{-1}_1
%D \mat{X}_1 \right)}^{-1}\pmatrix{\delta_{i0} \cr \i\delta_{i0}\cos(\theta)/n_1}.
%D \stopformula
%D vypočítám \type{R} a \type{T1}

if polarization==1 %TM polarization
    transfer=[-I W(:,:,1)*(I+X_b_a_1_X);1i*Z1 V(:,:,1)*(I-X_b_a_1_X)];
    R_T1=transfer\incident.';
    %R_T1=inv(transfer)*incident.';
    
    R=R_T1(1:number_of_orders);
    T1=R_T1(number_of_orders+1:size(R_T1));
    T=new_T*T1;
elseif polarization==2 %TE polarization
    transfer=[-I W(:,:,1)*(I+X_b_a_1_X);1i*Y1 V(:,:,1)*(I-X_b_a_1_X)];
    R_T1=transfer\incident.';
    %R_T1=inv(transfer)*incident.';
    
    R=R_T1(1:number_of_orders);
    T1=R_T1(number_of_orders+1:size(R_T1));
    T=new_T*T1;
else
    if change_matrix_base==1
        transfer=[-I zero; 1i*Y1 zero; zero -I; zero 1i*Z1];    
        transfer( : , 2*number_of_orders+1:4*number_of_orders)=f_g;    
        R_T1=transfer\incident.';    
        R_s=R_T1(1:number_of_orders);    
        R_p=R_T1(number_of_orders+1:2*number_of_orders);    
        T1=R_T1(2*number_of_orders+1:size(R_T1));    
    
        T=new_T*T1;    
        T_s=T(1:number_of_orders);    
        T_p=T(number_of_orders+1:size(T));
    elseif change_matrix_base==2
        velky_M=[new_I; 1i*K_y_n*inv_Y1*K_x, 1i*(Y1+K_y_n*inv_Y1*K_y_n); -1i*(Y1+K_x*inv_Y1*K_x), -1i*K_x*inv_Y1*K_y_n];
    
        transfer=[-velky_M, f_g];    
        R_T1=transfer\incident.';
    
        R_x=R_T1(1:number_of_orders);    
        R_y=R_T1(number_of_orders+1:2*number_of_orders);    
        T1=R_T1(2*number_of_orders+1:size(R_T1));    
    
        T=new_T*T1;    
        T_x=T(1:number_of_orders);    
        T_y=T(number_of_orders+1:size(T));
        R_z=(R_x.*k_x + R_y.*k_y)./k_1_z;
        T_z=-(T_x.*k_x + T_y.*k_y)./k_3_z;
        clear velky_M
    end
end
clear transfer 