if plot_field_component==1 % E_1
    E_1_harmonics2=reshape(E_1_harmonics(1:n_sq)',n_x,n_x)'; % reshape vector -> matrix
elseif plot_field_component==2 % E_2   
    E_2_harmonics2=reshape(E_2_harmonics(n_sq+1:2*n_sq)',n_x,n_x)';
elseif plot_field_component==3 % E_3, S_3=1i*inv(eps_33)*(K_x*U_2-K_y*U_1)   
    

    if xy_slice_view==1
        to_view=xy_slice_section_to_view;
    elseif xz_slice_view==1
        to_view=xz_slice_section_to_view;
    elseif yz_slice_view==1
        to_view=yz_slice_section_to_view;
    end
    
    F_1_h=eval(['inv_E_' num2str(matice_stejnosti_final(to_view,3))])*(alpha./k_0);    
    F_2_h=eval(['inv_E_' num2str(matice_stejnosti_final(to_view,3))])*(beta./k_0);
    
    %F_1_h=inv_E*(alpha./k_0);    
    %F_2_h=inv_E*(beta./k_0);

    H_1_harmonics2=reshape((F_2_h*H_1_harmonics(1:n_sq))',n_x,n_x)';
    H_2_harmonics2=reshape((F_1_h*H_2_harmonics(n_sq+1:2*n_sq))',n_x,n_x)';
    
    E_3_harmonics2=(H_2_harmonics2-H_1_harmonics2);
elseif plot_field_component==4 % H_1  
    H_1_harmonics2=reshape(H_1_harmonics(1:n_sq)',n_x,n_x)';
elseif plot_field_component==5 % H_2   
    H_2_harmonics2=reshape(H_2_harmonics(n_sq+1:2*n_sq)',n_x,n_x)';
elseif plot_field_component==6 % H_3, U_3=1i*inv(mu_33)*(K_x*S_2-K_y*S_1)

    G_1_h=beta./k_0;    
    G_2_h=alpha./k_0;

    E_1_harmonics2=reshape((G_1_h*E_1_harmonics(1:n_sq))',n_x,n_x)';
    E_2_harmonics2=reshape((G_2_h*E_2_harmonics(n_sq+1:2*n_sq))',n_x,n_x)';
    
    H_3_harmonics2=-(E_2_harmonics2-E_1_harmonics2);
end