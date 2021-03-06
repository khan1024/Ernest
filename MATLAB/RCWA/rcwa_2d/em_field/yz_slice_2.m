
z_points=z_points_2;
exp_alpha_y_2=zeros(length(y_points),length(n_K));

for y_p=1:1:length(y_points)         
    exp_alpha_y_2(y_p,:)=exp(-1i*beta_n*1E6*y_points(y_p));
end

profile_mode2=[];
for z_points_loop=1:1:length(z_points)
    yz_slice_coordinate_z=z_points(z_points_loop);    
    
    total_thickness(1)=0;
    for i=1:1:length(thickness)    
        total_thickness(i+1)=total_thickness(i)+thickness(i);    
        if yz_slice_coordinate_z>=total_thickness(i) && yz_slice_coordinate_z<total_thickness(i+1)        
            yz_slice_section_to_view=i;       
            break    
        end        
    end
    
    if z_points_loop==1
        first_yz_slice_section_to_view=yz_slice_section_to_view;
    end
    
    u_plus=c_p(:,:,yz_slice_section_to_view);
    d_minus=c_m(:,:,yz_slice_section_to_view);

    D_l_pred=total_thickness(yz_slice_section_to_view);
    D_l=D_l_pred+thickness(yz_slice_section_to_view);
    z_plus=yz_slice_coordinate_z-D_l_pred;
    z_minus=yz_slice_coordinate_z-D_l;


    x_points=1;
    if symmetry==1    
        % u_plus a d_minus jsou spocitany na horni a dolni casti vrstvicky
        %X_plus_u_plus=c_p(:,:,yz_slice_section_to_view_reduce).*exp(z_plus*1i*diag(gamma_new(:,:,matice_stejnosti_final(yz_slice_section_to_view,3))));        
        X_plus_u_plus=c_p(:,:,yz_slice_section_to_view_reduce).*exp(z_plus*1i*diag(eval(['gamma_new_' num2str(matice_stejnosti_final(yz_slice_section_to_view,3))])));        
        % pred eig je minus, ale <0 je take y_minus
        %X_minus_d_minus=c_m(:,:,yz_slice_section_to_view_reduce).*exp(-z_minus*1i*diag(gamma_new(:,:,matice_stejnosti_final(yz_slice_section_to_view,3))));
        X_minus_d_minus=c_m(:,:,yz_slice_section_to_view_reduce).*exp(-z_minus*1i*diag(eval(['gamma_new_' num2str(matice_stejnosti_final(yz_slice_section_to_view,3))])));
    elseif symmetry==2
        % u_plus a d_minus jsou spocitany na horni a dolni casti vrstvicky
        %X_plus_u_plus=c_p(:,:,yz_slice_section_to_view_reduce).*exp(z_plus*1i*diag(gamma(:,:,matice_stejnosti_final(yz_slice_section_to_view,3))));        
        X_plus_u_plus=c_p(:,:,yz_slice_section_to_view).*exp(-z_plus*1i*diag(eval(['gamma_2_' num2str(matice_stejnosti_final(yz_slice_section_to_view,3))])));        
        % pred eig je minus, ale <0 je take y_minus
        %X_minus_d_minus=c_m(:,:,yz_slice_section_to_view_reduce).*exp(-z_minus*1i*diag(gamma(:,:,matice_stejnosti_final(yz_slice_section_to_view,3))));
        X_minus_d_minus=c_m(:,:,yz_slice_section_to_view).*exp(-(-1)*z_minus*1i*diag(eval(['gamma_2_' num2str(matice_stejnosti_final(yz_slice_section_to_view,3))])));
    end    
    
    if plot_field_component==1 % E_1    
        E_1=zeros(length(x_points),length(y_points));
        
        if z_points_loop==1
            E_1_plus_eigvec=eigvec_E(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
        else
            if first_yz_slice_section_to_view==yz_slice_section_to_view
            else
                first_yz_slice_section_to_view=yz_slice_section_to_view;
                E_1_plus_eigvec=eigvec_E(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
            end
        end
        
        %E_1_plus_eigvec=eigvec_E(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));    
        %E_1_minus_eigvec=E_1_plus_eigvec;    
        %E_1_harmonics=E_1_plus_eigvec*X_plus_u_plus+E_1_minus_eigvec*X_minus_d_minus;
        
        E_1_harmonics=E_1_plus_eigvec*(X_plus_u_plus+X_minus_d_minus);
    elseif plot_field_component==2 % E_2    
        E_2=zeros(length(x_points),length(y_points));
        
        if z_points_loop==1
            E_2_plus_eigvec=eigvec_E(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
        else
            if first_yz_slice_section_to_view==yz_slice_section_to_view
            else
                first_yz_slice_section_to_view=yz_slice_section_to_view;
                E_2_plus_eigvec=eigvec_E(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
            end
        end
    
        %E_2_plus_eigvec=eigvec_E(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));    
        %E_2_minus_eigvec=E_2_plus_eigvec;    
        %E_2_harmonics=E_2_plus_eigvec*X_plus_u_plus+E_2_minus_eigvec*X_minus_d_minus;
        
        
        E_2_harmonics=E_2_plus_eigvec*(X_plus_u_plus+X_minus_d_minus);

    elseif plot_field_component==3 % E_3    
        E_3=zeros(length(x_points),length(y_points));  
        
        if z_points_loop==1
            H_1_plus_eigvec=eigvec_H(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
        else
            if first_yz_slice_section_to_view==yz_slice_section_to_view
            else
                first_yz_slice_section_to_view=yz_slice_section_to_view;
                H_1_plus_eigvec=eigvec_H(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
            end
        end
    
        %H_1_plus_eigvec=eigvec_H(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));    
        %H_1_minus_eigvec=-H_1_plus_eigvec;    
        %H_1_harmonics=H_1_plus_eigvec*X_plus_u_plus+H_1_minus_eigvec*X_minus_d_minus;
    
        %H_2_plus_eigvec=eigvec_H(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));    
        %H_2_minus_eigvec=-H_2_plus_eigvec;    
        %H_2_harmonics=H_2_plus_eigvec*X_plus_u_plus+H_2_minus_eigvec*X_minus_d_minus;
        
        H_1_harmonics=H_1_plus_eigvec*(X_plus_u_plus-X_minus_d_minus);
        H_2_harmonics=H_1_harmonics; %H_2_plus_eigvec*(X_plus_u_plus-X_minus_d_minus);
    elseif plot_field_component==4 % H_1    
        H_1=zeros(length(x_points),length(y_points));
        
        if z_points_loop==1
            H_1_plus_eigvec=eigvec_H(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
        else
            if first_yz_slice_section_to_view==yz_slice_section_to_view
            else
                first_yz_slice_section_to_view=yz_slice_section_to_view;
                H_1_plus_eigvec=eigvec_H(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
            end
        end
    
        %H_1_plus_eigvec=eigvec_H(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));    
        %H_1_minus_eigvec=-H_1_plus_eigvec;    
        %H_1_harmonics=H_1_plus_eigvec*X_plus_u_plus+H_1_minus_eigvec*X_minus_d_minus;
        
        H_1_harmonics=H_1_plus_eigvec*(X_plus_u_plus-X_minus_d_minus);
    elseif plot_field_component==5 % H_2    
        H_2=zeros(length(x_points),length(y_points));
        
        if z_points_loop==1
            H_2_plus_eigvec=eigvec_H(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
        else
            if first_yz_slice_section_to_view==yz_slice_section_to_view
            else
                first_yz_slice_section_to_view=yz_slice_section_to_view;
                H_2_plus_eigvec=eigvec_H(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
            end
        end
    
        %H_2_plus_eigvec=eigvec_H(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));    
        %H_2_minus_eigvec=-H_2_plus_eigvec;    
        %H_2_harmonics=H_2_plus_eigvec*X_plus_u_plus+H_2_minus_eigvec*X_minus_d_minus;
        
        H_2_harmonics=H_2_plus_eigvec*(X_plus_u_plus-X_minus_d_minus);
    elseif plot_field_component==6 % H_3    
        H_3=zeros(length(x_points),length(y_points));
        
        
       if z_points_loop==1
            E_1_plus_eigvec=eigvec_E(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
        else
            if first_yz_slice_section_to_view==yz_slice_section_to_view
            else
                first_yz_slice_section_to_view=yz_slice_section_to_view;
                E_1_plus_eigvec=eigvec_E(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));
            end
        end
        
        %E_1_plus_eigvec=eigvec_E(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));    
        %E_1_minus_eigvec=E_1_plus_eigvec;    
        %E_1_harmonics=E_1_plus_eigvec*X_plus_u_plus+E_1_minus_eigvec*X_minus_d_minus;    
    
        %E_2_plus_eigvec=eigvec_E(:,:,matice_stejnosti_final(yz_slice_section_to_view,3));    
        %E_2_minus_eigvec=E_2_plus_eigvec;    
        %E_2_harmonics=E_2_plus_eigvec*X_plus_u_plus+E_2_minus_eigvec*X_minus_d_minus;
        
        E_1_harmonics=E_1_plus_eigvec*(X_plus_u_plus+X_minus_d_minus);
        E_2_harmonics=E_1_harmonics; %E_2_plus_eigvec*(X_plus_u_plus+X_minus_d_minus);
    end    
    
    % indexy (i,j) pro N_X=2
    % (-2,-2), (-2,-1), (-2,0), (-2,1), (-2,2)
    % (-1,-2), (-1,-1), (-1,0), (-1,1), (-1,2)
    %  (0,-2),  (0,-1),  (0,0),  (0,1),  (0,2)
    %  (1,-2),  (1,-1),  (1,0),  (1,1),  (1,2)
    %  (2,-2),  (2,-1),  (2,0),  (2,1),  (2,2)

    if symmetry==1    
        if symmetry_switch==1 || symmetry_switch==2        
            % x-symmetry        
            components_x_symmetry;    
        elseif symmetry_switch==3 || symmetry_switch==4        
            % y_symmetry        
            components_y_symmetry;    
        elseif symmetry_switch==5 || symmetry_switch==6        
            % x_y symmetry        
            components_x_y_symmetry_1;    
        elseif symmetry_switch==7 || symmetry_switch==8        
            % x_y symmetry        
            error('not yet implemented')        
            %components_x_y_symmetry_2;    
        end        
    elseif symmetry==2    
        components_no_symmetry;
    end    
    
    x_p=1;
    
    exp_alpha_x=exp(-1i*alpha_m*yz_slice_coordinate_x);
    
    %{
    for y_p=1:1:length(y_points)
        
        %exp_alpha_y=exp(-1i*2*pi*n_K/Lambda_y*1E6*y_points(y_p));
        
        if plot_field_component==1 % E_1
            E_1(x_p,y_p)=sum((E_1_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x');
        elseif plot_field_component==2 % E_2
            E_2(x_p,y_p)=sum((E_2_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x');
        elseif plot_field_component==3 % E_3
            E_3(x_p,y_p)=sum((E_3_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x');
        elseif plot_field_component==4 % H_1
            H_1(x_p,y_p)=sum((H_1_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x');
        elseif plot_field_component==5 % H_2
            H_2(x_p,y_p)=sum((H_2_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x');
        elseif plot_field_component==6 % H_3
            H_3(x_p,y_p)=sum((H_3_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x');
        end        
    end
    %}
    
    %%{
    if plot_field_component==1 % E_1
        E_1(x_p,:)=sum(bsxfun(@times,(E_1_harmonics2*exp_alpha_y_2'),exp_alpha_x'));       
    elseif plot_field_component==2 % E_2
        E_2(x_p,:)=sum(bsxfun(@times,(E_2_harmonics2*exp_alpha_y_2'),exp_alpha_x'));      
    elseif plot_field_component==3 % E_3
        E_3(x_p,:)=sum(bsxfun(@times,(E_3_harmonics2*exp_alpha_y_2'),exp_alpha_x'));     
    elseif plot_field_component==4 % H_1
        H_1(x_p,:)=sum(bsxfun(@times,(H_1_harmonics2*exp_alpha_y_2'),exp_alpha_x'));       
    elseif plot_field_component==5 % H_2
        H_2(x_p,:)=sum(bsxfun(@times,(H_2_harmonics2*exp_alpha_y_2'),exp_alpha_x'));        
    elseif plot_field_component==6 % H_3
        H_3(x_p,:)=sum(bsxfun(@times,(H_3_harmonics2*exp_alpha_y_2'),exp_alpha_x'));       
    end
    %%}
    
    %H_3(x_p,:)=sum(bsxfun(@times,(H_3_harmonics2*exp_alpha_y_2'),exp_alpha_x')); 
    
    
    if plot_field_component==1
        profile_mode=zobrazeni_pole(E_1,field_function);
    elseif plot_field_component==2
        profile_mode=zobrazeni_pole(E_2,field_function);
    elseif plot_field_component==3
        profile_mode=zobrazeni_pole(E_3,field_function);
    elseif plot_field_component==4
        profile_mode=zobrazeni_pole(H_1,field_function);
    elseif plot_field_component==5    
        profile_mode=zobrazeni_pole(H_2,field_function); 
    elseif plot_field_component==6    
        profile_mode=zobrazeni_pole(H_3,field_function); 
    end

    profile_mode2=[profile_mode2; profile_mode];
end
%{
figure
set(gcf, 'Renderer', 'zbuffer');
axes('FontSize',13)
surf(1E6.*y_points,z_points,profile_mode2,'LineWidth',1);

shading flat
axis([coordinate_y_1(1) coordinate_y_1(end) z_points(1) z_points(end)]);
ylabel('z [\mu{}m]','FontSize',13);
xlabel('y [\mu{}m]','FontSize',13);
%view([90 90]);

if set_waveguide_framework==1
hold on;
waveguide_framework_yz;
end

%view([0 -90])

view([90 90])
%}
