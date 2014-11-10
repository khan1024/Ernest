
z_points=z_points_3;

exp_alpha_x_2=zeros(length(x_points),length(n_K));

for x_p=1:1:length(x_points)         
    exp_alpha_x_2(x_p,:)=exp(-1i*alpha_m*1E6*x_points(x_p));
end

if plot_field_component>=3
    T_h=-W_2(:,:,1)*T;
    T_x_h=T_h(1:n_sq);
    T_y_h=T_h(n_sq+1:2*n_sq);
end

profile_mode3=[];
for z_points_loop=1:1:length(z_points)
    xz_slice_coordinate_z=z_points(z_points_loop);
    z_minus=xz_slice_coordinate_z-sum(thickness);

    y_points=1;
    if symmetry==1    
        % u_plus a d_minus jsou spocitany na horni a dolni casti vrstvicky
        %X_plus_u_plus=c_p(:,:,xz_slice_section_to_view_reduce).*exp(z_plus*1i*diag(gamma_new(:,:,matice_stejnosti_final(xz_slice_section_to_view,3))));       
        X_plus_u_plus=c_p(:,:,xz_slice_section_to_view_reduce).*exp(z_plus*1i*diag(eval(['gamma_new_' num2str(matice_stejnosti_final(xz_slice_section_to_view,3))])));       
        % pred eig je minus, ale <0 je take y_minus
        %X_minus_d_minus=c_m(:,:,xz_slice_section_to_view_reduce).*exp(-z_minus*1i*diag(gamma_new(:,:,matice_stejnosti_final(xz_slice_section_to_view,3))));
        X_minus_d_minus=c_m(:,:,xz_slice_section_to_view_reduce).*exp(-z_minus*1i*diag(eval(['gamma_new_' num2str(matice_stejnosti_final(xz_slice_section_to_view,3))])));
    elseif symmetry==2    
        
    end    
    
    if plot_field_component==1 % E_1    
        E_1=zeros(length(x_points),length(y_points));
        E_1_harmonics3_v=T_x_Li2.*exp((1)*-z_minus*1i*diag(full(gamma_3)));
        E_1_harmonics3=reshape(E_1_harmonics3_v',n_x,n_x)';
    elseif plot_field_component==2 % E_2    
        E_2=zeros(length(x_points),length(y_points));    
        E_2_harmonics3_v=T_y_Li2.*exp((1)*-z_minus*1i*diag(full(gamma_3)));    
        E_2_harmonics3=reshape(E_2_harmonics3_v',n_x,n_x)';
    elseif plot_field_component==3 % E_3    
        E_3=zeros(length(x_points),length(y_points));    
        T_z_Li2=(alpha*T_x_Li2+beta*T_y_Li2)./diag(gamma_3);    
        E_3_harmonics3_v=T_z_Li2.*exp((1)*-z_minus*1i*diag(full(gamma_3)));    
        E_3_harmonics3=reshape(E_3_harmonics3_v',n_x,n_x)';
    elseif plot_field_component==4 % H_1    
        H_1=zeros(length(x_points),length(y_points));
        H_1_harmonics3_v=T_x_h.*exp((1)*-z_minus*1i*diag(full(gamma_3)));    
        H_1_harmonics3=reshape(H_1_harmonics3_v',n_x,n_x)';
    elseif plot_field_component==5 % H_2    
        H_2=zeros(length(x_points),length(y_points));
        H_2_harmonics3_v=(-1)*T_y_h.*exp((1)*-z_minus*1i*diag(full(gamma_3)));    
        H_2_harmonics3=reshape(H_2_harmonics3_v',n_x,n_x)';
    elseif plot_field_component==6 % H_3    
        H_3=zeros(length(x_points),length(y_points));
        T_z_h=-(alpha*T_x_h+beta*T_y_h)./diag(gamma_3);    
        H_3_harmonics3_v=T_z_h.*exp((1)*-z_minus*1i*diag(full(gamma_3)));    
        H_3_harmonics3=reshape(H_3_harmonics3_v',n_x,n_x)';   
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
        %components_no_symmetry;
    end
    
    y_p=1;
    
    exp_alpha_y=exp(-1i*beta_n*xz_slice_coordinate_y);
    
    % predpocitam si nasledujici vektory
    
    %{
    if plot_field_component==1 % E_1            
        harmonics2_vector=E_1_harmonics3*exp_alpha_y';
    elseif plot_field_component==2 % E_2            
        harmonics2_vector=E_2_harmonics3*exp_alpha_y';
    elseif plot_field_component==3 % E_3            
        harmonics2_vector=E_3_harmonics3*exp_alpha_y';
    elseif plot_field_component==4 % H_1            
        harmonics2_vector=H_1_harmonics3*exp_alpha_y';
    elseif plot_field_component==5 % H_2            
        harmonics2_vector=H_2_harmonics3*exp_alpha_y';
    elseif plot_field_component==6 % H_3            
        harmonics2_vector=H_3_harmonics3*exp_alpha_y';    
    end
    %}
    
    %{
    for x_p=1:1:length(x_points)
        
        %exp_alpha_x=exp(-1i*2*pi*m_K/Lambda_x*1E6*x_points(x_p));
        
        if plot_field_component==1 % E_1
            E_1(x_p,y_p)=sum(harmonics2_vector.*exp_alpha_x_2(x_p,:)');
        elseif plot_field_component==2 % E_2
            E_2(x_p,y_p)=sum(harmonics2_vector.*exp_alpha_x_2(x_p,:)');
        elseif plot_field_component==3 % E_3
            E_3(x_p,y_p)=sum(harmonics2_vector.*exp_alpha_x_2(x_p,:)');
        elseif plot_field_component==4 % H_1
            H_1(x_p,y_p)=sum(harmonics2_vector.*exp_alpha_x_2(x_p,:)');
        elseif plot_field_component==5 % H_2
            H_2(x_p,y_p)=sum(harmonics2_vector.*exp_alpha_x_2(x_p,:)');
        elseif plot_field_component==6 % H_3
            H_3(x_p,y_p)=sum(harmonics2_vector.*exp_alpha_x_2(x_p,:)');
        end
    end
    
    %}
    %%{
    if plot_field_component==1 % E_1
        E_1(:,y_p)=sum(bsxfun(@times,(E_1_harmonics3*exp_alpha_y'),exp_alpha_x_2'));      
    elseif plot_field_component==2 % E_2
        E_2(:,y_p)=sum(bsxfun(@times,(E_2_harmonics3*exp_alpha_y'),exp_alpha_x_2'));       
    elseif plot_field_component==3 % E_3
        E_3(:,y_p)=sum(bsxfun(@times,(E_3_harmonics3*exp_alpha_y'),exp_alpha_x_2'));      
    elseif plot_field_component==4 % H_1
        H_1(:,y_p)=sum(bsxfun(@times,(H_1_harmonics3*exp_alpha_y'),exp_alpha_x_2'));       
    elseif plot_field_component==5 % H_2
        H_2(:,y_p)=sum(bsxfun(@times,(H_2_harmonics3*exp_alpha_y'),exp_alpha_x_2'));       
    elseif plot_field_component==6 % H_3
        H_3(:,y_p)=sum(bsxfun(@times,(H_3_harmonics3*exp_alpha_y'),exp_alpha_x_2'));       
    end
    %%}
    
    if plot_field_component==1
        profile_mode=zobrazeni_pole(E_1,field_function);
    elseif plot_field_component==2
        profile_mode=zobrazeni_pole(E_2,field_function);
    elseif plot_field_component==3
        profile_mode=zobrazeni_pole(E_3,field_function);
    elseif plot_field_component==4
        profile_mode=zobrazeni_pole(-H_1,field_function); 
    elseif plot_field_component==5
        profile_mode=zobrazeni_pole(H_2,field_function); 
    elseif plot_field_component==6
        profile_mode=zobrazeni_pole(H_3,field_function); 
    end
    profile_mode3=[profile_mode3, profile_mode];
end

%{
figure
set(gcf, 'Renderer', 'zbuffer');
axes('FontSize',13)
surf(z_points,1E6.*x_points,profile_mode3,'LineWidth',1);

shading flat
axis([z_points(1) z_points(end) coordinate_x_1(1) coordinate_x_1(end)]);
ylabel('x [\mu{}m]','FontSize',13);
xlabel('z [\mu{}m]','FontSize',13);
%view([90 90]);

if set_waveguide_framework==1
hold on;
waveguide_framework_xz;
end

view([0 -90])
%}

