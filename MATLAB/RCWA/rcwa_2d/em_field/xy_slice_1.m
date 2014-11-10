
R_x_e=R(1:n_sq);
R_y_e=R(n_sq+1:2*n_sq);

H_inc=-W_2(:,:,end)*d_n1;
v_y=H_inc((n_x-1)/2+(n_x-1)*n_x/2+1+n_sq)*(-1); % ? (-1)
v_x=H_inc((n_x-1)/2+(n_x-1)*n_x/2+1);
v_z=(alpha_0*v_x-beta_0*v_y)/gamma_p00; % ? ...*v_x-beta_0*...

z_plus=xy_slice_coordinate_z;

if symmetry==1
    % u_plus a d_minus jsou spocitany na horni a dolni casti vrstvicky        
    %X_plus_u_plus=c_p(:,:,xy_slice_section_to_view_reduce).*exp(z_plus*1i*diag(gamma_new(:,:,matice_stejnosti_final(xy_slice_section_to_view,3))));        
    X_plus_u_plus=c_p(:,:,xy_slice_section_to_view_reduce).*exp(z_plus*1i*diag(eval(['gamma_new_' num2str(matice_stejnosti_final(xy_slice_section_to_view,3))])));        
    % pred eig je minus, ale <0 je take y_minus        
    %X_minus_d_minus=c_m(:,:,xy_slice_section_to_view_reduce).*exp(-z_minus*1i*diag(gamma_new(:,:,matice_stejnosti_final(xy_slice_section_to_view,3))));
    X_minus_d_minus=c_m(:,:,xy_slice_section_to_view_reduce).*exp(-z_minus*1i*diag(eval(['gamma_new_' num2str(matice_stejnosti_final(xy_slice_section_to_view,3))])));
elseif symmetry==2
    
end

x_points=x_point_b*1E-6:((x_point_e-x_point_b)*1E-6/resolution_x):(x_point_e*1E-6);
y_points=y_point_b*1E-6:((y_point_e-y_point_b)*1E-6/resolution_y):(y_point_e*1E-6);

if plot_field_component>=3
    R_h=W_2(:,:,end)*R;
    R_x_h=R_h(1:n_sq);
    R_y_h=R_h(n_sq+1:2*n_sq);
end

if plot_field_component==1 % E_1
    E_1=zeros(length(x_points),length(y_points));
    E_1_0=zeros(length(x_points),length(y_points));
    E_1_harmonics3_v=R_x_e.*exp(1i*z_plus*diag(full(gamma_1)));
    E_1_harmonics3=reshape(E_1_harmonics3_v',n_x,n_x)';
   
    E_1_harmonics_0=u_x/sqrt(normalize_2).*exp(-1i*z_plus*gamma_p00);
elseif plot_field_component==2 % E_2
    E_2=zeros(length(x_points),length(y_points));
    E_2_0=zeros(length(x_points),length(y_points));
    E_2_harmonics3_v=R_y_e.*exp(1i*z_plus*diag(full(gamma_1)));
    E_2_harmonics3=reshape(E_2_harmonics3_v',n_x,n_x)';
    
    E_2_harmonics_0=u_y/sqrt(normalize_2).*exp(-1i*z_plus*gamma_p00);
elseif plot_field_component==3 % E_3
    E_3=zeros(length(x_points),length(y_points));
    E_3_0=zeros(length(x_points),length(y_points));
    
    R_z_e=-(alpha*R_x_e+beta*R_y_e)./diag(gamma_1);
    
    E_3_harmonics3_v=R_z_e.*exp(1i*z_plus*diag(full(gamma_1)));
    E_3_harmonics3=reshape(E_3_harmonics3_v',n_x,n_x)';
    
    E_3_harmonics_0=-u_z/sqrt(normalize_2).*exp(-1i*z_plus*gamma_p00);
elseif plot_field_component==4 % H_1
    H_1=zeros(length(x_points),length(y_points));
    H_1_0=zeros(length(x_points),length(y_points));

    H_1_harmonics3_v=R_x_h.*exp(1i*z_plus*diag(full(gamma_1)));
    H_1_harmonics3=reshape(H_1_harmonics3_v',n_x,n_x)';
   
    H_1_harmonics_0=v_x.*exp(-1i*z_plus*gamma_p00);
elseif plot_field_component==5 % H_2
    H_2=zeros(length(x_points),length(y_points));
    H_2_0=zeros(length(x_points),length(y_points));
    
    H_2_harmonics3_v=(-1)*R_y_h.*exp(1i*z_plus*diag(full(gamma_1)));
    H_2_harmonics3=reshape(H_2_harmonics3_v',n_x,n_x)';
    
    H_2_harmonics_0=v_y.*exp(-1i*z_plus*gamma_p00);
elseif plot_field_component==6 % H_3
    H_3=zeros(length(x_points),length(y_points));
    H_3_0=zeros(length(x_points),length(y_points));
    
    R_z_h=-(alpha*R_x_h+beta*R_y_h)./diag(gamma_1);
    
    H_3_harmonics3_v=R_z_h.*exp(1i*z_plus*diag(full(gamma_1)));
    H_3_harmonics3=reshape(H_3_harmonics3_v',n_x,n_x)';
    
    H_3_harmonics_0=v_z.*exp(-1i*z_plus*gamma_p00);
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

exp_alpha_y_2_0=zeros(length(y_points),1);
exp_alpha_y_2=zeros(length(y_points),length(beta_n));
exp_alpha_x_2_0=zeros(length(x_points),1);
exp_alpha_x_2=zeros(length(x_points),length(alpha_m));

for y_p=1:1:length(y_points)
    exp_alpha_y_2_0(y_p,1)=exp(-1i*beta_0*1E6*y_points(y_p));
    exp_alpha_y_2(y_p,:)=exp(-1i*beta_n*1E6*y_points(y_p));
end

for x_p=1:1:length(x_points)
    exp_alpha_x_2_0(x_p,1)=exp(-1i*alpha_0*1E6*x_points(x_p));
    exp_alpha_x_2(x_p,:)=exp(-1i*alpha_m*1E6*x_points(x_p));
end

%{
for x_p=1:1:length(x_points)
    for y_p=1:1:length(y_points)
        
        %exp_alpha_x=exp(-1i*2*pi*m_K/Lambda_x*1E6*x_points(x_p));
        %exp_alpha_y=exp(-1i*2*pi*n_K/Lambda_y*1E6*y_points(y_p));
        
        if plot_field_component==1 % E_1
            E_1(x_p,y_p)=sum((E_1_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x_2(x_p,:)');
        elseif plot_field_component==2 % E_2
            E_2(x_p,y_p)=sum((E_2_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x_2(x_p,:)');
        elseif plot_field_component==3 % E_3
            E_3(x_p,y_p)=sum((E_3_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x_2(x_p,:)');
        elseif plot_field_component==4 % H_1
            H_1(x_p,y_p)=sum((H_1_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x_2(x_p,:)');
        elseif plot_field_component==5 % H_2
            H_2(x_p,y_p)=sum((H_2_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x_2(x_p,:)');
        elseif plot_field_component==6 % H_3
            H_3(x_p,y_p)=sum((H_3_harmonics2*exp_alpha_y_2(y_p,:)').*exp_alpha_x_2(x_p,:)');
        end
    end
end
%}

for y_p=1:1:length(y_points)
    
    %{
    if plot_field_component==1 % E_1            
        harmonics2_vector=E_1_harmonics3*exp_alpha_y_2(y_p,:)';
        harmonics2_vector_0=E_1_harmonics_0*exp_alpha_y_2_0(y_p,1)';
    elseif plot_field_component==2 % E_2            
        harmonics2_vector=E_2_harmonics3*exp_alpha_y_2(y_p,:)';
        harmonics2_vector_0=E_2_harmonics_0*exp_alpha_y_2_0(y_p,1)';
    elseif plot_field_component==3 % E_3            
        harmonics2_vector=E_3_harmonics3*exp_alpha_y_2(y_p,:)';
        harmonics2_vector_0=E_3_harmonics_0*exp_alpha_y_2_0(y_p,1)';
    elseif plot_field_component==4 % H_1            
        harmonics2_vector=H_1_harmonics3*exp_alpha_y_2(y_p,:)';
        harmonics2_vector_0=H_1_harmonics_0*exp_alpha_y_2_0(y_p,1)';
    elseif plot_field_component==5 % H_2            
        harmonics2_vector=H_2_harmonics3*exp_alpha_y_2(y_p,:)';
        harmonics2_vector_0=H_2_harmonics_0*exp_alpha_y_2_0(y_p,1)';
    elseif plot_field_component==6 % H_3            
        harmonics2_vector=H_3_harmonics3*exp_alpha_y_2(y_p,:)';
        harmonics2_vector_0=H_3_harmonics_0*exp_alpha_y_2_0(y_p,1)';
    end
    %}
    
    %{
    for x_p=1:1:length(x_points)
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
        E_1(:,y_p)=sum(bsxfun(@times,(E_1_harmonics3*exp_alpha_y_2(y_p,:)'),exp_alpha_x_2'));
        E_1_0(:,y_p)=(bsxfun(@times,(E_1_harmonics_0*exp_alpha_y_2_0(y_p,1)'),exp_alpha_x_2_0'));
        E_1(:,y_p)=E_1(:,y_p)+E_1_0(:,y_p);
    elseif plot_field_component==2 % E_2
        E_2(:,y_p)=sum(bsxfun(@times,(E_2_harmonics3*exp_alpha_y_2(y_p,:)'),exp_alpha_x_2'));
        E_2_0(:,y_p)=(bsxfun(@times,(E_2_harmonics_0*exp_alpha_y_2_0(y_p,1)'),exp_alpha_x_2_0'));
        E_2(:,y_p)=E_2(:,y_p)+E_2_0(:,y_p);
    elseif plot_field_component==3 % E_3
        E_3(:,y_p)=sum(bsxfun(@times,(E_3_harmonics3*exp_alpha_y_2(y_p,:)'),exp_alpha_x_2'));
        E_3_0(:,y_p)=(bsxfun(@times,(E_3_harmonics_0*exp_alpha_y_2_0(y_p,1)'),exp_alpha_x_2_0'));
        E_3(:,y_p)=E_3(:,y_p)+E_3_0(:,y_p);
    elseif plot_field_component==4 % H_1
        H_1(:,y_p)=sum(bsxfun(@times,(H_1_harmonics3*exp_alpha_y_2(y_p,:)'),exp_alpha_x_2'));
        H_1_0(:,y_p)=(bsxfun(@times,(H_1_harmonics_0*exp_alpha_y_2_0(y_p,1)'),exp_alpha_x_2_0'));
        H_1(:,y_p)=H_1(:,y_p)+H_1_0(:,y_p);
    elseif plot_field_component==5 % H_2
        H_2(:,y_p)=sum(bsxfun(@times,(H_2_harmonics3*exp_alpha_y_2(y_p,:)'),exp_alpha_x_2'));
        H_2_0(:,y_p)=(bsxfun(@times,(H_2_harmonics_0*exp_alpha_y_2_0(y_p,1)'),exp_alpha_x_2_0'));
        H_2(:,y_p)=H_2(:,y_p)+H_2_0(:,y_p);
    elseif plot_field_component==6 % H_3
        H_3(:,y_p)=sum(bsxfun(@times,(H_3_harmonics3*exp_alpha_y_2(y_p,:)'),exp_alpha_x_2'));
        H_3_0(:,y_p)=(bsxfun(@times,(H_3_harmonics_0*exp_alpha_y_2_0(y_p,1)'),exp_alpha_x_2_0'));
        H_3(:,y_p)=H_3(:,y_p)+H_3_0(:,y_p);
    end
    %%}
end

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
    profile_mode=zobrazeni_pole(-H_3,field_function); 
end

%clear exp_alpha_y_2_0 exp_alpha_y_2 exp_alpha_x_2_0 exp_alpha_x_2

%{
figure
set(gcf, 'Renderer', 'zbuffer');
axes('FontSize',13)
surf(1E6.*y_points,1E6.*x_points,profile_mode,'LineWidth',1);

shading flat
axis([coordinate_y_1(1) coordinate_y_1(end) coordinate_x_1(1) coordinate_x_1(end)]);
ylabel('x [\mu{}m]','FontSize',13);
xlabel('y [\mu{}m]','FontSize',13);
%view([90 90]);

if set_waveguide_framework==1
hold on;
waveguide_framework;
end

view([90 90])
%view([0 -90]);

%Mov(cylus_frame) = getframe;
%cylus_frame=cylus_frame+1;
%end

%movie(Mov,1,8)
%movie2avi(Mov, '1_2.avi','fps',4)
%}