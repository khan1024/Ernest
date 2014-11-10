% xz slice

plot_field_component=2; % 1-E_x, 2-E_y, 3-E_z, 4-H_x, 5-H_y, 6-H_z
field_function=1; % 1-abs, 2-angle, 3-real, 4-imag
set_waveguide_framework=2; % 1-yes, 2-no


xz_slice_coordinate_y=0.9;
resolution_x=400;
resolution_z=400;

x_point_b=0.0;%coordinate_x_1(1);
x_point_e=1.2;%coordinate_x_1(end);

z_point_b=-1;
z_point_e=2;

%-------------------------------------------------------------------------

z_points_b_e=z_point_b:(z_point_e-z_point_b)/resolution_z:z_point_e;
x_points=x_point_b*1E-6:((x_point_e-x_point_b)*1E-6/resolution_x):(x_point_e*1E-6);

z_points_input=find(z_points_b_e<=0);
z_points_grating=find(z_points_b_e>0 & z_points_b_e<sum(thickness));
z_points_output=find(z_points_b_e>=sum(thickness));

em_field=[];
if isempty(z_points_input)==0
    z_points_1=z_points_b_e(z_points_input(1):z_points_input(end));
    xz_slice_1;
    em_field=[em_field,profile_mode1];
else
    profile_mode1=[];
    em_field=[em_field,profile_mode1];
end

xz_slice_view=1;
if isempty(z_points_grating)==0
    z_points_2=z_points_b_e(z_points_grating(1):z_points_grating(end));
    xz_slice_2;
    em_field=[em_field,profile_mode2];
else
    profile_mode2=[];
    em_field=[em_field,profile_mode2];
end
xz_slice_view=2;

if isempty(z_points_output)==0
    z_points_3=z_points_b_e(z_points_output(1):z_points_output(end));
    xz_slice_3;
    em_field=[em_field,profile_mode3];
else
    profile_mode3=[];
    em_field=[em_field,profile_mode3];
end

figure
set(gcf, 'Renderer', 'zbuffer');
surf(z_points_b_e,1E6.*x_points,em_field); shading flat
ylabel('x [\mu{}m]','FontSize',13);
xlabel('z [\mu{}m]','FontSize',13);

axis([z_point_b z_point_e x_point_b x_point_e]);

view([90 90])