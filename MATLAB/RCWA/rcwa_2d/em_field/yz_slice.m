% yz slice

plot_field_component=2; % 1-E_x, 2-E_y, 3-E_z, 4-H_x, 5-H_y, 6-H_z
field_function=1; % 1-abs, 2-angle, 3-real, 4-imag
set_waveguide_framework=2; % 1-yes, 2-no


yz_slice_coordinate_x=0.6;
resolution_y=400;
resolution_z=400;


y_point_b=0.0;
y_point_e=1.2;

z_point_b=-1;
z_point_e=2;

%-------------------------------------------------------------------------

z_points_b_e=z_point_b:(z_point_e-z_point_b)/resolution_z:z_point_e;
y_points=y_point_b*1E-6:((y_point_e-y_point_b)*1E-6/resolution_y):(y_point_e*1E-6);

z_points_input=find(z_points_b_e<=0);
z_points_grating=find(z_points_b_e>0 & z_points_b_e<sum(thickness));
z_points_output=find(z_points_b_e>=sum(thickness));

em_field=[];
if isempty(z_points_input)==0
    z_points_1=z_points_b_e(z_points_input(1):z_points_input(end));
    yz_slice_1;
    em_field=[em_field;profile_mode1];
else
    profile_mode1=[];
    em_field=[em_field;profile_mode1];
end

yz_slice_view=1;
if isempty(z_points_grating)==0
    z_points_2=z_points_b_e(z_points_grating(1):z_points_grating(end));
    yz_slice_2;
    em_field=[em_field;profile_mode2];
else
    profile_mode2=[];
    em_field=[em_field;profile_mode2];
end
yz_slice_view=2;

if isempty(z_points_output)==0
    z_points_3=z_points_b_e(z_points_output(1):z_points_output(end));
    yz_slice_3;
    em_field=[em_field;profile_mode3];
else
    profile_mode3=[];
    em_field=[em_field;profile_mode3];
end


figure
set(gcf, 'Renderer', 'zbuffer');
surf(1E6.*y_points,z_points_b_e,em_field); shading flat
xlabel('y [\mu{}m]','FontSize',13);
ylabel('z [\mu{}m]','FontSize',13);

axis([y_point_b y_point_e z_point_b z_point_e]);

view([0 -90])