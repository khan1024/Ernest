% xy_slice

plot_field_component=2; % 1-E_x, 2-E_y, 3-E_z, 4-H_x, 5-H_y, 6-H_z
field_function=1; % 1-abs, 2-angle, 3-real, 4-imag
set_waveguide_framework=2; % 1-yes, 2-no

xy_slice_coordinate_z=0.5;
resolution_x=200;
resolution_y=200;

x_point_b=coordinate_x_1(1);
x_point_e=coordinate_x_1(end);

y_point_b=coordinate_y_1(1);
y_point_e=coordinate_y_1(end);

%-------------------------------------------------------------------------
xy_slice_view=1;
if xy_slice_coordinate_z<=0
    xy_slice_1; % input region
elseif xy_slice_coordinate_z>0 && xy_slice_coordinate_z<sum(thickness)
    xy_slice_2; % grating region
elseif xy_slice_coordinate_z>=sum(thickness)
    xy_slice_3; % output region
end
xy_slice_view=2;

figure
set(gcf, 'Renderer', 'zbuffer');
axes('FontSize',13)
surf(1E6.*y_points,1E6.*x_points,profile_mode,'LineWidth',1);

shading flat
axis([y_point_b y_point_e x_point_b x_point_e]);
ylabel('x [\mu{}m]','FontSize',13);
xlabel('y [\mu{}m]','FontSize',13);
%view([90 90]);

if set_waveguide_framework==1
hold on;
waveguide_framework;
end

view([90 90])