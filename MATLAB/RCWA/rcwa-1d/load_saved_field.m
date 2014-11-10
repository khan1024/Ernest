% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

clear all
addpath('lib');
addpath('measurements');
load_file='test';

plot_field_region_1=2; %1-yes, 2-no
plot_field_region_2=2; %1-yes, 2-no
plot_field_region_3=2; %1-yes, 2-no

field_component_operation=3;  % 1-abs, 2-angle, 3-real, 4-imag
set_view_grating_boundaries=1; % 1-yes, 2-no, you can see grating boundaries when field is plotted

slice=2; % 1-yes, 2-no, only for grating region 2, E field=E/377
slice_position=0.5; % times region 2 thickness
slice2=2; % 1-yes, 2-no, only for grating region 2, E field=E/377
slice2_position=0.5; % times region 2 thickness

%--------------------------------------
temp_plot_field_region_1=plot_field_region_1;
temp_plot_field_region_2=plot_field_region_2;
temp_plot_field_region_3=plot_field_region_3;

filename=[load_file '.mat'];
load(filename);

if plot_field_region_1<=temp_plot_field_region_1 || plot_field_region_2<=temp_plot_field_region_2 || plot_field_region_3<=temp_plot_field_region_3
    % OK
else
    error('current region is not available');
end

plot_field_region_1=temp_plot_field_region_1;
plot_field_region_2=temp_plot_field_region_2;
plot_field_region_3=temp_plot_field_region_3;

new_figure=1;
save_field=2;
if polarization==0
    field_graphs_conical;
else
    field_graphs;
end