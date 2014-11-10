% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

clear all
format long;

addpath('gratings');
addpath('model');
addpath('lib');
tic;

%---------------------------------------------------
% basic setup of the RCWA algorithm
number_of_orders=41; % odd number
matrix_algorithm=1; % 1-S-matrix, 2-T-matrix
change_matrix_base=1; % 1-s,p base, 2-x,y,z base
%---------------------------------------------------
% setup of the incident wave
lambda=0.4; % wavelength (um)
use_dispersion=1; % 1-yes, 2-no, see file setup_dispersion.m setup_dispersion.m
theta0=75; % incident angle [degree]
polarization=1; % 0-conical diffraction, 1-TM polarization, 2-TE polarization
faktorization=1; % 1-correct factorization

phi0=0; % conical angle [degree], 0=planar diffraction
psi0=0; % polarization angle [degree], 0-TM polarization, 90-TE polarization

%---------------------------------------------------
% setup of the plotted field
plot_field_region_1=2; %1-yes, 2-no
plot_field_region_2=2; %1-yes, 2-no
plot_field_region_3=2; %1-yes, 2-no
number_of_plotted_period=1; % how many periods
size_region_1=1; % times region 2 thickness
size_region_3=1; % times region 2 thickness
resolution_x=200; % x resolution of figure
resolution_z=200; % x resolution of figure
slice=2; % 1-yes, 2-no, slice only for grating region 2, E field=E/377
slice_position=0.5; % times region 2 thickness
slice2=2; % 1-yes, 2-no, only for grating region 2, E field=E/377
slice2_position=0.5; % times region 2 thickness
save_field=2; % 1 --- yes, 2 ---no
save_field_as='test';

% planar_field_component:
%    TM polarization TE polarization
%  1=  E_x			H_x
%  2=  H_y			E_y
%  3=  E_z			H_z
field_component_operation=1; % 0-(abs)^2, 1-abs, 2-angle, 3-real, 4-imag
planar_field_component=1;
conical_field_component=2; %1-H_x, 2-H_y, 3-H_z, 4-E_x, 5-E_y, 6-E_z
set_view_grating_boundaries=1; % 1-yes, 2-no, you can see grating boundaries when field is plotted

get_electric_intensity=2; % 1-yes, 2-no
get_local_absoption=2; % 1-yes, 2-no
get_poynting_vector=2; % 1-yes, 2-no
%---------------------------------------------------
% Setup of the grating
%  grating=
%  0 --- binary grating
%  1 --- blazed grating, type 1
%  2 --- blazed grating, typy 2
%  3 --- sinusoidal grating
%  4 --- triangular grating
%  5 --- parabolic grating
%  6 --- general triangular grating
%  7 --- trapezoidal grating
%  8 --- tetragonal grating
%  9 --- Gauss grating
%  10 --- arbitrary grating from file

grating=0; % switch grating
Lambda=.5; % grating period [um] (exception grating=10)
thickness_total=0.5; % total thickness [um] (exception grating=10)
plot_grating=2; % 1--yes, 2--no

% refractive indicies (sign for metals n_r-1i*n_i)
%---------------------------------------------------
n1=1; % homogenous incident region 1
% n3=0.22-6.71*1i; % homogenous region 3
n3=2;
ng=1; % Groove indicies, disabled when grating=10
nr=n3; % Ridge indices, disabled when grating=10
%---------------------------------------------------
see_diffraction_layout=2; % 1--yes, 2--no; view diffraction orders layout - planar diffraction
%---------------------------------------------------
switch grating
    case 0 % binary grating
        duty_cycle=.3; 
        shift=.5; 
    case 1 % blazed grating, type 1 <
        number_of_layers=20;
    case 2 % blazed grating, type 2 >
        number_of_layers=20;
    case 3 % sinusoidal grating
        number_of_layers=20;
    case 4 % triangular grating
        number_of_layers=20;
    case 5 % parabolic grating
        number_of_layers=10;
    case 6 % general triangular grating
        number_of_layers=20;
        a=.65; % distance of the triangular top, from left corner 0<a<1
        b=.7; % size of the triangular base 0<b<1
    case 7 % trapezoidal grating
        number_of_layers=20; 
        a_1=.3; % distance of the left top edge, from left corner
        a_2=.5; %  distance of the right top edge, from right corner, a_1+a_2<1
        b=.5; % size of the trapezoidal base
    case 8 % tetragonal grating
        number_of_layers=20;
        a_1=.25; % horizontal position of the left top edge
        a_2=.1; % horizontal position of the right top edge
        c_1=.5; % vertical position of the left top edge
        c_2=1; % vertical position of the right top edge
        b=.7; % size of the tetragonal base
    case 9 % Gaussian grating
        number_of_layers=20;
    case 10 % arbitrary grating from file
        grating_relative_dimensions=2; % 1-relative size, 2-absolute size
        
        input_grating_file='met_binary.m';
end

%-------------------------------------------------------     
% setup measurement
% template measurement
% 0 --- shows all diffraction efficiencies
% 1 --- dependence of the diffraction efficiencies on the number of orders
% 2 --- dependence of the diffraction efficiencies on the wavelength
% 3 --- dependence of the diffraction efficiencies on the total thickness
% 4 --- dependence of the diffraction efficiencies on the refractive index $n_1$ ($n_g$)
% 5 --- dependence of the diffraction efficiencies on the refractive index $n_3$ ($n_r$)
% 6 --- dependence of the diffraction efficiencies on the incident angle theta
% 7 --- dependence of the diffraction efficiencies on the number of layers
% 8 --- dependence of the diffraction efficiencies on the grating period Lambda
% 9 --- dependence of the diffraction efficiencies on the binary duty cycle
% 10 --- dependence of the diffraction efficiencies on the angle psi (only conical diffraction)
% 11 --- dependence of the diffraction efficiencies on the conical angle phi (only conical diffraction)
%D \stopitemize

measurement=0; % measurement 1 --- check studying_order

% !!! setup of measured diff. order measurement>0

diffraction_efficiencies_c=1; %1--D_R, 2--D_T
studying_order=0; %check minimum_number_of_orders
save_m=2; % 1 --- yes, 2 ---no
save_as='jin';

switch measurement
    case 0  % binary grating
        duty_cycle=0.3; 
        shift=.5;      
    case 1
        minimum_number_of_orders=3; % only odd number
        maximum_number_of_orders=101; % only odd number
        step=2; % step=2 
    case 2 %
        minimum_wavelength=1.24; %um
        maximum_wavelength=4.13; %um
        step=0.1;   
    case 3 % 
        minimum_thickness=0.05; %um
        maximum_thickness=1; %um
        step=0.005;
    case 4
        minimum_index=1; 
        maximum_index=3; 
        step=0.1;
        change_of_index=2; % 1 --- n1=ng, 2 --- n1 is changing, 3 --- ng is changing
    case 5
        minimum_index=1.4; 
        maximum_index=1.6; 
        step=.01;
        change_of_index=2; % 1 --- n3=nr, 2 --- n3 is changing, 3 --- nr is changing
    case 6
        minimum_angle=0; % theta [degree] 
        maximum_angle=20; % theta [degree]
        step=0.1;
    case 7
        minimum_layer=50; 
        maximum_layer=80;
        step=5;
    case 8
        minimum_Lambda=0.05; % um
        maximum_Lambda=5; % um
        step=0.01;
    case 9
        minimum_duty_cycle=0.02; % um
        maximum_duty_cycle=0.98; % um
        step=0.00096;
    case 10
        minimum_angle_psi=0; % [degree] 
        maximum_angle_psi=90; % [degree]
        step=0.5;        
    case 11
        minimum_angle_phi=0; % [degree]
        maximum_angle_phi=89; % [degree]
        step=0.5;
        
    otherwise
        error('it is not defined')
        
end
help_control_file;
save_field_function; % for field calculation
string_time=sprintf('time of calculation = %e s',toc);
disp(string_time);
disp(diff_R)
% beep;