% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% !!!! to speed-up this code compile (directory /mex_files/source): asr_fourier_2D_mex.c, 
% !!!! asr_fourier_2D_split_bt_mex.c, asr_fourier_2D_split_tb_mex.c, and
% !!!! alpha_beta_mex.c; and place these library into the /mex-files directory
% !!!! Matlab: mex filename.c, Octave: mkoctfile --mex filename.c
% !!!! currently, there are mex-files for Linux 32/64-bit and Windows 32/64-bit

clear all
format long;
time1=tic;
addpath('lib');
addpath('model');
addpath('structures');
addpath('mex_files');

%---------------------------------------------------
N_X=10; % N_X=20 takes approx. 8 GB of memory
N_Y=N_X;
use_fftw=2;
%---------------------------------------------------
Lambda_x=1; %  [um] only for grating 1 and 2
Lambda_y=Lambda_x; %  [um] only for grating 1 and 2
n_1=1.0; % superstrate
n_3=1.5; % substrate
mu=1;
%---------------------------------------------------
use_dispersion=2; % 1-yes, 2-no, see file setup_dispersion.m setup_dispersion.m
lambda=0.58; % [um]
theta0=0; % incident angle [degree], if theta=0 -> theta0+1E-10
phi0=0; % conical angle [degree]
psi0=90;  % polarization angle [degree], 0-TM polarization, 90-TE polarization

% switch grating
% 0 - grating from file, Li's factorization rules, there is a possibility
% to create case 1 and case 2 (zig-zag method)
% 1 - eliptic cylinder, rectangular grating, Lalanne's factorization rules (parameter_alpha),
% calculation of metallic grating is not very precise
% 2 - eliptic cylinder, hexagonal grating, Lalanne's factorization rules (parameter_alpha), 
% calculation of metallic grating is not very precise

grating=0;
switch grating
    case 0
        Lambda_y=0.55;
        Lambda_x=1.5*0.55;
        r_cylinder=0.04;
        setup_dispersion;
        open_grating_file='rectangular_grating_sensor_B_d.m';
        optimalization_grating_0=1; % 1-yes, 2-no
    case 1 % eliptic cylinder, rectangular grating
        number_of_layers=1;
        setup_dispersion;
        n_layer=1.5;
        n_cylinder=1;
        thickness=1; % [um]
        a_ellipse=0.25; % semi-axis [um]
        b_ellipse=0.25; % semi-axis [um]
    case 2 % eliptic cylinder, hexagonal grating, set Lambda_y=Lambda_x*sqrt(3)
        number_of_layers=1;
        setup_dispersion;
        n_layer=1;
        n_cylinder=2;
        thickness=0.5; % [um]
        a_ellipse=0.1; % semi-axis [um]
        b_ellipse=a_ellipse; % semi-axis [um]
end

parameter_alpha=0.5; % in range (0-1) alpha=1 classic calculation (without factorization)

% setup measurement
% template measurement
% 0 --- shows diffraction efficiency
% 1 --- dependence of the diffraction efficiencies on the number of orders
% 2 --- dependence of the diffraction efficiencies on the wavelength
% 3 --- dependence of the diffraction efficiencies on the grating thickness
% 4 --- dependence of the diffraction efficiencies on the incident angle theta
% 5 --- dependence of the diffraction efficiencies on the angle psi
% 6 --- dependence of the diffraction efficiencies on the conical angle phi

measurement=0;

diffraction_efficiency_order=[0,0]; 
save_m=0; % 1 --- yes, 2 ---no
save_m_as='test';

switch measurement
    case 0 % diffraction_efficiencies       
    case 1
        minimum_N_X_N_Y=1; %
        maximum_N_X_N_Y=15; % max. 20 (it takes 8 GB)
        step=1;
    case 2 
        minimum_wavelength=0.632; % [um]
        maximum_wavelength=0.633; % [um]
        step=0.0005;   
    case 3 % 
        minimum_thickness=0.4; % [um]
        maximum_thickness=1.0; % [um]
        step=0.02;
    case 4
        minimum_angle_theta=0; % [degree]
        maximum_angle_theta=89; % [degree]
        step=0.5;
    case 5
        minimum_angle_psi=0; % [degree] 
        maximum_angle_psi=90; % [degree]
        step=0.5;        
    case 6
        minimum_angle_phi=0; % [degree]
        maximum_angle_phi=89; % [degree]
        step=0.5;
    otherwise
        error('it is not defined')
end
help_control_file;
string_time=sprintf('\n time of calculation = %d s',toc(time1));
disp(string_time);
beep;