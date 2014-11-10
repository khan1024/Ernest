%Test case from J. Opt. Soc. Am. A/Vol. 14, No. 10/Oct. 1997, pp. 2758-2767, Example 3

% N_X=N_Y=10;
%order  [ 0.00,  0.00] R_0=0.3006980687 T_0=0.0129889964 sumR=0.5527021551 sumT=0.1252458023 R_90=0.5430724145 T_90=0.0178547750 sumR=0.6373546277 sumT=0.1153169719 R_psi=0.3006980687 T_psi=0.0129889964 sumR=0.5527021551 sumT=0.1252458023

number_of_layers=1;

n_1=1;
n_3=1.5;
lambda=1;
theta0=30;
phi0=30;
zeta0=30;

l=1; % layer 1
thickness(l)=1;

coordinate_x_1=2.*[0,0.25,0.75,1];
coordinate_y_1=2.*[0,0.25,0.75,1];

nb=1;
na=1-5*1i; % correct value according to [2] Li, J. Opt. A: Pure Appl. Opt., Vol.5, 2003 - Fourier modal method for crossed
% anisotropic gratings with arbitrary permittivity and permeability tensors

r_index_1=[na, na, na;
    na, nb, na;
    na, na, na];

