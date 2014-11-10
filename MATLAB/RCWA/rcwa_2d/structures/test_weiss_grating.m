%clear all

number_of_layers=1; % this parameter defines number of layers

Lambda_x=0.7; %  [um]
Lambda_y=0.7; %  [um]

dx=0.001; % discretization step [um] 

r_cylinder=0.15; % radius of the cylinder [um]
t_cylinder=0.05; % thickness of the cylinder [um]

n_cylinder=n_weiss_drude; % refractive index of the cylinder
n_layer=1; % refractive index of the surrounding layer

show_grating_mesh=2; % 1-yes, 2-no

%factorization_method=1; % 1-Li, 2-moje radial, 3-complex radial, 4-
%---------------------------------------------------------
addpath('structures_scripts')

setup_rectangular_grating; % it creates the discretization
load rectangular_grating_mat

%--------------------------------------------------------- 
% definition of layer number 1 
l=1;

thickness(l)=t_cylinder; % thickness of the layer

px(:,:,l)=p_x;
py(:,:,l)=p_y;

n_mat(:,:,l)=mesh_n_mat; % matrix of the refrative index grid 


