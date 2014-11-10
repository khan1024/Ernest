%clear all

number_of_layers=1; % this parameter defines number of layers

Lambda_x=0.3; %  [um]
Lambda_y=0.3; %  [um]
lambda=0.5;
theta0=60;

n_1=1;
%n_3=sqrt(2.138);
n_3=sqrt(-2.865-3.211*1i);

resolution_dx=512;

r_cylinder=0.1; % radius of the cylinder [um]
t_cylinder=0.1; % thickness of the cylinder [um]

n_cylinder=1; % refractive index of the cylinder
%n_layer=sqrt(2.138);%-2.865-3.211*1i); % refractive index of the surrounding layer
n_layer=sqrt(-2.865-3.211*1i);

show_grating_mesh=2; % 1-yes, 2-no

addpath('structures_scripts')
[p_x,p_y,mesh_n_mat]=rectangular_grid(Lambda_x,Lambda_y,resolution_dx,r_cylinder,n_cylinder,n_layer,show_grating_mesh);
%--------------------------------------------------------- 

% definition of layer number 1 
l=1;

nvm_method(l)=1; % 1=normal vector method is on (only in this layer)

factorization_method=1; % 1-radial, 2-complex radial, 3-
[N_x_input,N_y_input]=rectangular_grid_nvm(Lambda_x,Lambda_y,resolution_dx,factorization_method,r_cylinder);


thickness(l)=t_cylinder; % thickness of the layer

coordinate_x_1=p_x;
coordinate_y_1=p_y;

r_index_1=mesh_n_mat; % matrix of the refrative index grid 


