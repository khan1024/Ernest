%clear all

number_of_layers=2; % this parameter defines number of layers

Lambda_x=1; %  [um]
Lambda_y=1; %  [um]

resolution_dx=256; 

r_cylinder=0.25; % radius of the cylinder [um]
t_cylinder=1; % thickness of the cylinder [um]

n_cylinder=1; % refractive index of the cylinder
n_layer=1.5; % refractive index of the surrounding layer

show_grating_mesh=2; % 1-yes, 2-no

%---------------------------------------------------------
% generation of refractive index mesh
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

%--------------------------------------------------------- 
% definition of layer number 2 
l=2;

thickness(l)=2; % thickness of the layer

coordinate_x_2=[0,Lambda_x];
coordinate_y_2=[0,Lambda_y];

r_index_2=1.5; % matrix of the refrative index grid 

