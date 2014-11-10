clear all

number_of_layers=2; % this parameter defines number of layers

Lambda_x=.5; %  [um]
Lambda_y=.5; %  [um]

resolution_dx=256; 

r_cylinder=0.2; % radius of the cylinder [um]
t_cylinder=.1; % thickness of the cylinder [um]

n_cylinder=6; % refractive index of the cylinder
n_layer=1.33; % refractive index of the surrounding layer

r_cylinder_cover=r_cylinder+0.01;
n_cover=1.5;

show_grating_mesh=1; % 1-yes, 2-no

%---------------------------------------------------------
% generation of refractive index mesh
addpath('structures_scripts')
%[p_x,p_y,mesh_n_mat]=rectangular_grid(Lambda_x,Lambda_y,resolution_dx,r_cylinder,n_cylinder,n_layer,show_grating_mesh);
[p_x,p_y,mesh_n_mat]=cylinder_cover(Lambda_x,Lambda_y,resolution_dx,r_cylinder,r_cylinder_cover,n_cylinder,n_cover,n_layer,show_grating_mesh);
%--------------------------------------------------------- 

%{
% definition of layer number 1 
l=1;

nvm_method(l)=1; % 1=normal vector method is on (only in this layer)
factorization_method=2; % 1-radial, 2-complex radial, 3-
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
%}
