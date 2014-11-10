%clear all

number_of_layers=1;

Lambda_x=1;
Lambda_y=Lambda_x*sqrt(3);
resolution_dx=512;

r_cylinder=0.25; % radius
t_cylinder=0.5; % thickness

n_cylinder=sqrt(2.56);
n_layer=1;

show_grating_mesh=2; % 1-yes, 2-no

%---------------------------------------------------------
% generation of refractive index mesh
addpath('structures_scripts')
[p_x,p_y,mesh_n_mat]=hexagonal_grid(Lambda_x,Lambda_y,resolution_dx,r_cylinder,n_cylinder,n_layer,show_grating_mesh);
%--------------------------------------------------------- 
% definition of layer number 1 
l=1;
nvm_method(l)=1; % 1=normal vector method is on
factorization_method=1; % 0-none (Li), 1-radial
[N_x_input,N_y_input]=hexagonal_grid_nvm(Lambda_x,Lambda_y,resolution_dx,factorization_method);

thickness(l)=t_cylinder;

coordinate_x_1=p_x;
coordinate_y_1=p_y;

r_index_1=mesh_n_mat;

%{
% definition of layer number 2 
l=2;

thickness(l)=t_cylinder;

coordinate_x_2=[0,p_x(end)];
coordinate_y_2=[0,p_y(end)];

r_index_2=1.5;
%}