% !! For RCWA calculation: comment “clear all” and set plot_grating=2;

%clear all
plot_grating=2;

% 1-add homogeneous layer
% 2-add polygon layer
% 3-add sinus layer
% 4-add user's layer

n1=1;
n3=2;

Lambda=1;


layer_1=add_sin_layer(Lambda,Lambda,0.5,2,1,1); % sinus layer


create_grating;

% add_hom_layer(Lambda,layer_thickness,ref_index)
% add_sin_layer(Lambda,Lambda_object,layer_thickness,ref_index_in,ref_index_out,n_layers)
% add_tri_layer(Lambda,Lambda_object,layer_thickness,ref_index_in,ref_index_out,n_layers)
% add_pol_layer(Lambda,x,y,ref_index_in,ref_index_out,n_layers)

% create_grating; % script creates the input