%{
clear all

Lambda_x=0.3;
Lambda_y=Lambda_x*sqrt(3);

dx=0.01;
r_cylinder=0.1;
eps_cylinder=1;
eps_layer=3;

show_grating_mesh=1; % 1-yes, 2-no
%}

% nutnost přenastavit Lambda_x a Lambda_y podle aktualni rozliseni

%{
Lambda_x_old=Lambda_x;
Lambda_y_old=Lambda_y;

% 1) Lambda_x
fix_Lambda_x=fix(Lambda_x/dx);
Lambda_x=fix_Lambda_x*dx;

% 1) Lambda_y
fix_Lambda_y=fix(Lambda_y/dx);
Lambda_y=fix_Lambda_y*dx;
%}

function [p_x,p_y,mesh_n_mat]=cylinder_cover(Lambda_x,Lambda_y,resolution_dx,r_cylinder,r_cylinder_cover,n_cylinder,n_cover,n_layer,show_grating_mesh)

save_as='rectangular_grating_mat';

p_x=linspace(0,Lambda_x,resolution_dx+1);
p_y=linspace(0,Lambda_y,resolution_dx+1);

%p_x=0:dx:Lambda_x;
%p_y=0:dx:Lambda_y;

mesh_n_mat=zeros(length(p_y)-1,length(p_x)-1)+n_layer;
pixel_center=zeros(length(p_y)-1,length(p_x)-1,2); % 1-x coordinatex, 2-y coordinate

% nalézt středy jednotlivých pixelů
for i=1:1:length(p_y)-1
    for j=1:1:length(p_x)-1
        pixel_center(i,j,1)=(p_x(j)+p_x(j+1))/2;
        pixel_center(i,j,2)=(p_y(i)+p_y(i+1))/2;
    end
end

% pokud střed pixelu bude uvnitř nejakého prvku, permitivita pixelu se
% označí podle permitivity prvku
object(:,:,1)=[Lambda_x/2,Lambda_y/2,r_cylinder_cover,n_cover];
object(:,:,2)=[Lambda_x/2,Lambda_y/2,r_cylinder,n_cylinder];

[m,n,number_of_objects] = size(object);

for i=1:1:length(p_y)-1
    for j=1:1:length(p_x)-1
        
        x_coordinate=pixel_center(i,j,1);
        y_coordinate=pixel_center(i,j,2);
        
        for k=1:1:number_of_objects
            
            distance_origin_pixel=sqrt((x_coordinate-object(1,1,k))^2+(y_coordinate-object(1,2,k))^2);
            
            if distance_origin_pixel<=object(1,3,k)
                mesh_n_mat(i,j)=object(1,4,k);
            end
        
        end
        
    end
end
%

size_eps=size(mesh_n_mat);

n_of_rows=size_eps(1);
n_of_columns=size_eps(2);

if show_grating_mesh==1
    grating_mesh;
end

p_x_mat=p_x;
p_y_mat=p_y;

%{
savefile = [save_as '.mat'];
if factorization_method==0
    save(savefile, 'p_x_mat', 'p_y_mat', 'mesh_n_mat');
else
    save(savefile, 'p_x_mat', 'p_y_mat', 'mesh_n_mat', 'N_x_input', 'N_y_input');
end
%}