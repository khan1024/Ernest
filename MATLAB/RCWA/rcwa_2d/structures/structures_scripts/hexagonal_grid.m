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
function [p_x,p_y,mesh_n_mat]=hexagonal_grid(Lambda_x,Lambda_y,resolution_dx,r_cylinder,n_cylinder,n_layer,show_grating_mesh)

dx=Lambda_x/resolution_dx;

Lambda_x_old=Lambda_x;
Lambda_y_old=Lambda_y;

% 1) Lambda_x
fix_Lambda_x=fix(Lambda_x/dx);
Lambda_x=fix_Lambda_x*dx;

% 1) Lambda_y
fix_Lambda_y=fix(Lambda_y/dx);
Lambda_y=fix_Lambda_y*dx;

object(:,:,1)=[Lambda_x/2,Lambda_y/2,r_cylinder,n_cylinder];
object(:,:,2)=[0,0,r_cylinder,n_cylinder];
object(:,:,3)=[0,Lambda_y,r_cylinder,n_cylinder];
object(:,:,4)=[Lambda_x,0,r_cylinder,n_cylinder];
object(:,:,5)=[Lambda_x,Lambda_y,r_cylinder,n_cylinder];

save_as='hexagonal_grating_mat';

p_x=0:dx:Lambda_x;
p_y=0:dx:Lambda_y;

mesh_n_mat=zeros(length(p_y)-1,length(p_x)-1)+n_layer;
pixel_center=zeros(length(p_y)-1,length(p_x)-1,2); % 1-x coordinatex, 2-y coordinate

center_px=(p_x(1:end-1)+p_x(2:end))/2;
center_py=(p_y(1:end-1)+p_y(2:end))/2;

for i=1:1:length(p_y)-1
    pixel_center(i,:,1)=center_px;
end

for j=1:1:length(p_x)-1
    pixel_center(:,j,2)=center_py;
end

%{
% nalézt středy jednotlivých pixelů
for i=1:1:length(p_y)-1
    for j=1:1:length(p_x)-1
        pixel_center(i,j,1)=(p_x(j)+p_x(j+1))/2;
        pixel_center(i,j,2)=(p_y(i)+p_y(i+1))/2;
    end
end
%}

% pokud střed pixelu bude uvnitř nejakého prvku, permitivita pixelu se
% označí podle permitivity prvku

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


end

%savefile = [save_as '.mat'];
%save(savefile, 'p_x_mat', 'p_y_mat', 'mesh_n_mat');