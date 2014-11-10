% minus sign because imag(eps)<0
local_absoption_matrix=-(k_0^2/k_1_z((number_of_orders+1)/2)*z_res_/resolution_x).*bsxfun(@times,intensity_E,imag(epsilon_matrix_local/epsilon));

global_absorption=1-sum(D_R)-sum(D_T);

if polarization==1
    local_absoption_matrix=local_absoption_matrix/(sqrt(mu/epsilon)/n1)^2;
elseif polarization==2
    % ok
end

z_min=1;
z_max=absorption_in_layer_thickness(1);
layer_absorption=zeros(number_of_layers,1);
for i=1:1:number_of_layers
    layer_absorption(i)=sum(sum(local_absoption_matrix(:,z_min:z_max)));
    if i<number_of_layers
    z_min=z_min+absorption_in_layer_thickness(i);
    z_max=z_max+absorption_in_layer_thickness(i+1);
    end
end


integrated_absoption=sum(sum(local_absoption_matrix));
norm_constant=global_absorption/integrated_absoption;

integrated_absoption
global_absorption

size_field=size(intensity_E);
size_x=size_field(1);
size_z=size_field(2);
new_z=0:thickness/size_z:thickness-thickness/size_z;

figure;
set(gcf, 'Renderer', 'zbuffer');
surf(new_z.*1E6,x.*1E6,local_absoption_matrix); 
axis([0 thickness*1E6 x_min*1E6 x_max*1E6]);
title('local absoption');
xlabel('z [\mu{}m]','FontSize',14);
ylabel('x [\mu{}m]','FontSize',14);
shading flat;

if in_octave==1 && old_field_version==1
    view_octave([90 90]);
else
    view([90 90]);    
end
