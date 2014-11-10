size_field=size(intensity_E);
size_x=size_field(1);
size_z=size_field(2);
new_z=0:thickness/size_z:thickness-thickness/size_z;

figure;
set(gcf, 'Renderer', 'zbuffer');
surf(new_z.*1E6,x.*1E6,intensity_E); 
axis([0 thickness*1E6 x_min*1E6 x_max*1E6]);
title('electric intensity');
xlabel('z [\mu{}m]','FontSize',14);
ylabel('x [\mu{}m]','FontSize',14);
shading flat;

if in_octave==1 && old_field_version==1
    view_octave([90 90]);
else
    view([90 90]);    
end