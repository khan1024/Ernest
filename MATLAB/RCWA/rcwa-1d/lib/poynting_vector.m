% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic


plot_field_region_2=1;
if polarization==1                    
    for i=1:1:3
        planar_field_component=i;
        if i==1
            try
                field2_6;
                old_field_version=2;
            catch
                sprintf('switch to older version');
                old_field_version=1;
                field2_5;
            end                
            component_E_x=field_E_2;                
            clear field_E_2
        elseif i==2
            if old_field_version==1
                field2_5;
            else
                field2_6;
            end
            component_H_y=field_H_2;
            clear field_H_2
        elseif i==3    
            if old_field_version==1
                field2_5;
            else
                field2_6;
            end
            component_E_z=field_E_2;    
            clear field_E_2
        end
    end
elseif polarization==2
    for i=1:1:3
        planar_field_component=i;
        if i==1  
            try
                field2_6;
                old_field_version=2;
            catch
                sprintf('switch to older version');
                old_field_version=1;
                field2_5;
            end                  
            component_H_x=field_H_2;                
            clear field_H_2              
        elseif i==2
            if old_field_version==1
                field2_5;
            else
                field2_6;
            end
            component_E_y=field_E_2;
            clear field_E_2
        elseif i==3
            if old_field_version==1
                field2_5;
            else
                field2_6;
            end                
            component_H_z=field_H_2;                
            clear field_H_2
        end
    end
end

%-----------------------------------------------------------
%normovani=sqrt(mu/epsilon);
if polarization==1 
    S_x=real(-(component_E_z./3.767303134483109e+02).*conj(component_H_y));
    S_z=real((component_E_x./3.767303134483109e+02).*conj(component_H_y));
    magnitude_S=sqrt(S_x.*S_x+S_z.*S_z);
elseif polarization==2
    S_x=real(376*component_E_y.*conj(component_H_z));
    S_z=real(-376*component_E_y.*conj(component_H_x));
    magnitude_S=sqrt(S_x.*S_x+S_z.*S_z);
end

%delta_x=Lambda/resolution_x;
%X=0:delta_x:Lambda;
%Y=real(S_z(:,25));
%Z=trapz(X,Y);

%cast_X=0:delta_x:1.4+delta_x;
%cast_Y=real(S_z(1:2988,25));
%cast_Z=trapz(cast_X,cast_Y);

%flux=cast_Z/Z;

size_field=size(S_x);
size_x=size_field(1);
size_z=size_field(2);
new_z=0:thickness/size_z:thickness-thickness/size_z;

figure;
set(gcf, 'Renderer', 'zbuffer');
surf(new_z.*1E6,x.*1E6,S_x); 
axis([0 thickness*1E6 x_min*1E6 x_max*1E6]);
title('S_x');
xlabel('z [\mu{}m]','FontSize',14);
ylabel('x [\mu{}m]','FontSize',14);
shading flat;

if in_octave==1 && old_field_version==1
    view_octave([90 90]);
else
    view([90 90]);    
end

figure;
set(gcf, 'Renderer', 'zbuffer');
surf(new_z.*1E6,x.*1E6,S_z); 
axis([0 thickness*1E6 x_min*1E6 x_max*1E6]);
title('S_z');
xlabel('z [\mu{}m]','FontSize',14);
ylabel('x [\mu{}m]','FontSize',14);
shading flat;

if in_octave==1 && old_field_version==1
    view_octave([90 90]);
else
    view([90 90]);    
end

figure;
set(gcf, 'Renderer', 'zbuffer');
surf(new_z.*1E6,x.*1E6,magnitude_S); 
axis([0 thickness*1E6 x_min*1E6 x_max*1E6]);
%quiver(new_z.*1E6,x.*1E6,px,py)
title('magnitude S');
xlabel('z [\mu{}m]','FontSize',14);
ylabel('x [\mu{}m]','FontSize',14);
shading flat;

if in_octave==1 && old_field_version==1
    view_octave([90 90]);
else
    view([90 90]);    
end

figure;
set(gcf, 'Renderer', 'zbuffer');
contour(new_z.*1E6,x.*1E6,magnitude_S); 
axis([0 thickness*1E6 x_min*1E6 x_max*1E6]);
hold on
[px,py] = gradient(magnitude_S,z_rozliseni*1E6,x_rozliseni*1E6);
quiver(new_z(1:10:end).*1E6,x(1:10:end).*1E6,px(1:10:end,1:10:end),py(1:10:end,1:10:end))
%quiver(new_z.*1E6,x.*1E6,px,py)
title('direction S');
xlabel('z [\mu{}m]','FontSize',14);
ylabel('x [\mu{}m]','FontSize',14);
shading flat;

if in_octave==1 && old_field_version==1
    view_octave([90 90]);
else
    view([90 90]);    
end