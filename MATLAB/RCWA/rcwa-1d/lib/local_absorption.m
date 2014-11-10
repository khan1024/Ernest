% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

epsilon_matrix_local=[];
epsilon_matrix_local_=[];

plot_field_region_1=2;
plot_field_region_2=1;
plot_field_region_3=2;

if polarization==0
    for i=4:1:6
        conical_field_component=i;
        if i==4
            try
                field2_6_local_ab;
                old_field_version=2;
            catch
                sprintf('switch to older version');
                old_field_version=1;
                field2_5_local_ab;
            end                
            component_E_x=field_E_2;                
            clear field_E_2
        elseif i==5
            try
                field2_6_local_ab;
                old_field_version=2;
            catch
                sprintf('switch to older version');
                old_field_version=1;
                field2_5_local_ab;
            end                
            component_E_y=field_E_2;                
            clear field_E_2
        elseif i==6    
            if old_field_version==1
                field2_5_local_ab;
            else
                field2_6_local_ab;
            end
            component_E_z=field_E_2;    
            clear field_E_2
        end
    end
elseif polarization==1                    
    for i=1:2:3
        planar_field_component=i;
        if i==1
            try
                field2_6_local_ab;
                old_field_version=2;
            catch
                sprintf('switch to older version');
                old_field_version=1;
                field2_5_local_ab;
            end                
            component_E_x=field_E_2;                
            clear field_E_2
        elseif i==3    
            if old_field_version==1
                field2_5_local_ab;
            else
                field2_6_local_ab;
            end
            component_E_z=field_E_2;    
            clear field_E_2
        end
    end
elseif polarization==2
    for i=2:1:2
        planar_field_component=i;
        try
            field2_6_local_ab;
            old_field_version=2;
        catch
            sprintf('switch to older version');
            old_field_version=1;
            field2_5_local_ab;
        end
        component_E_y=field_E_2;
        clear field_E_2              
    end
end

%-----------------------------------------------------------
%normovani=sqrt(mu/epsilon);
if polarization==0
    intensity_E=(abs(component_E_x)).^2+(abs(component_E_y)).^2+(abs(component_E_z)).^2;
elseif polarization==1 
    intensity_E=((abs(component_E_x)).^2+(abs(component_E_z)).^2);
elseif polarization==2
    intensity_E=(abs(component_E_y)).^2;
end

if get_local_absoption==1
    integrated_absorption;
end

if get_electric_intensity==1
    electric_intensity;
end