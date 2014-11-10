% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

if new_figure==1
    figure;
end

if slice==1 || slice2==1
    display('slice is on');
    if plot_field_region_1==1
        display('plot_field_region_1=2');
    end
    
    if plot_field_region_3==1
        display('plot_field_region_1=2');
    end
    
    
    plot_field_region_1=2;
    plot_field_region_2=1;
    plot_field_region_3=2;
end

set(gcf, 'Renderer', 'zbuffer');
axes('FontSize',17)
if plot_field_region_1==1 && plot_field_region_2==1 && plot_field_region_3==1
    total_axis_z=[z_11 z_22 z_33];
    
    switch planar_field_component
        case 1
            if polarization==1
                total_field=view_field([field_E_1 field_E_2 field_E_3],field_component_operation);
            elseif polarization==2
                total_field=view_field([field_H_1 field_H_2 field_H_3],field_component_operation);
            end
        case 2
            if polarization==1
                total_field=view_field([field_H_1 field_H_2 field_H_3],field_component_operation);
            elseif polarization==2
                total_field=view_field([field_E_1 field_E_2 field_E_3],field_component_operation);
            end
        case 3
            if polarization==1
                total_field=view_field([field_E_1 field_E_2 field_E_3],field_component_operation);
            elseif polarization==2
                total_field=view_field([field_H_1 field_H_2 field_H_3],field_component_operation);
            end
    end
    change_scale_z;    
    surf(total_axis_z.*1E6,x.*1E6,total_field);
    axis([-size_region_1*thickness_total (1+size_region_3)*thickness_total x_min*1e6 x_max*1e6])
elseif plot_field_region_1==2 && plot_field_region_2==1 && plot_field_region_3==1
    total_axis_z=[z_22 z_33];
    switch planar_field_component
        case 1
            if polarization==1
                total_field=view_field([field_E_2 field_E_3],field_component_operation);
            elseif polarization==2
                total_field=view_field([field_H_2 field_H_3],field_component_operation);
            end
        case 2
            if polarization==1
                total_field=view_field([field_H_2 field_H_3],field_component_operation);
            elseif polarization==2
                total_field=view_field([field_E_2 field_E_3],field_component_operation);
            end
        case 3
            if polarization==1
                total_field=view_field([field_E_2 field_E_3],field_component_operation);
            elseif polarization==2
                total_field=view_field([field_H_2 field_H_3],field_component_operation);
            end
    end
    change_scale_z;
    surf(total_axis_z.*1E6,x.*1E6,total_field);
    axis([0 (1+size_region_3)*thickness_total x_min*1e6 x_max*1e6])
elseif plot_field_region_1==2 && plot_field_region_2==2 && plot_field_region_3==1
    set_view_grating_boundaries=2;
    total_axis_z=z_33;
    switch planar_field_component
        case 1
            if polarization==1
                total_field=view_field(field_E_3,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_H_3,field_component_operation);
            end
        case 2
            if polarization==1
                total_field=view_field(field_H_3,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_E_3,field_component_operation);
            end
        case 3
            if polarization==1
                total_field=view_field(field_E_3,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_H_3,field_component_operation);
            end
    end
    change_scale_z;
    surf(total_axis_z.*1E6,x.*1E6,total_field);
    axis([thickness_total (1+size_region_3)*thickness_total x_min*1e6 x_max*1e6])
elseif plot_field_region_1==1 && plot_field_region_2==2 && plot_field_region_3==2
    set_view_grating_boundaries=2;
    total_axis_z=z_11;
    switch planar_field_component
        case 1
            if polarization==1
                total_field=view_field(field_E_1,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_H_1,field_component_operation);
            end
        case 2
            if polarization==1
                total_field=view_field(field_H_1,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_E_1,field_component_operation);
            end
        case 3
            if polarization==1
                total_field=view_field(field_E_1,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_H_1,field_component_operation);
            end
    end
    change_scale_z;
    surf(total_axis_z.*1E6,x.*1E6,total_field);
    axis([-size_region_1*thickness_total 0 x_min*1e6 x_max*1e6])
elseif plot_field_region_1==2 && plot_field_region_2==1 && plot_field_region_3==2
    total_axis_z=z_22;
        switch planar_field_component
            case 1            
                if polarization==1                
                    total_field=view_field(field_E_2,field_component_operation);            
                elseif polarization==2                
                    total_field=view_field(field_H_2,field_component_operation);            
                end
            case 2            
                if polarization==1                
                    total_field=view_field(field_H_2,field_component_operation);           
                elseif polarization==2                
                    total_field=view_field(field_E_2,field_component_operation);           
                end
            case 3            
                if polarization==1                
                    total_field=view_field(field_E_2,field_component_operation);            
                elseif polarization==2                
                    total_field=view_field(field_H_2,field_component_operation);           
                end
        end

        change_scale_z;
        surf(total_axis_z.*1E6,x.*1E6,total_field);
        axis([0 thickness_total x_min*1e6 x_max*1e6])    
elseif plot_field_region_1==1 && plot_field_region_2==1 && plot_field_region_3==2
    total_axis_z=[z_11 z_22];
    switch planar_field_component
        case 1
            if polarization==1
                total_field=view_field([field_E_1 field_E_2],field_component_operation);
            elseif polarization==2
                total_field=view_field([field_H_1 field_H_2],field_component_operation);
            end
        case 2
            if polarization==1
                total_field=view_field([field_H_1 field_H_2],field_component_operation);
            elseif polarization==2
                total_field=view_field([field_E_1 field_E_2],field_component_operation);
            end
        case 3
            if polarization==1
                total_field=view_field([field_E_1 field_E_2],field_component_operation);
            elseif polarization==2
                total_field=view_field([field_H_1 field_H_2],field_component_operation);
            end
    end
    change_scale_z;
    surf(total_axis_z.*1E6,x.*1E6,total_field);
    axis([-size_region_1*thickness_total thickness_total x_min*1e6 x_max*1e6])
elseif plot_field_region_1==1 && plot_field_region_2==2 && plot_field_region_3==1
    set_view_grating_boundaries=2;
    total_axis_z=z_11;
    switch planar_field_component
        case 1
            if polarization==1
                total_field=view_field(field_E_1,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_H_1,field_component_operation);
            end
        case 2
            if polarization==1
                total_field=view_field(field_H_1,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_E_1,field_component_operation);
            end
        case 3
            if polarization==1
                total_field=view_field(field_E_1,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_H_1,field_component_operation);
            end
    end
    change_scale_z;
    surf(total_axis_z.*1E6,x.*1E6,total_field);
    axis([-size_region_1*thickness_total 0 x_min*1e6 x_max*1e6])
    
    xlabel('z [\mu{}m]','FontSize',17);
    ylabel('x [\mu{}m]','FontSize',17);
    
    try     
        colorbar('FontSize',17);
    catch        
        colorbar;
    end
    
    shading flat;
    
    if in_octave==1 && old_field_version==1
        view_octave([90 90]);
    else
        view([90 90]);
    end
    
    figure
    set(gcf, 'Renderer', 'zbuffer');
    axes('FontSize',17)
    
    total_axis_z=z_33;
    switch planar_field_component
        case 1
            if polarization==1
                total_field=view_field(field_E_3,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_H_3,field_component_operation);
            end
        case 2
            if polarization==1
                total_field=view_field(field_H_3,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_E_3,field_component_operation);
            end
        case 3
            if polarization==1
                total_field=view_field(field_E_3,field_component_operation);
            elseif polarization==2
                total_field=view_field(field_H_3,field_component_operation);
            end
    end
    change_scale_z;
    surf(total_axis_z.*1E6,x.*1E6,total_field);
    axis([thickness_total (1+size_region_3)*thickness_total x_min*1e6 x_max*1e6])
    
end

xlabel('z [\mu{}m]','FontSize',17);  
ylabel('x [\mu{}m]','FontSize',17);

try 
    colorbar('FontSize',17);
catch
    colorbar;
end

shading flat;

if set_view_grating_boundaries==1
    hold on;
    framework;
end

if in_octave==1 && old_field_version==1
    view_octave([90 90]);
else
    view([90 90]);    
end

if plot_field_region_1==2 && plot_field_region_2==1 && plot_field_region_3==2       
    if slice==1           
        if new_figure==1                
            figure;            
        end
        if slice_position==0                
            plot(x.*1E6,total_field(:,1));            
        else
            plot(x.*1E6,total_field(:,resolution_z*slice_position));            
        end
        xlabel('x [\mu{}m]','FontSize',14);        
    end
    if slice2==1
        if new_figure==1        
            figure;    
        end
        if slice2_position==0
            change_scale_z;
            plot(total_axis_z.*1E6,total_field(1,:));            
        else
            change_scale_z;
            plot(total_axis_z.*1E6,total_field(resolution_z*slice2_position,:));            
        end
        xlabel('z [\mu{}m]','FontSize',14);
   end
end

if save_field==1
    save_field_function;
end