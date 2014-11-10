% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

if save_field==1
    savefile=[save_field_as '.mat'];
    if plot_field_region_1==1 && plot_field_region_2==1 && plot_field_region_3==1
        
        if polarization==0
            if conical_field_component==1 || conical_field_component==2 || conical_field_component==3 % H
                save(savefile,'z_11','z_22','z_33')            
                save(savefile,'-append', 'field_H_1','field_H_2','field_H_3')                
            else
                save(savefile,'z_11','z_22','z_33')            
                save(savefile,'-append', 'field_E_1','field_E_2','field_E_3')
            end
        elseif polarization==1
            if planar_field_component==1 || planar_field_component==3 % E
                save(savefile,'z_11','z_22','z_33')            
                save(savefile,'-append', 'field_E_1','field_E_2','field_E_3')                
            else
                save(savefile,'z_11','z_22','z_33')            
                save(savefile,'-append', 'field_H_1','field_H_2','field_H_3')                
            end
        elseif polarization==2
            if planar_field_component==1 || planar_field_component==3 % H
                save(savefile,'z_11','z_22','z_33')            
                save(savefile,'-append', 'field_H_1','field_H_2','field_H_3')                
            else
                save(savefile,'z_11','z_22','z_33')            
                save(savefile,'-append', 'field_E_1','field_E_2','field_E_3')                
            end
        end
        
    elseif plot_field_region_1==2 && plot_field_region_2==1 && plot_field_region_3==1
        if polarization==0
            if conical_field_component==1 || conical_field_component==2 || conical_field_component==3 % H
                save(savefile,'z_22','z_33')            
                save(savefile,'-append','field_H_2','field_H_3')                
            else
                save(savefile,'z_22','z_33')            
                save(savefile,'-append','field_E_2','field_E_3')
            end
        elseif polarization==1
            if planar_field_component==1 || planar_field_component==3 % E
                save(savefile,'z_22','z_33')            
                save(savefile,'-append','field_E_2','field_E_3')                
            else
                save(savefile,'z_22','z_33')            
                save(savefile,'-append','field_H_2','field_H_3')                
            end
        elseif polarization==2
            if planar_field_component==1 || planar_field_component==3 % H
                save(savefile,'z_22','z_33')            
                save(savefile,'-append','field_H_2','field_H_3')                
            else
                save(savefile,'z_22','z_33')            
                save(savefile,'-append','field_E_2','field_E_3')                
            end
        end
    elseif plot_field_region_1==2 && plot_field_region_2==2 && plot_field_region_3==1
        if polarization==0
            if conical_field_component==1 || conical_field_component==2 || conical_field_component==3 % H
                save(savefile,'z_33')            
                save(savefile,'-append', 'field_H_3')                
            else
                save(savefile,'z_33')            
                save(savefile,'-append', 'field_E_3')
            end
        elseif polarization==1
            if planar_field_component==1 || planar_field_component==3 % E
                save(savefile,'z_33')            
                save(savefile,'-append', 'field_E_3')                
            else
                save(savefile,'z_33')            
                save(savefile,'-append','field_H_3')                
            end
        elseif polarization==2
            if planar_field_component==1 || planar_field_component==3 % H
                save(savefile,'z_33')            
                save(savefile,'-append', 'field_H_3')                
            else
                save(savefile,'z_33')            
                save(savefile,'-append','field_E_3')                
            end
        end
    elseif plot_field_region_1==1 && plot_field_region_2==2 && plot_field_region_3==2
        if polarization==0
            if conical_field_component==1 || conical_field_component==2 || conical_field_component==3 % H
                save(savefile,'z_11')            
                save(savefile,'-append', 'field_H_1')                
            else
                save(savefile,'z_11')            
                save(savefile,'-append', 'field_E_1')
            end
        elseif polarization==1
            if planar_field_component==1 || planar_field_component==3 % E
                save(savefile,'z_11')            
                save(savefile,'-append', 'field_E_1')                
            else
                save(savefile,'z_11')            
                save(savefile,'-append', 'field_H_1')                
            end
        elseif polarization==2
            if planar_field_component==1 || planar_field_component==3 % H
                save(savefile,'z_11')            
                save(savefile,'-append', 'field_H_1')                
            else
                save(savefile,'z_11')            
                save(savefile,'-append', 'field_E_1')                
            end
        end
    elseif plot_field_region_1==2 && plot_field_region_2==1 && plot_field_region_3==2
        if polarization==0
            if conical_field_component==1 || conical_field_component==2 || conical_field_component==3 % H
                save(savefile,'z_22')            
                save(savefile,'-append', 'field_H_2')                
            else
                save(savefile,'z_22')            
                save(savefile,'-append','field_E_2')
            end
        elseif polarization==1
            if planar_field_component==1 || planar_field_component==3 % E
                save(savefile,'z_22')            
                save(savefile,'-append','field_E_2')                
            else
                save(savefile,'z_22')            
                save(savefile,'-append','field_H_2')                
            end
        elseif polarization==2
            if planar_field_component==1 || planar_field_component==3 % H
                save(savefile,'z_22')            
                save(savefile,'-append','field_H_2')                
            else
                save(savefile,'z_22')            
                save(savefile,'-append','field_E_2')                
            end
        end
    elseif plot_field_region_1==1 && plot_field_region_2==1 && plot_field_region_3==2
        if polarization==0
            if conical_field_component==1 || conical_field_component==2 || conical_field_component==3 % H
                save(savefile,'z_11','z_22')            
                save(savefile,'-append', 'field_H_1','field_H_2')                
            else
                save(savefile,'z_11','z_22')            
                save(savefile,'-append', 'field_E_1','field_E_2')
            end
        elseif polarization==1
            if planar_field_component==1 || planar_field_component==3 % E
                save(savefile,'z_11','z_22')            
                save(savefile,'-append', 'field_E_1','field_E_2')                
            else
                save(savefile,'z_11','z_22')            
                save(savefile,'-append', 'field_H_1','field_H_2')                
            end
        elseif polarization==2
            if planar_field_component==1 || planar_field_component==3 % H
                save(savefile,'z_11','z_22')            
                save(savefile,'-append', 'field_H_1','field_H_2')                
            else
                save(savefile,'z_11','z_22')            
                save(savefile,'-append', 'field_E_1','field_E_2')                
            end
        end
    end
    
    if grating==10 
        framework_grating_10;
    else
        framework_grating_0_9;
    end 
    
    save(savefile,'-append', 'x')
    save(savefile,'-append', 'polarization', 'Lambda', 'thickness_total', 'number_of_plotted_period', 'size_region_1', 'size_region_3', 'resolution_x', 'resolution_z', 'planar_field_component', 'plot_field_region_1','plot_field_region_2','plot_field_region_3','in_octave','grating')
    
    if grating==10 
        save(savefile,'-append', 'data','n1','n3');
    else
        save(savefile,'-append', 'layer_thickness','duty_cycle','shift','ng','nr','n1','n3');
    end 
    
    fid = fopen([save_field_as '.sav'], 'wt');
    fprintf(fid, '\n planar_field_component=%0f ',planar_field_component);
    fprintf(fid, '\n plot_field_region_1=%0f ',plot_field_region_1);
    fprintf(fid, '\n plot_field_region_2=%0f ',plot_field_region_2);
    fprintf(fid, '\n plot_field_region_3=%0f ',plot_field_region_3);
    fprintf(fid, '\n number_of_plotted_period=%0f ',number_of_plotted_period);
    fprintf(fid, '\n size_region_1=%0f ',size_region_1);
    fprintf(fid, '\n size_region_3=%0f ',size_region_3);
    fprintf(fid, '\n resolution_x=%0f ',resolution_x);
    fprintf(fid, '\n resolution_z=%0f ',resolution_z);
    fprintf(fid, '\n slice=%0f ',slice);
    fprintf(fid, '\n slice_position=%0f ',slice_position);
    fprintf(fid, '\n number_of_orders=%0f ',number_of_orders);
    fprintf(fid, '\n theta0=%0f ',theta0);
    fprintf(fid, '\n n1=%0f-%0f*j ',real(n1),abs(imag(n1)));
    fprintf(fid, '\n n3=%0f-%0f*j ',real(n3),abs(imag(n3)));
    fprintf(fid, '\n lambda=%0f ',lambda);
    fprintf(fid, '\n ng=%0f-%0f*j ',real(ng),abs(imag(ng)));
    fprintf(fid, '\n nr=%0f-%0f*j ',real(nr),abs(imag(nr)));
    fprintf(fid, '\n polarization=%0f ',polarization);
    fprintf(fid, '\n faktorization=%0f ',faktorization);
    fprintf(fid, '\n Lambda=%0f ',Lambda);
    fprintf(fid, '\n thickness_total=%0f ',thickness_total);
    fprintf(fid, '\n grating=%0f ',grating);
    fprintf(fid, '\n number_of_layers=%0f ',number_of_layers);
    fclose(fid);
    movefile([save_field_as '.mat'],'measurements')
    movefile([save_field_as '.sav'],'measurements')
end