% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

size_field=size(total_field);
size_z=size_field(2);    
diff_z=length(total_axis_z)-size_z;    
if diff_z==0    
else    
    total_axis_new_z=total_axis_z(1:end-diff_z);        
    clear total_axis_z        
    total_axis_z=total_axis_new_z;        
    warning('z -> change scale')    
end