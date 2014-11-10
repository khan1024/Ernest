% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

function em_field_procedure=view_field(em_field,procedure_field)

if procedure_field==0
    em_field_procedure=abs(em_field).^2;
elseif procedure_field==1
    em_field_procedure=abs(em_field);
elseif procedure_field==2
    em_field_procedure=angle(em_field);
elseif procedure_field==3
    em_field_procedure=real(em_field);
elseif procedure_field==4
    em_field_procedure=imag(em_field);
end