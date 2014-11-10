% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%faktorization=1;
%save=1;
%save_as='muj.txt'
%measurement=9;
if save_m==1
switch measurement
    case 1
        copyfile('data1.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_number_of_orders=%0f ',minimum_number_of_orders);
        fprintf(fid, '\n maximum_number_of_orders=%0f ',maximum_number_of_orders);
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
        fprintf(fid, '\n psi0=%0f ',psi0);
        fprintf(fid, '\n phi0=%0f ',phi0);
        fclose(fid);
    case 2
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data2.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_wavelength=%0f ',minimum_wavelength);
        fprintf(fid, '\n maximum_wavelength=%0f ',maximum_wavelength);
        fprintf(fid, '\n step=%0f ',step);
        fprintf(fid, '\n number_of_orders=%0f ',number_of_orders);
        fprintf(fid, '\n theta0=%0f ',theta0);
        fprintf(fid, '\n n1=%0f-%0f*j ',real(n1),abs(imag(n1)));
        fprintf(fid, '\n n3=%0f-%0f*j ',real(n3),abs(imag(n3)));
        fprintf(fid, '\n ng=%0f-%0f*j ',real(ng),abs(imag(ng)));
        fprintf(fid, '\n nr=%0f-%0f*j ',real(nr),abs(imag(nr)));
        fprintf(fid, '\n polarization=%0f ',polarization);
        fprintf(fid, '\n faktorization=%0f ',faktorization);
        fprintf(fid, '\n Lambda=%0f ',Lambda);
        fprintf(fid, '\n thickness_total=%0f ',thickness_total);
        fprintf(fid, '\n grating=%0f ',grating);
        fprintf(fid, '\n number_of_layers=%0f ',number_of_layers);
        fprintf(fid, '\n psi0=%0f ',psi0);
        fprintf(fid, '\n phi0=%0f ',phi0);
        fclose(fid);
    case 3
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data3.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_thickness=%0f ',minimum_thickness);
        fprintf(fid, '\n maximum_thickness=%0f ',maximum_thickness);
        fprintf(fid, '\n step=%0f ',step);
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
        fprintf(fid, '\n grating=%0f ',grating);
        fprintf(fid, '\n number_of_layers=%0f ',number_of_layers);
        fprintf(fid, '\n psi0=%0f ',psi0);
        fprintf(fid, '\n phi0=%0f ',phi0);
        fclose(fid);
    case 4
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data4.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_index=%0f ',minimum_index);
        fprintf(fid, '\n maximum_index=%0f ',maximum_index);
        fprintf(fid, '\n step=%0f ',step);
        fprintf(fid, '\n zmena_indexu=%0f ',zmena_indexu);
        fprintf(fid, '\n number_of_orders=%0f ',number_of_orders);
        fprintf(fid, '\n theta0=%0f ',theta0);
        fprintf(fid, '\n n3=%0f-%0f*j ',real(n3),abs(imag(n3)));
        fprintf(fid, '\n lambda=%0f ',lambda);
        fprintf(fid, '\n nr=%0f-%0f*j ',real(nr),abs(imag(nr)));
        fprintf(fid, '\n polarization=%0f ',polarization);
        fprintf(fid, '\n faktorization=%0f ',faktorization);
        fprintf(fid, '\n Lambda=%0f ',Lambda);
        fprintf(fid, '\n thickness_total=%0f ',thickness_total);
        fprintf(fid, '\n grating=%0f ',grating);
        fprintf(fid, '\n number_of_layers=%0f ',number_of_layers);
        fprintf(fid, '\n psi0=%0f ',psi0);
        fprintf(fid, '\n phi0=%0f ',phi0);
        fclose(fid);
        save(savefile,'-append', 'zmena_indexu')
    case 5
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data5.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_index=%0f ',minimum_index);
        fprintf(fid, '\n maximum_index=%0f ',maximum_index);
        fprintf(fid, '\n step=%0f ',step);
        fprintf(fid, '\n zmena_indexu=%0f ',zmena_indexu);
        fprintf(fid, '\n number_of_orders=%0f ',number_of_orders);
        fprintf(fid, '\n theta0=%0f ',theta0);
        fprintf(fid, '\n n1=%0f-%0f*j ',real(n1),abs(imag(n1)));
        fprintf(fid, '\n lambda=%0f ',lambda);
        fprintf(fid, '\n ng=%0f-%0f*j ',real(ng),abs(imag(ng)));
        fprintf(fid, '\n polarization=%0f ',polarization);
        fprintf(fid, '\n faktorization=%0f ',faktorization);
        fprintf(fid, '\n Lambda=%0f ',Lambda);
        fprintf(fid, '\n thickness_total=%0f ',thickness_total);
        fprintf(fid, '\n grating=%0f ',grating);
        fprintf(fid, '\n number_of_layers=%0f ',number_of_layers);
        fprintf(fid, '\n psi0=%0f ',psi0);
        fprintf(fid, '\n phi0=%0f ',phi0);
        fclose(fid);
        save(savefile,'-append', 'zmena_indexu')
    case 6
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data6.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_angle=%0f ',minimum_angle);
        fprintf(fid, '\n maximum_angle=%0f ',maximum_angle);
        fprintf(fid, '\n step=%0f ',step);
        fprintf(fid, '\n number_of_orders=%0f ',number_of_orders);
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
        fprintf(fid, '\n psi0=%0f ',psi0);
        fprintf(fid, '\n phi0=%0f ',phi0);
        fclose(fid);
    case 7
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data7.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_layer=%0f ',minimum_layer);
        fprintf(fid, '\n maximum_layer=%0f ',maximum_layer);
        fprintf(fid, '\n step=%0f ',step);
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
        fprintf(fid, '\n psi0=%0f ',psi0);
        fprintf(fid, '\n phi0=%0f ',phi0);
        fclose(fid);
    case 8
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data8.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_Lambda=%0f ',minimum_Lambda);
        fprintf(fid, '\n maximum_Lambda=%0f ',maximum_Lambda);
        fprintf(fid, '\n step=%0f ',step);
        fprintf(fid, '\n number_of_orders=%0f ',number_of_orders);
        fprintf(fid, '\n theta0=%0f ',theta0);
        fprintf(fid, '\n n1=%0f-%0f*j ',real(n1),abs(imag(n1)));
        fprintf(fid, '\n n3=%0f-%0f*j ',real(n3),abs(imag(n3)));
        fprintf(fid, '\n lambda=%0f ',lambda);
        fprintf(fid, '\n ng=%0f-%0f*j ',real(ng),abs(imag(ng)));
        fprintf(fid, '\n nr=%0f-%0f*j ',real(nr),abs(imag(nr)));
        fprintf(fid, '\n polarization=%0f ',polarization);
        fprintf(fid, '\n faktorization=%0f ',faktorization);
        fprintf(fid, '\n thickness_total=%0f ',thickness_total);
        fprintf(fid, '\n grating=%0f ',grating);
        fprintf(fid, '\n number_of_layers=%0f ',number_of_layers);
        fprintf(fid, '\n psi0=%0f ',psi0);
        fprintf(fid, '\n phi0=%0f ',phi0);
        fclose(fid);
    case 9
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data9.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_duty_cycle=%0f ',minimum_duty_cycle);
        fprintf(fid, '\n maximum_duty_cycle=%0f ',maximum_duty_cycle);
        fprintf(fid, '\n step=%0f ',step);
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
        fprintf(fid, '\n psi0=%0f ',psi0);
        fprintf(fid, '\n phi0=%0f ',phi0);
        fclose(fid);
    case 10
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data10.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_angle_psi=%0f ',minimum_angle_psi);
        fprintf(fid, '\n maximum_angle_psi=%0f ',maximum_angle_psi);
        fprintf(fid, '\n step=%0f ',step);
        fprintf(fid, '\n number_of_orders=%0f ',number_of_orders);
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
        fprintf(fid, '\n phi0=%0f ',phi0);
        fclose(fid);
    case 11
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data10.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fprintf(fid, '\n diffraction_efficiencies_c=%0f ',diffraction_efficiencies_c);
        fprintf(fid, '\n studying_order=%0f ',studying_order);
        fprintf(fid, '\n minimum_angle_phi=%0f ',minimum_angle_phi);
        fprintf(fid, '\n maximum_angle_phi=%0f ',maximum_angle_phi);
        fprintf(fid, '\n step=%0f ',step);
        fprintf(fid, '\n number_of_orders=%0f ',number_of_orders);
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
        fprintf(fid, '\n psi0=%0f ',psi0);
        fclose(fid);
    otherwise
        savefile = [save_as '.mat'];
        save(savefile, 'data_record_2')
        copyfile('data_new.dat',[save_as '.dat']);
        fid = fopen([save_as '.sav'], 'wt');
        fprintf(fid, ' measurement=%0f ',measurement);
        fclose(fid);
        
end
if measurement==1
    movefile([save_as '.sav'],'measurements')
    movefile([save_as '.dat'],'measurements')
else
    movefile([save_as '.mat'],'measurements')
    movefile([save_as '.sav'],'measurements')
    movefile([save_as '.dat'],'measurements')
end

end