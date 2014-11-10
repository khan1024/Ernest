% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

addpath('measurements');
addpath('lib');
figure;

%----------------------------------------------
load test;
diffraction_efficiencies_c=1; %1--D_R, 2--D_T
studying_order=1;
export_studying_order=2; % 1-yes, 2-no;
save_export_as='conical_sinusoidal_grating_p1T.txt';
%----------------------------------------------

a=data_record_2;
measurement=a(1,1);

data_size=size(data_record_2);
number_of_orders=(data_size(2)-1)/2;
O_R=data_size(2)-number_of_orders-(number_of_orders-1)/2;
O_T=data_size(2)-(number_of_orders-1)/2;

if studying_order>(number_of_orders-1)/2
error('this order does not exists')
end

if diffraction_efficiencies_c==1
    position=O_R+studying_order;
else
    position=O_T+studying_order;
end

if diffraction_efficiencies_c==1
    symbol_R_or_T='R';
else
    symbol_R_or_T='T';
end

switch measurement
    case 1
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        xlabel('number of orders','FontSize',18);
        ylabel(['efficiency ' symbol_R_or_T num2str(studying_order)],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the number of orders'] ,'FontSize',20);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    case 2
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        xlabel('\lambda [\mum]','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the wavelength'],'FontSize',20);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    case 3
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        xlabel('thickness [\mum]','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the thickness'],'FontSize',20);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    case 4
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        switch zmena_indexu
            case 1
                xlabel('refractive index n1','FontSize',18);
                %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the refractive index n1'],'FontSize',20);
            case 2
                xlabel('refractive index n1','FontSize',18);
                %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the refractive index n1'],'FontSize',20);
            case 3
                xlabel('refractive index ng','FontSize',18);
                %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the refractive index ng'],'FontSize',20);
        end
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    case 5
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        switch zmena_indexu
            case 1
                xlabel('refractive index n3','FontSize',18);
                %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the refractive index n3'],'FontSize',20);
            case 2
                xlabel('refractive index n3','FontSize',18);
                %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the refractive index n3'],'FontSize',20);
            case 3
                xlabel('refractive index nr','FontSize',18);
                %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the refractive index nr'],'FontSize',20);
        end
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    case 6
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
                
        xlabel('angle \theta','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the incident angle \theta'],'FontSize',20);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    case 7
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        xlabel('number of layers','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the number of layers'],'FontSize',20);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    case 8
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        xlabel('grating period \Lambda [\mum]','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the grating period \Lambda'],'FontSize',20);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    case 9
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        xlabel('duty cycle','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the binary duty cycle'],'FontSize',20);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    case 10
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        xlabel('angle \psi','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the angle \psi'],'FontSize',20);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    case 11
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        xlabel('angle \phi','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the angle \phi'],'FontSize',20);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
    otherwise
        axes('FontSize',14);
        plot(a(2:data_size(1),1),a(2:data_size(1),position),'k','Marker','+','MarkerSize',3);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        axis([a(2,1) a(data_size(1),1) min(a(2:data_size(1),position)) max(a(2:data_size(1),position))])
end

if export_studying_order==1
    fid = fopen(save_export_as, 'wt');
    for j=2:1:data_size(1)
        fprintf(fid, '\n%0f %0.12f', a(j,1), a(j,position));
    end
    fclose(fid);
end
