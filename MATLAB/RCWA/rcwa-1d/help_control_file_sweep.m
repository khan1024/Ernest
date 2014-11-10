% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% you are in octave or matlab? 
matlab_octave_version;

% this is not useful
vykresleni_urciteho_radu=2; %1-yes, 2-no
rad=0;

% check if studying_order is available 
order_index=(number_of_orders+1)/2+rad;
if order_index<=number_of_orders
    %ok
else
    error('increase number_of_orders')
end

addpath('lib');
addpath('model');

% plot field -> measurement!=0 -> measurement=0
set_calculation_field=2;
if plot_field_region_1==1 || plot_field_region_2==1 || plot_field_region_3==1
    set_calculation_field=1;
    if measurement==0
        %ok
    else
        measurement=0;
        %disp('switch to measurement=0');
    end
end

% x,y,z base matrix instabilities
if polarization==0 && change_matrix_base==2
    old_change_matrix_base=2;
elseif polarization==0 && change_matrix_base==1
    old_change_matrix_base=1;
end

switch measurement
    case 0
        if use_dispersion==1
            setup_dispersion;
        end
        grating_setup;
        
        rcwa_procedure;

        %disp('reflected orders')
        diff_R=[n(imag(k_1_z)==0),D_R(imag(k_1_z)==0)];
        
        try
            diff_R_angles=[n(imag(k_1_z)==0), D_R_angles];
        catch
        end
        
        if isempty(diff_R)==1
            %disp('only evanescent orders')
        else
            %disp('   order               efficiency')
            %disp(diff_R); %THIS IS WHAT YOU WANT
        end
        
        %disp('transmitted orders')
        diff_T=[n(imag(k_3_z)==0),D_T(imag(k_3_z)==0)];
        
        try
            diff_T_angles=[n(imag(k_3_z)==0), D_T_angles];
        catch
        end
        
        if isempty(diff_T)==1
            %disp('only evanescent orders')
        else
            %disp('   order               efficiency')
            %disp(diff_T); %THIS IS WHAT YOU WANT
        end
        
        %diff_eff=[n D_R D_T]
        sum_eff=[sum(D_R) sum(D_T) sum(D_R)+sum(D_T)];
%------------------------------------------------------------------        
    case 1
        fid = fopen('data1.dat', 'wt');
        for loop=minimum_number_of_orders:step:maximum_number_of_orders
            number_of_orders=loop;
            if use_dispersion==1
                setup_dispersion;
            end
            grating_setup;
            
            rcwa_procedure;
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record=[number_of_orders diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data1.dat;
        a=data1;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
        xlabel('number of orders','FontSize',18);
        ylabel(['efficiency ' symbol_R_or_T num2str(studying_order)],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the number of orders'] ,'FontSize',20);
        %delete('data1.dat')
%------------------------------------------------------------------
    case 2 % set number_of_orders
        data_record_2=[measurement -(number_of_orders-1)/2:1:(number_of_orders-1)/2 -(number_of_orders-1)/2:1:(number_of_orders-1)/2];
        fid = fopen('data2.dat', 'wt');
        for loop=minimum_wavelength:step:maximum_wavelength
            lambda=loop;
            
            if use_dispersion==1
                setup_dispersion;
            end
            grating_setup;
            
            rcwa_procedure;
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record_c=[loop D_R' D_T'];
            data_record_2=[data_record_2;data_record_c];
            
            data_record=[loop diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data2.dat;
        a=data2;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
        xlabel('\lambda [\mum]','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the wavelength'],'FontSize',20);
        delete('data2.dat')
%------------------------------------------------------------------    
    case 3 % 
        data_record_2=[measurement -(number_of_orders-1)/2:1:(number_of_orders-1)/2 -(number_of_orders-1)/2:1:(number_of_orders-1)/2];
        fid = fopen('data3.dat', 'wt');
        for loop=minimum_thickness:step:maximum_thickness
            if use_dispersion==1
                setup_dispersion;
            end
            grating_setup;            
            switch grating
                case 0
                    layer_thickness=loop*1E-6;
                case 1
                    thickness_total=loop;
                    blazed_grating_1;
                case 2
                    thickness_total=loop;
                    blazed_grating_2;
                case 3
                    thickness_total=loop;
                    sinusoidal_grating_y;
                case 4
                    thickness_total=loop;
                    triangular_grating;
                case 5
                    thickness_total=loop;
                    parabolic_grating;
                case 6
                    thickness_total=loop;
                    general_triangular_grating;
                case 7
                    thickness_total=loop;
                    trapezoidal_grating;
                case 8
                    thickness_total=loop;
                    tetragonal_grating;
                case 9
                    thickness_total=loop/1000;
                    gaussian_grating;
                case 10
                    thickness_total=loop;
                    thickness=thickness_total*1E-6;
                case 11
                    thickness_total=loop;
                    Kat_grating_1;

                otherwise
                    error('undefined switch')
            end
            
            rcwa_procedure;
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record_c=[loop D_R' D_T'];
            data_record_2=[data_record_2;data_record_c];
            data_record=[loop diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data3.dat;
        a=data3;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
        xlabel('thickness [\mum]','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the thickness'],'FontSize',20);
        delete('data3.dat')
%------------------------------------------------------------------    
    case 4
        data_record_2=[measurement -(number_of_orders-1)/2:1:(number_of_orders-1)/2 -(number_of_orders-1)/2:1:(number_of_orders-1)/2];
        fid = fopen('data4.dat', 'wt');
        for loop=minimum_index:step:maximum_index
            
            grating_setup;
            
            switch change_of_index
                case 1
                    n1=loop;
                    ng=loop;
                case 2
                    n1=loop;
                case 3
                    ng=loop;
            end
            
            rcwa_procedure;
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record_c=[loop D_R' D_T'];
            data_record_2=[data_record_2;data_record_c];
            data_record=[loop diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
            minimum_index=minimum_index+step;
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data4.dat;
        a=data4;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
        switch zmena_indexu
            case 1
                xlabel('refractive index n1','FontSize',18);
                %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the refractive index n1'],'FontSize',20);
            case 2
                xlabel('refractive index','FontSize',18);
                %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the refractive index n1'],'FontSize',20);
            case 3
                xlabel('refractive index','FontSize',18);
                %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the refractive index ng'],'FontSize',20);
        end
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        delete('data4.dat')
%------------------------------------------------------------------    
    case 5
        data_record_2=[measurement -(number_of_orders-1)/2:1:(number_of_orders-1)/2 -(number_of_orders-1)/2:1:(number_of_orders-1)/2];
        fid = fopen('data5.dat', 'wt');
        for loop=minimum_index:step:maximum_index
            
            grating_setup;
            
            switch change_of_index
                case 1
                    n3=loop;
                    nr=loop;
                case 2
                    n3=loop;
                case 3
                    nr=loop;
            end
            
            rcwa_procedure;
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record_c=[loop D_R' D_T'];
            data_record_2=[data_record_2;data_record_c];
            data_record=[loop diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data5.dat;
        a=data5;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
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
        delete('data5.dat')
%------------------------------------------------------------------    
    case 6
        data_record_2=[measurement -(number_of_orders-1)/2:1:(number_of_orders-1)/2 -(number_of_orders-1)/2:1:(number_of_orders-1)/2];
        fid = fopen('data6.dat', 'wt');
        for loop=minimum_angle:step:maximum_angle
            theta0=loop;
            if use_dispersion==1
                setup_dispersion;
            end
            grating_setup;
            
            rcwa_procedure;
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record_c=[loop D_R' D_T'];
            data_record_2=[data_record_2;data_record_c];
            data_record=[loop diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data6.dat;
        a=data6;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
        xlabel('incident angle \theta','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the incident angle \theta'],'FontSize',20);
        delete('data6.dat')
%------------------------------------------------------------------    
    case 7
        data_record_2=[measurement -(number_of_orders-1)/2:1:(number_of_orders-1)/2 -(number_of_orders-1)/2:1:(number_of_orders-1)/2];
        fid = fopen('data7.dat', 'wt');
        for loop=minimum_layer:step:maximum_layer
            if use_dispersion==1
                setup_dispersion;
            end
            grating_setup;
            
            number_of_layers=loop;
            switch grating
                case 1
                    blazed_grating_1;
                case 2
                    blazed_grating_1;
                case 3
                    sinusoidal_grating_y;
                case 4
                    triangular_grating;
                case 5
                    parabolic_grating;
                case 6
                    general_triangular_grating;
                case 7
                    trapezoidal_grating;
                case 8
                    tetragonal_grating;
                case 9
                    gaussian_grating;
                case 11
                    Kat_grating_1;
                otherwise
                    error('undefined switch')
            end
            
            rcwa_procedure;          
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record_c=[loop D_R' D_T'];
            data_record_2=[data_record_2;data_record_c];
            data_record=[loop diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data7.dat;
        a=data7;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
        xlabel('number of layers','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the number of layers'],'FontSize',20);
        delete('data7.dat')
%------------------------------------------------------------------    
    case 8
        data_record_2=[measurement -(number_of_orders-1)/2:1:(number_of_orders-1)/2 -(number_of_orders-1)/2:1:(number_of_orders-1)/2];
        fid = fopen('data8.dat', 'wt');
        for loop=minimum_Lambda:step:maximum_Lambda
            if use_dispersion==1
                setup_dispersion;
            end
            grating_setup;
            
            Lambda=loop;
            perioda=Lambda*1E-6;
                   
            switch grating                     
                       
                case 0                           
                    %thickness_total=Lambda*0.4;                           
                    %theta0=180/pi*asin(lambda/(2000*Lambda));rcwa;                       
                case 1                           
                    blazed_grating_1;                       
                case 2                           
                    blazed_grating_2;                       
                case 3                           
                    sinusoidal_grating_y;                       
                case 4                           
                    triangular_grating;                       
                case 5                           
                    parabolic_grating;                       
                case 6                           
                    general_triangular_grating;                       
                case 7                           
                    trapezoidal_grating;                       
                case 8                           
                    tetragonal_grating;                       
                case 9                           
                    gaussian_grating;                       
                case 10         
                case 11
                    Kat_grating_1;

                       
                otherwise
                    error('undefined switch')                   
            end
            
            rcwa_procedure;
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record_c=[loop D_R' D_T'];
            data_record_2=[data_record_2;data_record_c];
            data_record=[loop diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data8.dat;
        a=data8;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
        xlabel('grating period \Lambda [\mum]','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the grating period \Lambda'],'FontSize',20);
        delete('data8.dat')
%----------------------------------------------------------------------    
    case 9
        data_record_2=[measurement -(number_of_orders-1)/2:1:(number_of_orders-1)/2 -(number_of_orders-1)/2:1:(number_of_orders-1)/2];
        fid = fopen('data9.dat', 'wt');
        for loop=minimum_duty_cycle:step:maximum_duty_cycle
            if use_dispersion==1
                setup_dispersion;
            end
            grating_setup;
            
            duty_cycle=loop;
            
            rcwa_procedure;
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record_c=[loop D_R' D_T'];
            data_record_2=[data_record_2;data_record_c];
            data_record=[loop diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data9.dat;
        a=data9;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
        xlabel('duty cycle','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the binary duty cycle'],'FontSize',20);
        delete('data9.dat')
    %------------------------------------------------------------------    
    case 10
        data_record_2=[measurement -(number_of_orders-1)/2:1:(number_of_orders-1)/2 -(number_of_orders-1)/2:1:(number_of_orders-1)/2];
        fid = fopen('data10.dat', 'wt');
        for loop=minimum_angle_psi:step:maximum_angle_psi
            if use_dispersion==1
                setup_dispersion;
            end
            grating_setup;
            
            psi0=loop;
            
            rcwa_procedure;
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record_c=[loop D_R' D_T'];
            data_record_2=[data_record_2;data_record_c];
            data_record=[loop diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data10.dat;
        a=data10;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
        xlabel('angle \psi','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the angle \psi'],'FontSize',20);
        delete('data10.dat')
    %------------------------------------------------------------------    
    case 11
        data_record_2=[measurement -(number_of_orders-1)/2:1:(number_of_orders-1)/2 -(number_of_orders-1)/2:1:(number_of_orders-1)/2];
        fid = fopen('data11.dat', 'wt');
        for loop=minimum_angle_phi:step:maximum_angle_phi
            if use_dispersion==1
                setup_dispersion;
            end
            grating_setup;
            phi0=loop;
            
            rcwa_procedure;
            
            if diffraction_efficiencies_c==1
                diffraction_efficiencies=D_R;
                symbol_R_or_T='R';
            else
                diffraction_efficiencies=D_T;                
                symbol_R_or_T='T';
            end
            data_record_c=[loop D_R' D_T'];
            data_record_2=[data_record_2;data_record_c];
            data_record=[loop diffraction_efficiencies(p+studying_order)];
            fprintf(fid, '%0f %0f\n', data_record);
        end
        fclose(fid);
        save_all_measurements_results;
        figure;
        load data11.dat;
        a=data11;
        axes('FontSize',14);
        plot(a(:,1),a(:,2),'k','Marker','+','MarkerSize',3);
        xlabel('angle \phi','FontSize',18);
        ylabel(['efficiency ' num2str(studying_order) symbol_R_or_T ],'FontSize',18);
        %title(['dependence of the diffraction efficiency ' num2str(studying_order) symbol_R_or_T  ' on the angle \phi'],'FontSize',20);
        delete('data11.dat')
end