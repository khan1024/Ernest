% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

switch measurement
    case 0
        if use_dispersion==1
            setup_dispersion;
        end
        initialization_input_file;
        rcwa_procedure;
        save('temp.mat','D_R_0', 'D_T_0','D_R_90','D_T_90','D_R','D_T','N_X','N_Y');
        
        position_diff_efficiency=1+(diffraction_efficiency_order(2)+N_Y)+(2*N_Y+1)*(diffraction_efficiency_order(1)+N_Y);
        record=[diffraction_efficiency_order(1) diffraction_efficiency_order(2) D_R_0(position_diff_efficiency) D_T_0(position_diff_efficiency) sum(D_R_0) sum(D_T_0) D_R_90(position_diff_efficiency) D_T_90(position_diff_efficiency) sum(D_R_90) sum(D_T_90) D_R(position_diff_efficiency) D_T(position_diff_efficiency) sum(D_R) sum(D_T)];
        fprintf(1,'order  [%5.2f,%6.2f] R_0=%12.10f T_0=%12.10f sumR=%12.10f sumT=%12.10f R_90=%12.10f T_90=%12.10f sumR=%12.10f sumT=%12.10f R_psi=%12.10f T_psi=%12.10f sumR=%12.10f sumT=%12.10f\n', record);
    case 1
        fid = fopen('data1.txt', 'wt');
        for cycle=minimum_N_X_N_Y:step:maximum_N_X_N_Y            
            N_X=cycle;
            N_Y=N_X;
            if use_dispersion==1
                setup_dispersion;
            end
            
            time2=tic;
            
            initialization_input_file;
            rcwa_procedure;
            
            position_diff_efficiency=1+(diffraction_efficiency_order(2)+N_Y)+(2*N_Y+1)*(diffraction_efficiency_order(1)+N_Y);
            record=[cycle toc(time2) D_R_0(position_diff_efficiency) D_T_0(position_diff_efficiency) sum(D_R_0) sum(D_T_0) D_R_90(position_diff_efficiency) D_T_90(position_diff_efficiency) sum(D_R_90) sum(D_T_90) D_R(position_diff_efficiency) D_T(position_diff_efficiency) sum(D_R) sum(D_T)];
            fprintf(fid, '%0f %0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record);
            fprintf('%0f %0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record);
        end
        fclose(fid);
        
        if save_m==1
            copyfile('data1.txt',save_m_as)
        end
    case 2
        record_0=[];
        record_90=[];
        record_psi=[];
        fid = fopen('data2.txt', 'wt');
        for cycle=minimum_wavelength:step:maximum_wavelength            
            lambda=cycle;
            if use_dispersion==1
                setup_dispersion;
            end
            
            initialization_input_file;
            
            rcwa_procedure;
            
            save_all_results;
            
            position_diff_efficiency=1+(diffraction_efficiency_order(2)+N_Y)+(2*N_Y+1)*(diffraction_efficiency_order(1)+N_Y);
            record=[cycle D_R_0(position_diff_efficiency) D_T_0(position_diff_efficiency) sum(D_R_0) sum(D_T_0) D_R_90(position_diff_efficiency) D_T_90(position_diff_efficiency) sum(D_R_90) sum(D_T_90) D_R(position_diff_efficiency) D_T(position_diff_efficiency) sum(D_R) sum(D_T)];            fprintf(fid, '%0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record);
            fprintf('%0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record);
        end
        fclose(fid);
        
        if save_m==1
            copyfile('data2.txt',save_m_as)
        end
    case 3
        record_0=[];
        record_90=[];
        record_psi=[];
        fid = fopen('data3.txt', 'wt');
        for cycle=minimum_thickness:step:maximum_thickness            
            thickness=cycle;            
            if use_dispersion==1
                setup_dispersion;
            end
            initialization_input_file;
            rcwa_procedure;
            
            save_all_results;
            
            position_diff_efficiency=1+(diffraction_efficiency_order(2)+N_Y)+(2*N_Y+1)*(diffraction_efficiency_order(1)+N_Y);
            record=[cycle D_R_0(position_diff_efficiency) D_T_0(position_diff_efficiency) sum(D_R_0) sum(D_T_0) D_R_90(position_diff_efficiency) D_T_90(position_diff_efficiency) sum(D_R_90) sum(D_T_90) D_R(position_diff_efficiency) D_T(position_diff_efficiency) sum(D_R) sum(D_T)];
            fprintf(fid, '%0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record);
            fprintf('%0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record);
        end
        fclose(fid);
        
        if save_m==1
            copyfile('data3.txt',save_m_as)
        end
        
        a=load(save_m_as);
        plot(a(:,1),a(:,2));
    case 4
        record_0=[];
        record_90=[];
        record_psi=[];
        fid = fopen('data4.txt', 'wt');
        for cycle=minimum_angle_theta:step:maximum_angle_theta            
            theta0=cycle;            
            if use_dispersion==1
                setup_dispersion;
            end
            rcwa_procedure;
            initialization_input_file;
            save_all_results;
            
            position_diff_efficiency=1+(diffraction_efficiency_order(2)+N_Y)+(2*N_Y+1)*(diffraction_efficiency_order(1)+N_Y);
            record=[cycle D_R_0(position_diff_efficiency) D_T_0(position_diff_efficiency) sum(D_R_0) sum(D_T_0) D_R_90(position_diff_efficiency) D_T_90(position_diff_efficiency) sum(D_R_90) sum(D_T_90) D_R(position_diff_efficiency) D_T(position_diff_efficiency) sum(D_R) sum(D_T)];
            fprintf(fid, '%0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record);
            fprintf('%0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record);
        end
        fclose(fid);
        
        if save_m==1
            copyfile('data4.txt',save_m_as)
        end
        
        a=load('data4.txt');
        plot(a(:,1),a(:,3));
    case 5
        record_0=[];
        record_90=[];
        record_psi=[];
        psi0=minimum_angle_psi;
        if use_dispersion==1
            setup_dispersion; 
        end
        initialization_input_file;
        rcwa_procedure;
        
        fid = fopen('data5.txt', 'wt');
        for cycle=minimum_angle_psi:step:maximum_angle_psi            
            psi0=cycle;  
            
            save_all_results;
            diffraction_efficiency;
            position_diff_efficiency=1+(diffraction_efficiency_order(2)+N_Y)+(2*N_Y+1)*(diffraction_efficiency_order(1)+N_Y);
            record=[cycle D_R(position_diff_efficiency) D_T(position_diff_efficiency) sum(D_R) sum(D_T)];
            fprintf(fid, '%0f %12.10f %12.10f %12.10f %12.10f\n', record);
            fprintf('%0f %12.10f %12.10f %12.10f %12.10f\n', record);
        end
        fclose(fid);
        
        if save_m==1
            copyfile('data5.txt',save_m_as)
        end
        
        a=load('data5.txt');
        plot(a(:,1),a(:,3));
    case 6
        record_0=[];
        record_90=[];
        record_psi=[];
        fid = fopen('data6.txt', 'wt');
        for cycle=minimum_angle_phi:step:maximum_angle_phi            
            phi0=cycle;            
            if use_dispersion==1
                setup_dispersion;
            end
            initialization_input_file;
            rcwa_procedure;
            
            save_all_results;
            
            position_diff_efficiency=1+(diffraction_efficiency_order(2)+N_Y)+(2*N_Y+1)*(diffraction_efficiency_order(1)+N_Y);
            record=[cycle D_R_0(position_diff_efficiency) D_T_0(position_diff_efficiency) sum(D_R_0) sum(D_T_0) D_R_90(position_diff_efficiency) D_T_90(position_diff_efficiency) sum(D_R_90) sum(D_T_90) D_R(position_diff_efficiency) D_T(position_diff_efficiency) sum(D_R) sum(D_T)];
            fprintf(fid, '%0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record);
            fprintf('%0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record);
        end
        fclose(fid);
        
        if save_m==1
            copyfile('data6.txt',save_m_as)
        end
        
        a=load('data6.txt');
        plot(a(:,1),a(:,3));

end