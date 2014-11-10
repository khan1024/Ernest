clear all

load_file='test_13.mat';
load_psi_efficiencies=1; % 1=0, 2=90, 3=psi 
order_R_or_T=2; % 1-reflection orders, 2-transmision orders
diffraction_efficiency_order=[0,0]; 
export_diffraction_efficiency_order=1; % 1-ano, 2-ne;
save_m_export_as='pokus_export.txt';

%---------------------------------------------------------
load(load_file);

if load_psi_efficiencies==1
    record=record_0;
elseif load_psi_efficiencies==2
    record=record_90;
elseif load_psi_efficiencies==3
    record=record_psi;
end

n_sq=(2*N_X+1)*(2*N_Y+1);

position_diff_efficiency_R=1+1+(diffraction_efficiency_order(2)+N_Y)+(2*N_Y+1)*(diffraction_efficiency_order(1)+N_Y); % 1+    
position_diff_efficiency_T=n_sq+1+1+(diffraction_efficiency_order(2)+N_Y)+(2*N_Y+1)*(diffraction_efficiency_order(1)+N_Y); % n_sq+1

axis_x=record(1,:);

if order_R_or_T==1
    osa_y=record(position_diff_efficiency_R,:);
else
    osa_y=record(position_diff_efficiency_T,:);
    %osa_y=sum(record_0(n_sq+2:end,:));
end

figure;
axes('FontSize',17);
plot(axis_x,osa_y,'Marker','.','MarkerEdgeColor','r','MarkerFaceColor','b','MarkerSize',15);

if export_diffraction_efficiency_order==1
    fid = fopen(save_m_export_as, 'wt');
    for j=1:1:length(axis_x)
        record_export=[axis_x(j) record_0(position_diff_efficiency_R,j) record_0(position_diff_efficiency_T,j) sum(record_0(2:n_sq+1,j)) sum(record_0(n_sq+2:end,j)) record_90(position_diff_efficiency_R,j) record_90(position_diff_efficiency_T,j) sum(record_90(2:n_sq+1,j)) sum(record_90(n_sq+2:end,j)) record_psi(position_diff_efficiency_R,j) record_psi(position_diff_efficiency_T,j) sum(record_psi(2:n_sq+1,j)) sum(record_psi(n_sq+2:end,j))];
        fprintf(fid, '%0f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', record_export);
    end
    fclose(fid);
end

