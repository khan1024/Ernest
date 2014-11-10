% % % % % Note: to make the contour graph this must be in the same folder as contourcreator1.m

% Copyright 2010 George F. Burkhard, Eric T. Hoke, Stanford University

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% This code has been heavily modified by Kat Han 1-25-14

% This program calculates the reflectance at in interface that contains one or more thin film layers of ARC

% The procedure was adapted from J. Appl. Phys Vol 86, No. 1 (1999) p.487
% and JAP 93 No. 7 p. 3693.

% George Burkhard and Eric Hoke February 5, 2010
% When citing this work, please refer to:
%
% G. F. Burkhard, E. T. Hoke, M. D. McGehee, Adv. Mater., 22, 3293.
% Accounting for Interference, Scattering, and Electrode Absorption to Make
% Accurate Internal Quantum Efficiency Measurements in Organic and Other 
% Thin Solar Cells


%------------BEGIN USER INPUT PARAMETERS SPECITIFCATION---------------
%

for target_three_quarter_wave = 200:50:200
    for number_layers = 28:2:28 % make this one higher than the number of mirror layers you want to include substrate. Air will be included later, too
        for n_top = 2
            for n_small = 1.5
                for n_large = 2.5
                    for n_bottom = 1
                        
% target_three_quarter_wave = a;
% number_layers = b;
% n_top = c;
% n_small = d;
% n_large = e;
% n_bottom = f;

startlambda = 400;
endlambda = 1300;
lambdadelta = 1;
lambda=startlambda:lambdadelta:endlambda; % Wavelengths over which reflectance is calculated
% The layer thicknesses are in nanometers

%  total_thickness = 1400:50:1400;

total_thickness = [0];

theta = 0:pi/64:31*pi/64;

n=zeros(number_layers + 1);   %Initializes the matrix

n(number_layers + 1) = n_top;
n(1) = n_bottom;

for x = 1:floor(number_layers/2)
    n(2*x) = n_large;
end

for y = 1:floor(number_layers/2)-1
    n(2*y+1) = n_small;
end

% n_step_x = 0:pi/(2*number_layers):pi/2;
% n_step=zeros(length(n_step_x));
% 
% for x=1:length(n_step_x)
%     n_step(x)=(sin(n_step_x(x)))^4;
%     n(number_layers +2-x)=n_min+n_step(x)*(n_max-n_min);
% end


% thicknesses = zeros(number_layers + 1, length(total_thickness));  
thicknesses=zeros(number_layers + 1);   %Initializes the matrix
for u = 1:floor(number_layers/2)
    thicknesses(2*u,1) = target_three_quarter_wave/n_large;
end

for v = 1:floor(number_layers/2)-1
    thicknesses(2*v+1,1) = target_three_quarter_wave/n_small;
end


% for y = 1:length(total_thickness)
% 
%     for x = 2:number_layers
%         layer_thickness = total_thickness(y) / number_layers;
%         thicknesses(x,y) = layer_thickness;
%     end
% end
% thickness of each corresponding layer in nm (thickness of the first layer is irrelivant)





%------------END USER INPUT PARAMETERS SPECIFICATION-------------------

% %assume non-wavelength dependent index of refraction
% n1=1.52; %index of substrate, glass
% n2=1.52; %index of first layer of ARC
% n3=1.52; %index of second layer of ARC
% n4=1.22; %index of third layer of ARC
% n5=1.0; %index of fourth layer of ARC
% n6=1.0; %index of fifth layer of ARC
% n7=1;
% % n8=1; %index of air
% 
% RI = [n1 n2 n3 n4 n5 n6 n7];
% 
% n=[];   %Initializes a blank matrix
% for index = 1:length(lambda) 
%     n(:,index) = RI;
% end



% Calculate transfer matrices, and field at each wavelength and position
% t(1)=0;
Rs = zeros(length(lambda),length(theta),length(total_thickness));
Rp = zeros(length(lambda),length(theta),length(total_thickness));
R = zeros(length(lambda),length(theta),length(total_thickness));

for t = 1:length(total_thickness)    
    for l = 1:length(lambda)
        for th = 1:length(theta)
        % Calculate transfer matrices for incoherent reflection and transmission at the first interface
        Ss=I_mat_s_polarized(n(1),n(2),theta(th));
        % Continue to add phase angle matrix and intensity matrix for each
        % additional layer
            for matindex=2:(length(thicknesses(:,1))-1)
                Ss=Ss*L_mat_aoi(n(matindex), thicknesses(matindex,t),lambda(l),theta(th))*I_mat_s_polarized(n(matindex),n(matindex+1),theta(th));
            end
            %Calculate the reflectance and transmittance for each wavelength
        Rs(l,th,t)=(abs(Ss(2,1)/Ss(1,1)))^2; %JAP Vol 86 p.487 Eq 9 Power Reflection from layers other than substrate
         % Calculate transfer matrices for incoherent reflection and transmission at the first interface
        Sp=I_mat_p_polarized(n(1),n(2),theta(th));
        % Continue to add phase angle matrix and intensity matrix for each
        % additional layer
            for matindex=2:(length(thicknesses)-1)
                Sp=Sp*L_mat_aoi(n(matindex), thicknesses(matindex,t),lambda(l),theta(th))*I_mat_p_polarized(n(matindex),n(matindex+1),theta(th));
            end
            %Calculate the reflectance and transmittance for each wavelength
        Rp(l,th,t)=(abs(Sp(2,1)/Sp(1,1)))^2; %JAP Vol 86 p.487 Eq 9 Power Reflection from layers other than substrate
        R(l,th,t)=(Rp(l,th,t)+Rs(l,th,t))/2;
        end
    end
     disp('Help!  Im trapped in a MATLAB program!');
end

WeightedR = zeros(length(theta),length(total_thickness));

load('AM15.mat')
StartIndex = find(AM15(:,1)==startlambda);
EndIndex = find(AM15(:,1)==endlambda);
totalAM1_5 = sum(AM15(StartIndex:EndIndex,3));
AM1_5_WeightedR = zeros (length(theta),length(total_thickness));
R_1 = zeros(length(lambda),length(theta),length(total_thickness));
Sum_R = zeros(length(theta),length(total_thickness));

for th = 1:length(theta)
    for t = 1:length(total_thickness)
        for l = 1:length(lambda)
            R_1(l,th,t) = R(l,th,t).*AM15(l+StartIndex-1,3);
% AM1_5_WeightedR (th,t) = sum((R(:,th,t).*AM1_5(:)))./totalAM1_5;
        end
        Sum_R(th,t)=sum(R_1(:,th,t));
        AM1_5_WeightedR (th,t) = Sum_R(th,t)./totalAM1_5;
    end
end

close all


% close all

opengl software
% Plots the total reflection expected from the interface
% (useful for comparing with experimentally measured reflection spectrum)


% %Experimental data for TMOSb path
% strTMOSbPath = 'C:\Users\Kat\Documents\0-Chang Research\Data\Graph of all TMOS b depositions.xlsx';
% 
% %Reads in all experimental data
% [ExperimentalTMOSData,~,~] = xlsread(strTMOSbPath,'Sheet1');
% 
% %Extracts out the data we're interested in (reflectance data, experimental)
% ExperimentalReflectance = ExperimentalTMOSData(8:408,156);
% 
% %Extracts out the experimental data's X-axis
% ExperimentalWaveLength = ExperimentalTMOSData(8:408,1);

% figure
% hold on
% plot(lambda,R(:,1)*100,'-k','LineWidth',2)
% plot(lambda,R(:,2)*100,'-b','LineWidth',2)
% plot(lambda,R(:,3)*100,'-g','LineWidth',2)
% plot(lambda,R(:,17)*100,'-r','LineWidth',2)

% % plot(ExperimentalWaveLength, ExperimentalReflectance,'-g','LineWidth',2) %Experimental data

% title('Percent of Light absorbed or reflected');
% xlabel('Wavelength (nm)');
% ylabel('Light Intensity Percent');
% % ylim([0 5])

theta_deg=theta.*180./pi;

% figure
% % mesh(theta_deg,lambda,R)
% mesh(theta_deg,total_thickness,AM1_5_WeightedR'*100)
% title('Average %R Weighted by AM1.5 for sin^4 profile, n=1.38 to n=2');
% xlabel('Angle of Incidence (deg)');
% ylabel('Total Thickness of ARC (nm)');
% zlabel('Average %R in AM1.5');

% contourf(theta_deg,total_thickness,AM1_5_WeightedR'*100)
% title('Average %R Weighted by AM1.5 for sin^4 profile, n=1 to n=2');
% xlabel('Angle of Incidence (deg)');
% ylabel('Total Thickness of ARC (nm)');
% zlabel('Average %R in AM1.5');
% colorbar
%------------------
% contourcreator1(theta_deg,total_thickness,AM1_5_WeightedR'*100)

% Rlambda = zeros(length(lambda),length(total_thickness));
% for lam = 1:length(lambda)
%     for thick = 1:length(total_thickness)
%     Rlambda(lam,thick) = R(lam,1,thick);
%     end
% end
% mesh(Rlambda*100)

%---------------------


% figure
% % mesh(theta_deg,lambda,R)
% mesh(lambda,theta_deg,R(:,:,1)'*100)
% title('%R for sin^4 profile, n=1.38 to n=2, 99 nm thick, 32 layers');
% xlabel('Wavelength (nm)');
% ylabel('Angle of Incidence (deg)');
% zlabel('%R');

% figure
% % mesh(theta_deg,lambda,R)
% plot(theta_deg,R(1,:,1)'*100)
% title('%R for sin^4 profile, n=1.38 to n=2, 99 nm thick, 32 layers');
% xlabel('Angle of Incidence (deg)');
% ylabel('%R');
increment = 0:100;

overall_thickness = sum(thicknesses(:,:));

h1=figure;
mesh(theta_deg,lambda,R(:,:,1)*100)
set(gca,'FontSize',14)
colorbar
xlim([0 90])
ylim([400 1300])
zlim([0 100])
xlabel('Angle of Incidence (deg)','FontSize',16);
ylabel('Wavelength (nm)','FontSize',16);
zlabel('%R','FontSize',16);
title({['Reflectance for Dichroic Mirror of ' num2str(number_layers -1) ...
    ' layers, n=' num2str(n_large) ' and n=' num2str(n_small)],[...
    ' at ' num2str(round(thicknesses(2,1))) ' and ' num2str(round(thicknesses(3,1))) ...
    ' nm layers and overall thickness of '...
    num2str(round(overall_thickness)) ' nm from n=' num2str(n_top) ' to n=' num2str(n_bottom)']},'FontSize',14)
saveas(h1,['.\figures\Reflectance for Dichroic Mirror of ' num2str(number_layers -1) ...
    ' layers, n= ' num2str(floor(n_large)) '_' num2str(100*(n_large - floor(n_large)))...
    ' and n= ' num2str(floor(n_small)) '_' num2str(100*(n_small - floor(n_small)))...
    ' at ' num2str(round(thicknesses(2,1))) ' and ' num2str(round(thicknesses(3,1))) ' nm layers and overall thickness of '...
    num2str(round(overall_thickness)) ' nm from n=' num2str(floor(n_top)) '_' num2str(100*(n_top - floor(n_top)))...
    ' to n=' num2str(floor(n_bottom)) '_' num2str(100*(n_bottom - floor(n_bottom))) '.png'])

h2=figure;
contourf(theta_deg,lambda,R(:,:,1)*100,increment,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlim([0 90])
ylim([400 1300])
zlim([0 100])
xlabel('Angle of Incidence (deg)','FontSize',16);
ylabel('Wavelength (nm)','FontSize',16);
zlabel('%R','FontSize',16);
title({['Reflectance for Dichroic Mirror of ' num2str(number_layers -1) ...
    ' layers, n= ' num2str(n_large) ' and n= ' num2str(n_small)],[...
    ' at ' num2str(round(thicknesses(2,1))) ' and ' num2str(round(thicknesses(3,1))) ...
    ' nm layers and overall thickness of '...
    num2str(round(overall_thickness)) ' nm from n=' num2str(n_top) ' to n=' num2str(n_bottom)]},'FontSize',14)
saveas(h2,['.\figures\Reflectance for Dichroic Mirror of ' num2str(number_layers -1) ...
    ' layers, n= ' num2str(floor(n_large)) '_' num2str(100*(n_large - floor(n_large)))...
    ' and n= ' num2str(floor(n_small)) '_' num2str(100*(n_small - floor(n_small)))...
    ' at ' num2str(round(thicknesses(2,1))) ' and ' num2str(round(thicknesses(3,1))) ' nm layers and overall thickness of '...
    num2str(round(overall_thickness)) ' nm from n=' num2str(floor(n_top)) '_' num2str(100*(n_top - floor(n_top)))...
    ' to n=' num2str(floor(n_bottom)) '_' num2str(100*(n_bottom - floor(n_bottom))) 'contour.png'])

%  clear all
clear target_three_quarter_wave
clear number_layers 
clear n_top 
clear n_small 
clear n_large 
clear n_bottom 
clear n
clear thicknesses

                    end
                end
            end
        end
    end
end


% figure
% % mesh(theta_deg,lambda,R)
% contour(lambda,total_thickness,R(:,1,:)*100)
% contour(Rlambda.*100)


% figure
% % mesh(theta_deg,lambda,R)
% mesh(lambda,theta_deg,R(:,:,3)'*100)
% 
% figure
% % mesh(theta_deg,lambda,R)
% mesh(lambda,theta_deg,R(:,:,6)'*100)


%Function for contour graph
