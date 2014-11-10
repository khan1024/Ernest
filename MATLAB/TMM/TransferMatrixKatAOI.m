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
lambda=700:1200; % Wavelengths over which reflectance is calculated

% The layer thicknesses are in nanometers.



thicknesses = [0 200 0];  % thickness of each corresponding layer in nm (thickness of the first layer is irrelivant)

% theta = [0:60/2/pi];
theta = [0:pi/64:31*pi/64];

%------------END USER INPUT PARAMETERS SPECIFICATION-------------------

%assume non-wavelength dependent index of refraction
n1=1.5; %index of substrate, glass
n2=1.4; %index of first layer of ARC
% n3=1.1; %index of second layer of ARC
% n4=1.2; %index of third layer of ARC
% n5=1.3; %index of fourth layer of ARC
% n6=1.4; %index of fifth layer of ARC
n3=1;
% n8=1; %index of air

RI = [n1 n2 n3];

n=[];   %Initializes a blank matrix
for index = 1:length(lambda) 
    n(:,index) = RI;
end



% Calculate transfer matrices, and field at each wavelength and position
t(1)=0;
Rs = zeros(length(lambda),length(theta));
Rp = zeros(length(lambda),length(theta));
R = zeros(length(lambda),length(theta));

    
for l = 1:length(lambda)
    for th = 1:length(theta)
    % Calculate transfer matrices for incoherent reflection and transmission at the first interface
    Ss=I_mat_s_polarized(n(1,l),n(2,l),theta(th));
    % Continue to add phase angle matrix and intensity matrix for each
    % additional layer
        for matindex=2:(length(thicknesses)-1)
            Ss=Ss*L_mat_aoi(n(matindex,l), thicknesses(1,matindex),lambda(l),theta(th))*I_mat_s_polarized(n(matindex,l),n(matindex+1,l),theta(th));
        end
        %Calculate the reflectance and transmittance for each wavelength
    Rs(l,th)=(abs(Ss(2,1)/Ss(1,1)))^2; %JAP Vol 86 p.487 Eq 9 Power Reflection from layers other than substrate
     % Calculate transfer matrices for incoherent reflection and transmission at the first interface
    Sp=I_mat_p_polarized(n(1,l),n(2,l),theta(th));
    % Continue to add phase angle matrix and intensity matrix for each
    % additional layer
        for matindex=2:(length(thicknesses)-1)
            Sp=Sp*L_mat_aoi(n(matindex,l), thicknesses(1,matindex),lambda(l),theta(th))*I_mat_p_polarized(n(matindex,l),n(matindex+1,l),theta(th));
        end
        %Calculate the reflectance and transmittance for each wavelength
    Rp(l,th)=(abs(Sp(2,1)/Sp(1,1)))^2; %JAP Vol 86 p.487 Eq 9 Power Reflection from layers other than substrate
    R(l,th)=(Rp(l,th)+Rs(l,th))/2;
    end
end


% close all
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

figure
hold on
plot(lambda,R(:,1)*100,'-k','LineWidth',2)
% plot(lambda,R(:,2)*100,'-b','LineWidth',2)
% plot(lambda,R(:,3)*100,'-g','LineWidth',2)
% plot(lambda,R(:,17)*100,'-r','LineWidth',2)

% plot(ExperimentalWaveLength, ExperimentalReflectance,'-g','LineWidth',2) %Experimental data

title('Percent of Light absorbed or reflected');
xlabel('Wavelength (nm)');
ylabel('Light Intensity Percent');
% ylim([0 5])

theta_deg=theta.*180./pi;
figure
% mesh(theta_deg,lambda,R)
mesh(lambda, theta_deg,R')

