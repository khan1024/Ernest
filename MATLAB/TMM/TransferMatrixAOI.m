
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

