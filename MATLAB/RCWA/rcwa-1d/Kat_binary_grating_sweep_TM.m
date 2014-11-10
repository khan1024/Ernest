%Runs 1-D grating RCWA (2D simulation) called main_binary.m
%by Katherine Han 7-26-14
pause
close all
clear all

%Input variables

lambda_min = 0.75; %sets minimum wavelength for sweep in micrometers
lambda_max = 0.75;   %sets maximum wavelength for sweep in micrometers
lambda_delta = 0.75; %sets delta for wavelength steps for sweep in micrometers

theta0_min = 0; %sets minimum angle of incidence for sweep in degrees
theta0_max = 0;   %sets maximum angle of incidence for sweep in degrees
theta0_delta = 5; %sets delta for angle of incidence steps for sweep in degrees

polarization_min = 1; %sets minimum polarization for sweep (1=TM)
polarization_max = 1;   %sets maximum polarization for sweep (2=TE)
polarization_delta = 1; %sets delta for polarization steps (must be 1)

Lambda_min = 0.4; %sets minimum grating period for sweep in micrometers
Lambda_max = 2.4;   %sets maximum grating period for sweep in micrometers
Lambda_delta = 0.1; %sets delta for grating period steps for sweep in micrometers

thickness_total_min = 0.75; %sets minimum grating thickness for sweep in micrometers
thickness_total_max = 0.75;   %sets maximum grating thickness for sweep in micrometers
thickness_total_delta = 0.5; %sets delta for grating thickness steps for sweep in micrometers

dutycycle_min = 0.5; %sets minimum wavelength for sweep in micrometers
dutycycle_max = 0.5;   %sets maximum wavelength for sweep in micrometers
dutycycle_delta = 0.4; %sets delta for wavelength steps for sweep in micrometers

n1 = 2; %sets index of refraction for superstrate (usually air) and grating troughs
n3 = 1;   %sets index of refraction for substrate and grating teeth material

lambda = lambda_min:lambda_delta:lambda_max;
theta0 = theta0_min:theta0_delta:theta0_max;
polarization = polarization_min:polarization_delta:polarization_max;
Lambda = Lambda_min:Lambda_delta:Lambda_max;
thickness_total = thickness_total_min:thickness_total_delta:thickness_total_max;
dutycycle = dutycycle_min:dutycycle_delta:dutycycle_max;

for l = 1:length(lambda)
    for t = 1:length(theta0)
        for p = 1:length(polarization)
            for L = 1:length(Lambda)
                for thick = 1:length(thickness_total)
                    for d = 1:length(dutycycle)
                         [diff_R, diff_R_angles, diff_T, diff_T_angles] = main_binary(lambda(l),theta0(t),polarization(p),Lambda(L),thickness_total(thick),dutycycle(d),n1,n3);
                          %output results to these matrices, notice curly
                          %brackets allow one to output a matrix into each cell
                         R_Orders_Output{l,t,p,L,thick,d} = diff_R;
                         T_Orders_Output{l,t,p,L,thick,d} = diff_T;
                         R_Orders_Angles{l,t,p,L,thick,d} = diff_R_angles;
                         T_Orders_Angles{l,t,p,L,thick,d} = diff_T_angles;
                    end
                end
            end
        end
    end
    l
end


%%%-------------Put your results somewhere useful----------%%%

% put all reflected zeroth order efficiencies for the first vallue of 
% wavelength  (lambda) = 0.4 um for a 1.2 um thick grating for TE polarization into
% a matrix that sweeps over AOI (theta0) and period (Lambda)

% for t = 1:length(theta0)
%     for L = 1:length(Lambda)
%         l = 1; %400 nm light
%         p = 1; %TE polarization
%         thick = 3; %1.2 um thickness grating
%         row = find(R_Orders_Output{l,t,p,L,thick}(:,1) == 0);%gives row number of order 0 in the first matrix calculated
%         R400TE1_2_0th(t,L) = R_Orders_Output{l,t,p,L,thick}(row,2);
%     end
% end

% for t = 1:length(theta0)
%     for L = 1:length(Lambda)
%         l = 1; %400 nm light
%         p = 1; %TE polarization
%         thick = 1; %1.2 um thickness grating
%         row = find(R_Orders_Angles{l,t,p,L,thick}(:,1) == 0);%gives row number of order 0 in the first matrix calculated
%         R400TE1_2_0th_Angles(t,L) = R_Orders_Angles{l,t,p,L,thick}(row,2);
%     end
% end

% 
% for t = 1:length(theta0)
%     for L = 1:length(Lambda)
%         l = 1; %400 nm light
%         p = 1; %TE polarization
%         thick = 3; %1.2 um thickness grating
%         d = 2; % duty cycle is 0.4
%         row = find(R_Orders_Output{l,t,p,L,thick,d}(:,1) == -1);%gives row number of order 0 in the first matrix calculated
%         R400TE1_2_neg1st(t,L) = R_Orders_Output{l,t,p,L,thick,d}(row,2);
%     end
% end
% 
% R400TE1_2_pos1st = zeros(length(theta0),length(Lambda));
% for t = 1:length(theta0)
%     for L = 1:length(Lambda)
%         l = 1; %400 nm light
%         p = 1; %TE polarization
%         thick = 3; %1.2 um thickness grating
%         d = 2; %duty cycle 0.4
%         row = find(R_Orders_Output{l,t,p,L,thick,d}(:,1) == 1);%gives row number of order 0 in the first matrix calculated
%         if (isempty(row) == true)
%             R400TE1_2_pos1st(t,L) = 0;
%         else
%             R400TE1_2_pos1st(t,L) = R_Orders_Output{l,t,p,L,thick,d}(row,2);
%         end
%     end
% end
    
%make a vector that is an independent axis for all of the diffraction
%orders that  you might want to graph
Order_Angles = -90:1:90;

% set up data to graph efficiencies of each order at their angle of
% diffraction for reflected and transmitted orders while sweeping along
% AOI, thickness, period, wavelength, or duty cycle.

% To sweep over angles of incidence for one set of the other parameters

R_All_Orders_Theta0_Sweep = zeros(length(Order_Angles),length(theta0));
R_All_Orders_Lambda_Sweep = zeros(length(Order_Angles),length(Lambda));
R_All_Orders_lambda_Sweep = zeros(length(Order_Angles),length(lambda));
R_All_Orders_thickness_total_Sweep = zeros(length(Order_Angles),length(thickness_total));
R_All_Orders_dutycycle_Sweep = zeros(length(Order_Angles),length(dutycycle));

T_All_Orders_Theta0_Sweep = zeros(length(Order_Angles),length(theta0));
T_All_Orders_Lambda_Sweep = zeros(length(Order_Angles),length(Lambda));
T_All_Orders_lambda_Sweep = zeros(length(Order_Angles),length(lambda));
T_All_Orders_thickness_total_Sweep = zeros(length(Order_Angles),length(thickness_total));
T_All_Orders_dutycycle_Sweep = zeros(length(Order_Angles),length(dutycycle));


% make all the angles negative, because the way they are plotted is
% confusing otherwise.  how about positive angles are headed toward the
% right, negative angles are headed toward the left.  To do this multiply
% the angles by -1 as they are put into these matrices

%% to sweep over period
for L=1:length(Lambda)
    t = 1; % AOI=0
    l = 1; %750 nm light 
    p = 1; %TE polarization
    thick = 1; %0.75 um thickness
    d = 1; %0.5 duty cycle
    for x=1:length(R_Orders_Output{l,t,p,L,thick,d}(:,1))
        tempmatrix(x,1)= (-1*(R_Orders_Angles{l,t,p,L,thick,d}(x,2)));
        tempmatrix(x,2)= R_Orders_Output{l,t,p,L,thick,d}(x,2);
        tempmatrix(x,1)=round(tempmatrix(x,1));
        anglerow = find(Order_Angles == tempmatrix(x,1));
        R_All_Orders_Lambda_Sweep(anglerow,L) = tempmatrix(x,2);
    end
    disp(tempmatrix)
    clear tempmatrix
end
    

for L=1:length(Lambda)
    t = 1; %AOI=0
    l = 1; %750 nm light
    p = 1; %TE polarization
    d = 1; %0.75 um thickness grating
    thick = 1;
    for x=1:length(T_Orders_Output{l,t,p,L,thick,d}(:,1))
        tempmatrix(x,1)= (-1*(T_Orders_Angles{l,t,p,L,thick,d}(x,2)));
        tempmatrix(x,2)= T_Orders_Output{l,t,p,L,thick,d}(x,2);
        tempmatrix(x,1)=round(tempmatrix(x,1));
        anglerow = find(Order_Angles == tempmatrix(x,1));
        T_All_Orders_Lambda_Sweep(anglerow,L) = tempmatrix(x,2);
    end
    disp (tempmatrix)
        clear tempmatrix
end

%% to sweep over thickness of grating
for thick=1:length(thickness_total)
    t = 1; % AOI=0
    l = 1; %750 nm light 
    p = 1; %TE polarization
    L = 1; %  um thickness
    d = 1; %0.5 duty cycle
    for x=1:length(R_Orders_Output{l,t,p,L,thick,d}(:,1))
        tempmatrix(x,1)= (-1*(R_Orders_Angles{l,t,p,L,thick,d}(x,2)));
        tempmatrix(x,2)= R_Orders_Output{l,t,p,L,thick,d}(x,2);
        tempmatrix(x,1)=round(tempmatrix(x,1));
        anglerow = find(Order_Angles == tempmatrix(x,1));
        R_All_Orders_thickness_total_Sweep(anglerow,thick) = tempmatrix(x,2);
    end
    disp(tempmatrix)
    clear tempmatrix
end
    

for thick=1:length(thickness_total)
    t = 1; %AOI=0
    l = 1; %750 nm light
    p = 1; %TE polarization
    d = 1; %0.75 um thickness grating
    L = 1;
    for x=1:length(T_Orders_Output{l,t,p,L,thick,d}(:,1))
        tempmatrix(x,1)= (-1*(T_Orders_Angles{l,t,p,L,thick,d}(x,2)));
        tempmatrix(x,2)= T_Orders_Output{l,t,p,L,thick,d}(x,2);
        tempmatrix(x,1)=round(tempmatrix(x,1));
        anglerow = find(Order_Angles == tempmatrix(x,1));
        T_All_Orders_thickness_total_Sweep(anglerow,thick) = tempmatrix(x,2);
    end
    disp (tempmatrix)
        clear tempmatrix
end

%% to sweep over angle of incidence of grating
for t=1:length(theta0)
    thick = 1; % thickness of grating
    l = 1; %750 nm light 
    p = 1; %TE polarization
    L = 1; %  um thickness
    d = 1; % duty cycle
    for x=1:length(R_Orders_Output{l,t,p,L,thick,d}(:,1))
        tempmatrix(x,1)= (-1*(R_Orders_Angles{l,t,p,L,thick,d}(x,2)));
        tempmatrix(x,2)= R_Orders_Output{l,t,p,L,thick,d}(x,2);
        tempmatrix(x,1)=round(tempmatrix(x,1));
        anglerow = find(Order_Angles == tempmatrix(x,1));
        R_All_Orders_thata0_Sweep(anglerow,t) = tempmatrix(x,2);
    end
    disp(tempmatrix)
    clear tempmatrix
end
    

for t=1:length(theta0)
    thick = 1; %thickness of grating
    l = 1; %750 nm light
    p = 1; %TE polarization
    d = 1; % duty cycle
    L = 1; %period
    for x=1:length(T_Orders_Output{l,t,p,L,thick,d}(:,1))
        tempmatrix(x,1)= (-1*(T_Orders_Angles{l,t,p,L,thick,d}(x,2)));
        tempmatrix(x,2)= T_Orders_Output{l,t,p,L,thick,d}(x,2);
        tempmatrix(x,1)=round(tempmatrix(x,1));
        anglerow = find(Order_Angles == tempmatrix(x,1));
        T_All_Orders_theta0_Sweep(anglerow,t) = tempmatrix(x,2);
    end
    disp (tempmatrix)
        clear tempmatrix
end

%% to sweep over input wavelength
for l=1:length(lambda)
    thick = 1; % thickness of grating
    t = 1; 
    p = 1; %TE polarization
    L = 1; %  um thickness
    d = 1; % duty cycle
    for x=1:length(R_Orders_Output{l,t,p,L,thick,d}(:,1))
        tempmatrix(x,1)= (-1*(R_Orders_Angles{l,t,p,L,thick,d}(x,2)));
        tempmatrix(x,2)= R_Orders_Output{l,t,p,L,thick,d}(x,2);
        tempmatrix(x,1)=round(tempmatrix(x,1));
        anglerow = find(Order_Angles == tempmatrix(x,1));
        R_All_Orders_lambda_Sweep(anglerow,l) = tempmatrix(x,2);
    end
    disp(tempmatrix)
    clear tempmatrix
end
    

for l=1:length(lambda)
    thick = 1; %thickness of grating
    t = 1; 
    p = 1; %TE polarization
    d = 1; % duty cycle
    L = 1; %period
    for x=1:length(T_Orders_Output{l,t,p,L,thick,d}(:,1))
        tempmatrix(x,1)= (-1*(T_Orders_Angles{l,t,p,L,thick,d}(x,2)));
        tempmatrix(x,2)= T_Orders_Output{l,t,p,L,thick,d}(x,2);
        tempmatrix(x,1)=round(tempmatrix(x,1));
        anglerow = find(Order_Angles == tempmatrix(x,1));
        T_All_Orders_lambda_Sweep(anglerow,l) = tempmatrix(x,2);
    end
    disp (tempmatrix)
        clear tempmatrix
end

%% to sweep over duty cycle of grating
for d=1:length(dutycycle)
    thick = 1; % thickness of grating
    t = 1; 
    p = 1; %TE polarization
    L = 1; %  um thickness
    l = 1; %
    for x=1:length(R_Orders_Output{l,t,p,L,thick,d}(:,1))
        tempmatrix(x,1)= (-1*(R_Orders_Angles{l,t,p,L,thick,d}(x,2)));
        tempmatrix(x,2)= R_Orders_Output{l,t,p,L,thick,d}(x,2);
        tempmatrix(x,1)=round(tempmatrix(x,1));
        anglerow = find(Order_Angles == tempmatrix(x,1));
        R_All_Orders_dutycycle_Sweep(anglerow,d) = tempmatrix(x,2);
    end
    disp(tempmatrix)
    clear tempmatrix
end
    

for d=1:length(dutycycle)
    thick = 1; %thickness of grating
    t = 1; 
    p = 1; %TE polarization
    l = 1; 
    L = 1; %period
    for x=1:length(T_Orders_Output{l,t,p,L,thick,d}(:,1))
        tempmatrix(x,1)= (-1*(T_Orders_Angles{l,t,p,L,thick,d}(x,2)));
        tempmatrix(x,2)= T_Orders_Output{l,t,p,L,thick,d}(x,2);
        tempmatrix(x,1)=round(tempmatrix(x,1));
        anglerow = find(Order_Angles == tempmatrix(x,1));
        T_All_Orders_dutycycle_Sweep(anglerow,d) = tempmatrix(x,2);
    end
    disp (tempmatrix)
        clear tempmatrix
end

% figure
% mesh (Order_Angles,Lambda,T_All_Orders_Theta_Sweep')
% xlabel ('Transmitted Diffraction Angle (deg)')
% ylabel ('Grating Period (micrometers)')
% zlabel ('Transmision Efficiency')

%% make a bunch of nice contour graphs that show all the orders of reflection and transmission for each sweep

v = 0:0.01:1;

% -----------Figures that sweep over period, Lambda
figure
contourf (Order_Angles,Lambda,T_All_Orders_Theta_Sweep',v,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlabel ('Transmitted Diffraction Angle (deg)','FontSize',18)
ylabel ('Grating Period (micrometers)','FontSize',18)
zlabel ('Transmision Efficiency','FontSize',18)
title ({'Transmitted Diffraction Efficiency for 1D Binary Grating for TE Light','with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1'},'FontSize',18)
print(gcf,'-dpng','-r600',['Transmitted Diffraction Efficiency for 1D Binary Grating for TE Light with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1.png'])

figure
contourf (Order_Angles,Lambda,R_All_Orders_Theta_Sweep',v,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlabel ('Reflected Diffraction Angle (deg)','FontSize',18)
ylabel ('Grating Period (micrometers)','FontSize',18)
zlabel ('Reflection Efficiency','FontSize',18)
title ({'Reflected Diffraction Efficiency for 1D Binary Grating for TE Light','with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1'},'FontSize',18)
print(gcf,'-dpng','-r600',['Reflected Diffraction Efficiency for 1D Binary Grating for TE Light with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1.png'])

% ------------Figures that sweep over grating thickness, thickness_total
figure
contourf (Order_Angles,Lambda,T_All_Orders_Theta_Sweep',v,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlabel ('Transmitted Diffraction Angle (deg)','FontSize',18)
ylabel ('Grating Period (micrometers)','FontSize',18)
zlabel ('Transmision Efficiency','FontSize',18)
title ({'Transmitted Diffraction Efficiency for 1D Binary Grating for TE Light','with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1'},'FontSize',18)
print(gcf,'-dpng','-r600',['Transmitted Diffraction Efficiency for 1D Binary Grating for TE Light with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1.png'])

figure
contourf (Order_Angles,Lambda,R_All_Orders_Theta_Sweep',v,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlabel ('Reflected Diffraction Angle (deg)','FontSize',18)
ylabel ('Grating Period (micrometers)','FontSize',18)
zlabel ('Reflection Efficiency','FontSize',18)
title ({'Reflected Diffraction Efficiency for 1D Binary Grating for TE Light','with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1'},'FontSize',18)
print(gcf,'-dpng','-r600',['Reflected Diffraction Efficiency for 1D Binary Grating for TE Light with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1.png'])


% ------------Figures that sweep over wavelength input (lambda)
figure
contourf (Order_Angles,Lambda,T_All_Orders_Theta_Sweep',v,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlabel ('Transmitted Diffraction Angle (deg)','FontSize',18)
ylabel ('Grating Period (micrometers)','FontSize',18)
zlabel ('Transmision Efficiency','FontSize',18)
title ({'Transmitted Diffraction Efficiency for 1D Binary Grating for TE Light','with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1'},'FontSize',18)
print(gcf,'-dpng','-r600',['Transmitted Diffraction Efficiency for 1D Binary Grating for TE Light with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1.png'])

figure
contourf (Order_Angles,Lambda,R_All_Orders_Theta_Sweep',v,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlabel ('Reflected Diffraction Angle (deg)','FontSize',18)
ylabel ('Grating Period (micrometers)','FontSize',18)
zlabel ('Reflection Efficiency','FontSize',18)
title ({'Reflected Diffraction Efficiency for 1D Binary Grating for TE Light','with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1'},'FontSize',18)
print(gcf,'-dpng','-r600',['Reflected Diffraction Efficiency for 1D Binary Grating for TE Light with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1.png'])

% ------------Figures that sweep over grating duty cycle, dutycycle
figure
contourf (Order_Angles,Lambda,T_All_Orders_Theta_Sweep',v,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlabel ('Transmitted Diffraction Angle (deg)','FontSize',18)
ylabel ('Grating Period (micrometers)','FontSize',18)
zlabel ('Transmision Efficiency','FontSize',18)
title ({'Transmitted Diffraction Efficiency for 1D Binary Grating for TE Light','with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1'},'FontSize',18)
print(gcf,'-dpng','-r600',['Transmitted Diffraction Efficiency for 1D Binary Grating for TE Light with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1.png'])

figure
contourf (Order_Angles,Lambda,R_All_Orders_Theta_Sweep',v,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlabel ('Reflected Diffraction Angle (deg)','FontSize',18)
ylabel ('Grating Period (micrometers)','FontSize',18)
zlabel ('Reflection Efficiency','FontSize',18)
title ({'Reflected Diffraction Efficiency for 1D Binary Grating for TE Light','with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1'},'FontSize',18)
print(gcf,'-dpng','-r600',['Reflected Diffraction Efficiency for 1D Binary Grating for TE Light with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1.png'])

% ------------Figures that sweep over angle of incidence, theta0
figure
contourf (Order_Angles,Lambda,T_All_Orders_Theta_Sweep',v,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlabel ('Transmitted Diffraction Angle (deg)','FontSize',18)
ylabel ('Grating Period (micrometers)','FontSize',18)
zlabel ('Transmision Efficiency','FontSize',18)
title ({'Transmitted Diffraction Efficiency for 1D Binary Grating for TE Light','with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1'},'FontSize',18)
print(gcf,'-dpng','-r600',['Transmitted Diffraction Efficiency for 1D Binary Grating for TE Light with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1.png'])

figure
contourf (Order_Angles,Lambda,R_All_Orders_Theta_Sweep',v,'edgecolor','none')
set(gca,'FontSize',14)
colorbar
xlabel ('Reflected Diffraction Angle (deg)','FontSize',18)
ylabel ('Grating Period (micrometers)','FontSize',18)
zlabel ('Reflection Efficiency','FontSize',18)
title ({'Reflected Diffraction Efficiency for 1D Binary Grating for TE Light','with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1'},'FontSize',18)
print(gcf,'-dpng','-r600',['Reflected Diffraction Efficiency for 1D Binary Grating for TE Light with0.75 um thickness 0.5 duty cycle for 750 nm light from n=2 to n=1.png'])
