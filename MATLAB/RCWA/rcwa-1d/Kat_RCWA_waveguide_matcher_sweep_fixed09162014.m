%Runs 1-D grating RCWA (2D simulation) called main_binary.m
%by Katherine Han 7-26-14
% pause
close all
clear all
clc

%% Input variables

gratingcase = 6;
% 0 = 'Binary' 
% 1 = 'Blazed type 1' 
% 2 = 'Blazed type 2' 
% 3 = 'Sinusoidal' 
% 4 = 'Triangular' 
% 5 = 'Parabolic' 
% 6 = 'General Triangular' 
% 7 = 'Trapezoidal' 
% 8 = 'Tetragonal' 
% 9 = 'Gaussian' 
% 10 = 'Arbitrary'
gratingcaselabel = ['Binary';'Blaze1';'Blaze2';'Sin   ';'Triang';'Parabo';'GenTri';'Trapez';'Tetrag';'Gaussi';'Arbitr'];

lambda_min = 0.7; %sets minimum wavelength for sweep in micrometers
lambda_max = 1.2;   %sets maximum wavelength for sweep in micrometers
lambda_delta = 0.05; %sets delta for wavelength steps for sweep in micrometers

theta0_min = 0; %sets minimum angle of incidence for sweep in degrees
theta0_max = 0;   %sets maximum angle of incidence for sweep in degrees
theta0_delta = 1; %sets delta for angle of incidence steps for sweep in degrees

polarization_min = 2; %sets minimum polarization for sweep (1=TM)
polarization_max = 2;   %sets maximum polarization for sweep (2=TE)
polarization_delta = 1; %sets delta for polarization steps (must be 1)

Lambda_min = 1.65; %sets minimum grating period for sweep in micrometers///////////////
Lambda_max = 1.65;   %sets maximum grating period for sweep in micrometers
Lambda_delta = 0.1; %sets delta for grating period steps for sweep in micrometers

thickness_total_min = 1.65;
thickness_total_max = 1.65;   %sets maximum grating thickness for sweep in micrometers
thickness_total_delta = 0.1; %sets delta for grating thickness steps for sweep in micrometers

dutycycle_min = 1; %sets minimum wavelength for sweep in micrometers
dutycycle_max = 1;   %sets maximum wavelength for sweep in micrometers
dutycycle_delta = 0.2; %sets delta for wavelength steps for sweep in micrometers/////////////

nair = 1;
n1 = 2; %sets index of refraction for superstrate (usually air) and grating troughs
n3 = 67;   %sets index of refraction for substrate and grating teeth material

WGd = 6.05; %waveguide thickness in micrometers//////////
WGn = n1; %waveguide index of refraction for grating on bottom////////

% NOTE: plottype must match the variable over which you are seeping with this run
plottype = 3;
% plottype = 1 for Lambda sweep
% plottype = 2 for thickness_total sweep
% plottype = 3 for lambda sweep
% plottype = 4 for dutycycle sweep
% plottype = 5 for theta0 sweep

%% 
lambda = lambda_min:lambda_delta:lambda_max;
theta0 = theta0_min:theta0_delta:theta0_max;
polarization = polarization_min:polarization_delta:polarization_max;
Lambda = Lambda_min:Lambda_delta:Lambda_max;
thickness_total = thickness_total_min:thickness_total_delta:thickness_total_max;
dutycycle = dutycycle_min:dutycycle_delta:dutycycle_max;


total_iterations = length(dutycycle)*length(theta0)*length(polarization)*length(Lambda)*length(thickness_total)*length(lambda);
current_iteration = 1;

for  d = 1:length(dutycycle)
    for t = 1:length(theta0)
        for p = 1:length(polarization)
            for L = 1:length(Lambda)
                for thick = 1:length(thickness_total)
                    for l = 1:length(lambda)
                         [diff_R, diff_R_angles, diff_T, diff_T_angles] = main_binary(lambda(l),theta0(t),polarization(p),Lambda(L),thickness_total(thick),dutycycle(d),n1,n3,gratingcase);
                          %output results to these matrices, notice curly
                          %brackets allow one to output a matrix into each cell
                         R_Orders_Output{l,t,p,L,thick,d} = diff_R;
                         T_Orders_Output{l,t,p,L,thick,d} = diff_T;
                         R_Orders_Angles{l,t,p,L,thick,d} = diff_R_angles;
                         T_Orders_Angles{l,t,p,L,thick,d} = diff_T_angles;
                         
                         my_progress = current_iteration/total_iterations*100;
                         my_progress_output = sprintf('%2.3f percent done',my_progress);
                         disp(my_progress_output);
                         current_iteration = current_iteration+1;
                    end
                end
            end
        end
    end
end



%make a vector that is an independent axis for all of the diffraction
%orders that  you might want to graph
Order_Angles = -90:1:90;

% set up data to graph efficiencies of each order at their angle of
% diffraction for reflected and transmitted orders while sweeping along
% AOI, thickness, period, wavelength, or duty cycle.

% To sweep over angles of incidence for one set of the other parameters


T_All_Orders_lambda_Sweep = zeros(length(Order_Angles),length(lambda));


% make all the transmitted angles negative, because the way they are plotted is
% confusing otherwise.  how about positive angles are headed toward the
% right, negative angles are headed toward the left (incoming light is from the left).  To do this multiply
% the transmitted angles by -1 as they are put into these matrices


for  d = 1:length(dutycycle)
    for t = 1:length(theta0)
        for p = 1:length(polarization)
            for L = 1:length(Lambda)
                for thick = 1:length(thickness_total)
                    R_All_Orders_lambda_Sweep = zeros(length(Order_Angles),length(lambda));

                    for l = 1:length(lambda)
                        for x=1:length(R_Orders_Output{l,t,p,L,thick,d}(:,1))

                            tempmatrix(x,1)= ((R_Orders_Angles{l,t,p,L,thick,d}(x,2)));
                            tempmatrix(x,2)= R_Orders_Output{l,t,p,L,thick,d}(x,2);
                            tempmatrix(x,1)=round(tempmatrix(x,1));
                            anglerow = find(Order_Angles == tempmatrix(x,1));
                            R_All_Orders_lambda_Sweep(anglerow,l) = tempmatrix(x,2);
                        end
                        R_Orders_lambda_all_param_sweep{t,p,L,thick,d} = R_All_Orders_lambda_Sweep;
                        clear tempmatrix
                    end
                end
            end
        end
    end
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
   % disp (tempmatrix)
        clear tempmatrix
end


D = zeros(1,length(WGd));
F = zeros(1,length(WGd));
Criteria = zeros(length(theta0),length(polarization),length(thickness_total),length(dutycycle),length(WGd));

for  d = 1:length(dutycycle)
    for t = 1:length(theta0)
        for p = 1:length(polarization)
            for L = 1:length(Lambda)
                for thick = 1:length(thickness_total)
                       for WGthickness = 1:length(WGd)
                            WG_Modes = zeros(length(Order_Angles),length(lambda));
                            Product = zeros(length(Order_Angles),length(lambda));
                            
                            for lam = 1:length(lambda)
                                M = 2*WGd(WGthickness)*WGn/lambda(lam)*(n1^2-nair^2)^0.5;
                                for m = 1:floor(M)
                                    theta_m=90-180/pi*asin(m*lambda(lam)/2/WGd(WGthickness)/WGn);
                                    angle2 = ceil(theta_m);
                                    a = find(Order_Angles == angle2);
                                    Diff = a - 91;
                                    b = 91- Diff;
                                    WG_Modes(a,lam)=2;
                                    WG_Modes(a+1,lam) = 1;
                                    WG_Modes(a-1,lam) = 1;
                                    WG_Modes(b,lam)=2;
                                    WG_Modes(b+1,lam) = 1;
                                    WG_Modes(b-1,lam) = 1;
                                end
                            end
                            Product = R_Orders_lambda_all_param_sweep{t,p,L,thick,d} .* WG_Modes;
                            D(t,p,L,thick,d,WGthickness) = sum(sum(Product));
                            F(t,p,L,thick,d,WGthickness) = length(lambda);
                            Criteria(t,p,L,thick,d,WGthickness) = D(t,p,L,thick,d,WGthickness)./F(t,p,L,thick,d,WGthickness);
                       end 
                    end
              end
        end
    end
end

[value, location] = max(Criteria(:))
clear a1 b1 c1 d1 e1 f1 
[a1,b1,c1,d1,e1,f1] = ind2sub(size(Criteria),location);
angle_of_incidence = theta0(a1)
TM1TE2 = polarization(b1)
grating_period = Lambda(c1)
grating_thickness = thickness_total(d1)
duty_cycle = dutycycle(e1)
Wave_Guide_thickness = WGd(f1)

v = 0:0.01:1;
w = 0:0.5:2;

graph1 = figure;
 contourf (Order_Angles,lambda.*1000,WG_Modes',w,'edgecolor','none')
        set(gca,'FontSize',14)
        colorbar
        xlabel ('Mode Ray angle (degrees)','FontSize',16)
        ylabel ('Wavelength (nm)','FontSize',16)
        zlabel ('Modes','FontSize',16)
        title (['Waveguide mode ray angles for thickness of ' num2str(WGd(f1))...
        ' um and index of ' num2str(n1)],'FontSize',14)
    saveas(graph1,['.\figures\Waveguide mode ray angles for thickness of ' num2str(WGd(f1))...
        ' um and index of ' num2str(n1) '.png'])

%         title ({['Transmitted Diffraction Efficiency for 1D ' gratingcaselabel(gratingcase+1,:) ' Grating Period Sweep for TE Light with '],[ num2str(thickness_total(1,thick)) ' um thickness ' num2str(dutycycle(1,d)) ' duty cycle for ' num2str(lambda(1,l)) ' um light from n=' num2str(n1) ' to n=' num2str(n3) ', AOI=' num2str(theta0(1,t))]},'FontSize',14)




%% make a bunch of nice contour graphs that show all the orders of reflection and transmission for each sweep



switch plottype

    case 3
        
        % ------------Figures that sweep over wavelength input (lambda) h5 and h6
%         h5=figure;
%         contourf (Order_Angles,lambda,T_All_Orders_lambda_Sweep',v,'edgecolor','none')
%         set(gca,'FontSize',14)
%         colorbar
%         xlabel ('Transmitted Diffraction Angle (deg)','FontSize',16)
%         ylabel ('Wavelength (micrometers)','FontSize',16)
%         zlabel ('Transmision Efficiency','FontSize',16)
%         title ({['Transmitted Diffraction Efficiency for 1D ' gratingcaselabel(gratingcase+1,:) ' Grating Wavelength Sweep for TE Light with '],[ num2str(thickness_total(1,thick)) ' um thickness, ' num2str(dutycycle(1,d)) ' duty cycle, ' num2str(Lambda(1,L)) ' um period from n=' num2str(n1) ' to n=' num2str(n3) ', AOI=' num2str(theta0(1,t))]},'FontSize',14)
%         saveas(h5,['.\figures\Transmitted Diffraction Efficiency for 1D ' gratingcaselabel(gratingcase+1,:) ' Grating Wavelength Sweep for TE Light with ' num2str(thickness_total(1,thick)) ' um thickness, ' num2str(dutycycle(1,d)) ' duty cycle, ' num2str(Lambda(1,L)) ' um period from n=' num2str(n1) ' to n=' num2str(n3)  ', AOI=' num2str(theta0(1,t)) '.png'])

        h6=figure;
        contourf (Order_Angles,lambda,R_All_Orders_lambda_Sweep',v,'edgecolor','none')
        set(gca,'FontSize',14)
        colorbar
        xlabel ('Reflected Diffraction Angle (deg)','FontSize',16)
        ylabel ('Wavelength (micrometers)','FontSize',16)
        zlabel ('Reflection Efficiency','FontSize',16)
        title ({['Reflected Diffraction Efficiency for 1D ' gratingcaselabel(gratingcase+1,:)...
            ' Grating Wavelength Sweep for TE Light with '],[ num2str(thickness_total(1,thick))...
            ' um thickness, ' num2str(dutycycle(1,d)) ' duty cycle, ' num2str(Lambda(1,L))...
            ' um period from n=' num2str(n1) ' to n=' num2str(n3) ', AOI=' num2str(theta0(1,t))]},'FontSize',14)
        saveas(h6,['.\figures\Reflected Diffraction Efficiency for 1D ' gratingcaselabel(gratingcase+1,:) ...
            ' Grating Wavelength Sweep for TE Light with ' num2str(thickness_total(1,thick)) ' um thickness, '...
            num2str(dutycycle(1,d)) ' duty cycle, ' num2str(Lambda(1,L)) ' um period from n=' num2str(n1) ' to n=' ...
            num2str(n3)  ', AOI=' num2str(theta0(1,t)) '.png'])


    otherwise
        warning('Unexpected plot type. No plot was created.')
end

