%% Initializes the script

%Clears the screen
clc;

%Clears all the variables
clear all;

%Closes all windows 
close all;

%% File Path Info

%Reads the Excel spreadsheet for parse locations and file settings
[~,~,RawParse] = xlsread('DataPath.xlsx','Raw FDTD Data');

%Constructs arrDataPath from Excel spreadsheet
arrDataPath = RawParse(2:(size(RawParse,1)),:);

%Determines the size of the data to parse
intDataParseSize = size(arrDataPath,1);

%Prints starting message
fprintf('\n\nGetting started!\n\nMake sure to close all other unnecessary programs...this will take some memory.\nClick on the command window and press Ctrl+C to terminate the job.\n');

%Loads the AM 1.5 data
load('AM1.5.mat')

%Loops through the sets of data to parse
for intDataParseIndex = 1:size(arrDataPath,1)
    %Prints progress of overall parsing
    fprintf('\n%i of %i datasets started.\n\n',intDataParseIndex, size(arrDataPath,1));
 
    
    %Length of data
    intDataXLength = cell2mat(arrDataPath(intDataParseIndex,9));
    intDataYLength = cell2mat(arrDataPath(intDataParseIndex,10));

    %Where to store the data
    strAnalyzedDataPath = [char(arrDataPath(intDataParseIndex,12)) '\'];

    %"Reflected" file name stub
    strFileReflected = ['\' char(arrDataPath(intDataParseIndex,13))];

    %"Transmitted" file name stub
    strFileTransmitted = ['\' char(arrDataPath(intDataParseIndex,14))];
    
    %Loads the data
    load([strAnalyzedDataPath char(arrDataPath(intDataParseIndex,2)) '.mat']);
    
    %Creates the wavelength vector
    c = 299792458;
    WaveLength = c./freq.*1e9;

    %Figures out the minimum wavelength
    minWave = min(WaveLength);
    maxWave = max(WaveLength);

    %Truncates the AM1.5 data to just the pertinent wavelengths
    AMTruncated = find(AM15(:,1) > minWave);
    AMStart = AMTruncated(1);
    AMTruncated = find(AM15(:,1) < maxWave);
    AMEnd = max(AMTruncated);
    AM15 = AM15(AMStart:AMEnd,:);

    %Clears useless variables
    clear 'c' 'AMEnd' 'AMStart' 'AMTruncated' 'min*' 'max*'; 

    %Builds the transmitted data array
    txArray(intDataParseIndex,:) = txAvg;
    
    %Builds the received data array
    refArray(intDataParseIndex,:) = refAvg;
    
    %Calculates the percent transmitted
    percentTxArray(intDataParseIndex,:) = txAvg ./ txAvgFS;
    
    %Calculates the percent reflected
    percentReflArray(intDataParseIndex,:) = refAvg ./ txAvgFS;
    
    %Calculates the average percent reflected 
    avgReflArray (intDataParseIndex) = sum(refAvg ./ txAvgFS) / length(txAvgFS);
    
    %===================================================
    % Performs weighted average
    %===================================================
    
    %Calculates the percent reflected
    PercentReflected = refAvg ./ txAvgFS;

    %Flips the data
    percentRefl = fliplr(PercentReflected);

    %Performs linear interpolation
    percentReflInterp = interp1(WaveLength,percentRefl,AM15(:,1));

    %Multiplies data by AM1.5 spectrum
    AMScaled = AM15(:,3) .* percentReflInterp;
    AMReference = AM15(:,3);
    WaveLenAM15 = AM15(:,1);
    
    %Calculates the weighted average percent reflected 
    avgWeightedReflArray (intDataParseIndex) = sum(AMScaled) / sum(AM15(:,3));
    
    %Builds the legend
    strLegend(intDataParseIndex) = arrDataPath(intDataParseIndex,2);
end

%Clears the misc. data
clear ('arrDataPath', 'intDataParseIndex', 'strAnalyzedDataPath','refAvg','txAvg');

%% =======================================
%  Makes figures
%=========================================

%Saves the data
%save (['M:\Kat FDTD Data\Analyzed Data\Frequency Response\Cumulative Data, Made on ' char(date)])
% % %=========================375 nm Cylinders and thin film================================
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2) % Bare 1.5 material
% plot(WaveLength, percentReflArray(186,:).*100,'--g','LineWidth',2) % test
% 
% title ('5 Layer Thin Film With and Without Imperfections yj')
% % legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['test.png'])
% %print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\test.png'])
% 
%% Plots data
% %=========================375 nm Cylinders and thin film================================
figure
hold on
plot(WaveLength, percentReflArray(4,:).*100,'-k','LineWidth',2) % Bare 1.5 material
plot(WaveLength, percentReflArray(5,:).*100,'--g','LineWidth',2) % 5 layer thin film, 100 nm/layer, 1.1, 1.15, 1.2, 1.25, 1.3
% plot(WaveLength, percentReflArray(169,:).*100,'--b','LineWidth',2)% 5 layer thin film with 50 random imperfections 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(170,:).*100,'-b','LineWidth',2)% 5 layer thin film with 10 random imperfections 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(171,:).*100,':b','LineWidth',2)% 5 layer thin film with 100 random imperfections 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(172,:).*100,':r','LineWidth',2)% 5 layer thin film with 100 spherical random imperfections 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(173,:).*100,'--m','LineWidth',2)% 5 layer thin film with 100 spherical random imperfections 50 nm diameter 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(174,:).*100,'-r','LineWidth',2)% 5 layer thin film with 200 spherical random imperfections 100 nm diameter 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(175,:).*100,'-g','LineWidth',2)% 5 layer thin film with 100 spherical random imperfections 100 nm diameter 1.1 to 1.5 RI

title ('5 Layer Thin Film With and Without Imperfections')
% legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
 ylim([0 5])

%--legend('Free Space','Transmitted','Reflected');
xlabel ('Wavelength (nm)');
ylabel ('Percent Reflected');

%Sets the font to 16, bold and saving the figure
set(findall(gcf,'-property','FontSize'),'FontSize',16) 
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
print(gcf,'-dpng','-r600',['testing20141004.png'])
%print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\5LayerImperfections_3-13-14.png'])

% % %=========================5 layer thin film with 100 spherical imperfections Different Sizes================================
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2) % Bare 1.5 material
% plot(WaveLength, percentReflArray(168,:).*100,'--k','LineWidth',2) % 5 layer thin film, 100 nm/layer, no imperfections 1.1, 1.15, 1.2, 1.25, 1.3
% plot(WaveLength, percentReflArray(173,:).*100,'-b','LineWidth',2)% 5 layer thin film with 100 spherical random imperfections 50 nm diameter 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(172,:).*100,'-c','LineWidth',2)% 5 layer thin film with 100 spherical random imperfections 100 nm diameter 1.1 to 1.3 RI
% % plot(WaveLength, percentReflArray(175,:).*100,'--r','LineWidth',2)% 5 layer thin film with 100 spherical random imperfections 100 nm diameter 1.1 to 1.5 RI
% % plot(WaveLength, percentReflArray(176,:).*100,'-b','LineWidth',2)% 5 layer thin film with 100 spherical random imperfections 150 nm diameter 1.1 to 1.5 RI
% plot(WaveLength, percentReflArray(177,:).*100,'-g','LineWidth',2)% 5 layer thin film with 100 spherical random imperfections 150 nm diameter 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(178,:).*100,'-r','LineWidth',2)% 5 layer thin film with 150 spherical random imperfections 50, 100, and 150 nm diameter 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(179,:).*100,'--r','LineWidth',2)% 5 layer thin film with 175 spherical random imperfections 50, 100, 150, and 200 nm diameter 1.1 to 1.3 RI
% % plot(WaveLength, percentReflArray(180,:).*100,':r','LineWidth',2)% 5 105 nm layer thin film (1.19 to 1.35) with 175 cube random imperfections 50, 100, 150, and 200 nm diameter 1.1 to 1.5 RI
% plot(WaveLength, percentReflArray(181,:).*100,':r','LineWidth',2)% 5 105 nm layer thin film (1.19 to 1.35) with 150 cube random imperfections 50, 100, and 150 nm diameter 1.1 to 1.5 RI
% plot(WaveLength, percentReflArray(182,:).*100,'-.r','LineWidth',2)% 5 105 nm layer thin film (1.19 to 1.35) with 125 cube random imperfections 50, 100, and 150 nm diameter 1.1 to 1.5 RI
% plot(WaveLength, percentReflArray(183,:).*100,':m','LineWidth',2)% 5 105 nm layer thin film (1.19 to 1.35) with 100 cube random imperfections 50, 100 nm diameter 1.1 to 1.5 RI
% plot(WaveLength, percentReflArray(184,:).*100,'-m','LineWidth',2)% 5 105 nm layer thin film (1.19 to 1.35) with no imperfections 
% 
% title ('5 Layer Thin Film With and Without 100 Spherical Imperfections')
% % legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
%  ylim([0 5])
%  xlim([400 800])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\5Layer100SphericalImperfections_n_1_to_1-3_3-13-14.png'])
% % %=========================5 Layer Thin Film With and Without 100 Spherical Imperfections================================
% 
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
% 
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2) % Bare 1.5 material
% plot(WaveLength, percentReflArray(184,:).*100,'--k','LineWidth',2)% 5 105 nm layer thin film (1.19 to 1.35) with no imperfections 
% plot(WaveLength, percentReflArray(183,:).*100,'-b','LineWidth',2)% 5 105 nm layer thin film (1.19 to 1.35) with 100 cube random imperfections 50, 100 nm diameter 1.1 to 1.5 RI
% plot(WaveLength, percentReflArray(182,:).*100,'-c','LineWidth',2)% 5 105 nm layer thin film (1.19 to 1.35) with 125 cube random imperfections 50, 100, and 150 nm diameter 1.1 to 1.5 RI
% plot(WaveLength, percentReflArray(181,:).*100,'-g','LineWidth',2)% 5 105 nm layer thin film (1.19 to 1.35) with 150 cube random imperfections 50, 100, and 150 nm diameter 1.1 to 1.5 RI
% plot(WaveLength, percentReflArray(185,:).*100,':g','LineWidth',2)% 5 105 nm layer thin film (1.19 to 1.35) with 150 sphere random imperfections 50, 100, and 150 nm diameter 1.1 to 1.5 RI
% 
% plot(ExperimentalWaveLength, ExperimentalReflectance,'-r','LineWidth',2) %Experimental data
% 
% title ('5 Layer Thin Film With and Without 100 Spherical Imperfections')
% % legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
% ylim([0 5])
% xlim([400 800])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\5Layer_n119_n135_GradientCubeImperfections_n_1_to_1-5_3-13-14.png'])
% % %=========================5 layer thin film with 100 nm sphere defects different numbers================================
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2) % Bare 1.5 material
% plot(WaveLength, percentReflArray(168,:).*100,'--k','LineWidth',2) % 5 layer thin film, 100 nm/layer, 1.1, 1.15, 1.2, 1.25, 1.3
% plot(WaveLength, percentReflArray(172,:).*100,':r','LineWidth',2)% 5 layer thin film with 100 spherical random imperfections 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(174,:).*100,'-b','LineWidth',2)% 5 layer thin film with 200 spherical random imperfections 100 nm diameter 1.1 to 1.3 RI
% plot(WaveLength, percentReflArray(175,:).*100,'-r','LineWidth',2)% 5 layer thin film with 100 spherical random imperfections 100 nm diameter 1.1 to 1.5 RI
% 
% title ('5 Layer Thin Film With and Without Imperfections')
% % legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\5LayerImperfections100nmSpheres_3-13-14.png'])
% %=========================375 nm Cylinders and thin film================================
% figure
% hold on
% plot(WaveLength, percentReflArray(53,:).*100,'-g','LineWidth',2) % Cylinder full width 200 nm wide
% plot(WaveLength, percentReflArray(54,:).*100,'--g','LineWidth',2) % Cylinder full width 100 nm wide
% plot(WaveLength, percentReflArray(55,:).*100,':g','LineWidth',2)% Cylinder full width 50 nm wide
% plot(WaveLength, percentReflArray(56,:).*100,'-r','LineWidth',2)% Cylinder 80% width 200 nm wide
% plot(WaveLength, percentReflArray(57,:).*100,'--r','LineWidth',2)% Cylinder 80% width 100 nm wide
% plot(WaveLength, percentReflArray(58,:).*100,':r','LineWidth',2) % Cylinder 80% width 50 nm wide
% % plot(WaveLength, percentReflArray(59,:).*100,'-b','LineWidth',2) % Cylinder half width 200 nm wide
% plot(WaveLength, percentReflArray(60,:).*100,'--b','LineWidth',2) % Cylinder half width 100 nm wide
% plot(WaveLength, percentReflArray(61,:).*100,':b','LineWidth',2) % Cylinder half width 50 nm wide
% plot(WaveLength, percentReflArray(52,:).*100,'-c','LineWidth',2) % Thin film 375 nm n=1.22
% 
% 
% title ('Cylinders 375 nm long and thin film, n=1.22, 375 nm thick')
% % legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Cylinders375nm.png'])

% %=========================125 nm Cylinders and thin film================================
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2) % no ARC
% plot(WaveLength, percentReflArray(52,:).*100,'-c','LineWidth',2) % Thin film 375 nm n=1.22
% plot(WaveLength, percentReflArray(53,:).*100,'-g','LineWidth',2) % Cylinder full width 200 nm wide
% plot(WaveLength, percentReflArray(54,:).*100,'--g','LineWidth',2) % Cylinder full width 100 nm wide
% plot(WaveLength, percentReflArray(55,:).*100,':g','LineWidth',2)% Cylinder full width 50 nm wide
% plot(WaveLength, percentReflArray(56,:).*100,'-r','LineWidth',2)% Cylinder 80% width 200 nm wide
% plot(WaveLength, percentReflArray(57,:).*100,'--r','LineWidth',2)% Cylinder 80% width 100 nm wide
% plot(WaveLength, percentReflArray(58,:).*100,':r','LineWidth',2) % Cylinder 80% width 50 nm wide
% % plot(WaveLength, percentReflArray(59,:).*100,'-b','LineWidth',2) % Cylinder half width 200 nm wide
% plot(WaveLength, percentReflArray(60,:).*100,'--b','LineWidth',2) % Cylinder half width 100 nm wide
% plot(WaveLength, percentReflArray(61,:).*100,':b','LineWidth',2) % Cylinder half width 50 nm wide
% 
% 
% title ('Cylinders 375 nm long and thin film, n=1.22, 375 nm thick')
% % legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Cylinders375nm.png'])

% %=======================All Cones=======================================
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2)
% plot(WaveLength, percentReflArray(98,:).*100,'b-','LineWidth',2) % Cones 200 nm wide 100 nm long
% plot(WaveLength, percentReflArray(99,:).*100,'g-','LineWidth',2) % Cones 200 nm wide 200 nm long
% plot(WaveLength, percentReflArray(100,:).*100,'r-','LineWidth',2) % Cones 200 nm wide 400 nm long
% plot(WaveLength, percentReflArray(101,:).*100,'c-','LineWidth',2) % Cones 200 nm wide 800 nm long
% 
% plot(WaveLength, percentReflArray(102,:).*100,'b--','LineWidth',2)% Cones 100 nm wide 100 nm long
% plot(WaveLength, percentReflArray(103,:).*100,'g--','LineWidth',2)% Cones 100 nm wide 200 nm long
% plot(WaveLength, percentReflArray(104,:).*100,'r--','LineWidth',2)% Cones 100 nm wide 400 nm long
% plot(WaveLength, percentReflArray(105,:).*100,'c--','LineWidth',2)% Cones 100 nm wide 800 nm long
% 
% plot(WaveLength, percentReflArray(106,:).*100,'b:','LineWidth',2)% Cones 50 nm wide 100 nm long
% plot(WaveLength, percentReflArray(107,:).*100,'g:','LineWidth',2)% Cones 50 nm wide 200 nm long
% plot(WaveLength, percentReflArray(108,:).*100,'r:','LineWidth',2)% Cones 50 nm wide 400 nm long
% plot(WaveLength, percentReflArray(109,:).*100,'c:','LineWidth',2)% Cones 50 nm wide 800 nm long
% title ('Cones')
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Cones_color_boldReflected.png'])
% 
% %=======================All Cones 5 nm mesh=======================================
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2)
% plot(WaveLength, percentReflArray(86,:).*100,'b-','LineWidth',2) % Cones 200 nm wide 100 nm long
% plot(WaveLength, percentReflArray(87,:).*100,'g-','LineWidth',2) % Cones 200 nm wide 200 nm long
% plot(WaveLength, percentReflArray(88,:).*100,'r-','LineWidth',2) % Cones 200 nm wide 400 nm long
% plot(WaveLength, percentReflArray(89,:).*100,'c-','LineWidth',2) % Cones 200 nm wide 800 nm long
% 
% plot(WaveLength, percentReflArray(90,:).*100,'b--','LineWidth',2)% Cones 100 nm wide 100 nm long
% plot(WaveLength, percentReflArray(91,:).*100,'g--','LineWidth',2)% Cones 100 nm wide 200 nm long
% plot(WaveLength, percentReflArray(92,:).*100,'r--','LineWidth',2)% Cones 100 nm wide 400 nm long
% plot(WaveLength, percentReflArray(93,:).*100,'c--','LineWidth',2)% Cones 100 nm wide 800 nm long
% 
% plot(WaveLength, percentReflArray(94,:).*100,'b:','LineWidth',2)% Cones 50 nm wide 100 nm long
% plot(WaveLength, percentReflArray(95,:).*100,'g:','LineWidth',2)% Cones 50 nm wide 200 nm long
% plot(WaveLength, percentReflArray(96,:).*100,'r:','LineWidth',2)% Cones 50 nm wide 400 nm long
% plot(WaveLength, percentReflArray(97,:).*100,'c:','LineWidth',2)% Cones 50 nm wide 800 nm long
% title ('Cones 5 nm Mesh')
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Cones_color_boldReflected5nmMesh.png'])
% %================== Pyramids 5 nm mesh=======================================
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2)
% plot(WaveLength, percentReflArray(74,:).*100,'b-','LineWidth',2) % Pyramids 200 nm wide 100 nm long
% plot(WaveLength, percentReflArray(75,:).*100,'g-','LineWidth',2) % Pyramids 200 nm wide 200 nm long
% plot(WaveLength, percentReflArray(76,:).*100,'r-','LineWidth',2) % Pyramids 200 nm wide 400 nm long
% plot(WaveLength, percentReflArray(77,:).*100,'c-','LineWidth',2) % Pyramids 200 nm wide 800 nm long
% 
% plot(WaveLength, percentReflArray(78,:).*100,'b--','LineWidth',2)% Pyramids 100 nm wide 100 nm long
% plot(WaveLength, percentReflArray(79,:).*100,'g--','LineWidth',2)% Pyramids 100 nm wide 200 nm long
% plot(WaveLength, percentReflArray(80,:).*100,'r--','LineWidth',2)% Pyramids 100 nm wide 400 nm long
% plot(WaveLength, percentReflArray(81,:).*100,'c--','LineWidth',2)% Pyramids 100 nm wide 800 nm long
% 
% plot(WaveLength, percentReflArray(82,:).*100,'b:','LineWidth',2)% Pyramids 50 nm wide 100 nm long
% plot(WaveLength, percentReflArray(83,:).*100,'g:','LineWidth',2)% Pyramids 50 nm wide 200 nm long
% plot(WaveLength, percentReflArray(84,:).*100,'r:','LineWidth',2)% Pyramids 50 nm wide 400 nm long
% plot(WaveLength, percentReflArray(85,:).*100,'c:','LineWidth',2)% Pyramids 50 nm wide 800 nm long
% 
% title ('Pyramids 5 nm mesh')
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\PyramidsWideTallSweepsReflected5nmMesh.png'])
% 
% %================== Pyramids=======================================
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2)
% plot(WaveLength, percentReflArray(62,:).*100,'b:','LineWidth',2) % Pyramids 200 nm wide 100 nm long
% plot(WaveLength, percentReflArray(63,:).*100,'g:','LineWidth',2) % Pyramids 200 nm wide 200 nm long
% plot(WaveLength, percentReflArray(64,:).*100,'r:','LineWidth',2) % Pyramids 200 nm wide 400 nm long
% plot(WaveLength, percentReflArray(65,:).*100,'c:','LineWidth',2) % Pyramids 200 nm wide 800 nm long
% 
% plot(WaveLength, percentReflArray(66,:).*100,'b--','LineWidth',2)% Pyramids 100 nm wide 100 nm long
% plot(WaveLength, percentReflArray(67,:).*100,'g--','LineWidth',2)% Pyramids 100 nm wide 200 nm long
% plot(WaveLength, percentReflArray(68,:).*100,'r--','LineWidth',2)% Pyramids 100 nm wide 400 nm long
% plot(WaveLength, percentReflArray(69,:).*100,'c--','LineWidth',2)% Pyramids 100 nm wide 800 nm long
% 
% plot(WaveLength, percentReflArray(70,:).*100,'b-','LineWidth',2)% Pyramids 50 nm wide 100 nm long
% plot(WaveLength, percentReflArray(71,:).*100,'g-','LineWidth',2)% Pyramids 50 nm wide 200 nm long
% plot(WaveLength, percentReflArray(72,:).*100,'r-','LineWidth',2)% Pyramids 50 nm wide 400 nm long
% plot(WaveLength, percentReflArray(73,:).*100,'c-','LineWidth',2)% Pyramids 50 nm wide 800 nm long
% 
% title ('Pyramids')
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\PyramidsWideTallSweepsReflected.png'])
% 
% % %=========================125 nm Cylinders and thin film================================
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2) % no ARC
% plot(WaveLength, percentReflArray(35,:).*100,'-m','LineWidth',2) % Thin film 125 nm n=1.22
% plot(WaveLength, percentReflArray(36,:).*100,'-g','LineWidth',2) % Cylinder full width 200 nm wide
% plot(WaveLength, percentReflArray(37,:).*100,'--g','LineWidth',2) % Cylinder full width 100 nm wide
% plot(WaveLength, percentReflArray(38,:).*100,':g','LineWidth',2)% Cylinder full width 50 nm wide
% plot(WaveLength, percentReflArray(39,:).*100,'-r','LineWidth',2)% Cylinder 80% width 200 nm wide
% plot(WaveLength, percentReflArray(40,:).*100,'--r','LineWidth',2)% Cylinder 80% width 100 nm wide
% plot(WaveLength, percentReflArray(41,:).*100,':r','LineWidth',2) % Cylinder 80% width 50 nm wide
% plot(WaveLength, percentReflArray(42,:).*100,'-c','LineWidth',2) % Cylinder half width 200 nm wide
% plot(WaveLength, percentReflArray(43,:).*100,'--c','LineWidth',2) % Cylinder half width 100 nm wide
% plot(WaveLength, percentReflArray(44,:).*100,':c','LineWidth',2) % Cylinder half width 50 nm wide
% % plot(WaveLength, percentReflArray(45,:).*100,'-b','LineWidth',2)% Cylinder quarter width 200 nm wide
% % plot(WaveLength, percentReflArray(46,:).*100,'--b','LineWidth',2)% Cylinder quarter width 100 nm wide
% % plot(WaveLength, percentReflArray(47,:).*100,':b','LineWidth',2)% Cylinder quarter width 50 nm wide
% 
% 
% title ('Cylinders 125 nm long and thin film, n=1.22, 125 nm thick')
% % legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Cylinders125nm.png'])
% 
% %=========================Hollow Spheres================================
% figure
% hold on
% plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2) % no ARC
% plot(WaveLength, percentReflArray(130,:).*100,'-b','LineWidth',2) % 1 layer hollow spheres 50 nm diameter 5nm shell thickness
% plot(WaveLength, percentReflArray(133,:).*100,'-g','LineWidth',2) % 2 layer hollow spheres 50 nm diameter 5nm shell thickness
% plot(WaveLength, percentReflArray(136,:).*100,'-r','LineWidth',2) % 4 layer hollow spheres 50 nm diameter 5nm shell thickness
% title ('Hollow Spheres')
% % legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['N:\Kat FDTD Data\Figures\HollowSpheresSimple.png'])
% 

% 
% %=========================Cylinders and Thin Film================================
% figure
% hold on
% plot(WaveLength, percentReflArray(1,:).*100,'-k','LineWidth',2) % no ARC
% plot(WaveLength, percentReflArray(17,:).*100,'-b','LineWidth',2) % 1% coverage 200 nm
% plot(WaveLength, percentReflArray(18,:).*100,'--b','LineWidth',2) % 1% coverage 100 nm
% plot(WaveLength, percentReflArray(19,:).*100,':b','LineWidth',2) % 1% coverage 50 nm
% plot(WaveLength, percentReflArray(20,:).*100,'-r','LineWidth',2)% 5% coverage 200 nm
% plot(WaveLength, percentReflArray(21,:).*100,'--r','LineWidth',2)% 5% coverage 100 nm
% plot(WaveLength, percentReflArray(22,:).*100,':r','LineWidth',2)% 5% coverage 50 nm
% plot(WaveLength, percentReflArray(23,:).*100,'-c','LineWidth',2)% 20% coverage 200 nm
% plot(WaveLength, percentReflArray(24,:).*100,'--c','LineWidth',2)% 20% coverage 100 nm
% plot(WaveLength, percentReflArray(25,:).*100,':c','LineWidth',2)% 20% coverage 50 nm
% plot(WaveLength, percentReflArray(26,:).*100,'-m','LineWidth',2) % 79% coverage 200 nm
% plot(WaveLength, percentReflArray(27,:).*100,'--m','LineWidth',2) % 79% coverage 100 nm
% plot(WaveLength, percentReflArray(28,:).*100,':m','LineWidth',2) % 79% coverage 50 nm
% plot(WaveLength, percentReflArray(44,:).*100,'g','LineWidth',2)% Thin film, 800 nm, 1.25=n
% title ('Cylinders and Thin Film')
% % legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\CylindersAndThinFilm800Reflected.png'])
% 
% %=========================Cylinders 800 nm and Thin Film================================
% figure
% hold on
% plot(WaveLength, percentReflArray(1,:).*100,'-k','LineWidth',2) % no ARC
% plot(WaveLength, percentReflArray(17,:).*100,'-b','LineWidth',2) % 1% coverage 800 nm cylinder
% plot(WaveLength, percentReflArray(20,:).*100,'-r','LineWidth',2)% 5% coverage 800 nm cylinder
% plot(WaveLength, percentReflArray(23,:).*100,'-c','LineWidth',2)% 20% coverage 800 nm cylinder
% plot(WaveLength, percentReflArray(26,:).*100,'-m','LineWidth',2) % 79% coverage 800 nm cylinder
% plot(WaveLength, percentReflArray(44,:).*100,'-g','LineWidth',2)% Thin film, 800 nm, 1.25=n
% title ('Cylinders and Thin Film')
% % legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Cylinders100x00andThinFilm800Reflected.png'])


%=======================All Cones=======================================
% figure
% hold on
% plot(WaveLength, percentReflArray(1,:).*100,'-k','LineWidth',2)
% plot(WaveLength, percentReflArray(2:5,:).*100,':','LineWidth',2) % 200 nm bases
% plot(WaveLength, percentReflArray(6:9,:).*100,'--','LineWidth',2) % 100 nm bases
% plot(WaveLength, percentReflArray(10:13,:).*100,'-','LineWidth',2) % 50 nm bases
% title ('Cones')
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Cones_color_boldReflected.png'])
% 
% %=============Dense aspect ratio 1 shapes=======================================
% figure
% hold on
% 
% plot(WaveLength, percentReflArray(1,:).*100,'-k','LineWidth',2)
% plot(WaveLength, percentReflArray(51,:).*100,':b','LineWidth',2) % Checkergboard 200
% plot(WaveLength, percentReflArray(66,:).*100,'--b','LineWidth',2)% Dense Cylinders 200 nm
% plot(WaveLength, percentReflArray(29,:).*100,'-b','LineWidth',2)% Hemispheres 200
% plot(WaveLength, percentReflArray(52,:).*100,':g','LineWidth',2) % Checkergboard 100
% plot(WaveLength, percentReflArray(67,:).*100,'--g','LineWidth',2)% Dense Cylinders100 nm
% plot(WaveLength, percentReflArray(30,:).*100,'-g','LineWidth',2)% Hemispheres 100
% plot(WaveLength, percentReflArray(53,:).*100,':r','LineWidth',2) % Checkergboard 50
% plot(WaveLength, percentReflArray(68,:).*100,'--r','LineWidth',2)% Dense Cylinders 50 nm
% plot(WaveLength, percentReflArray(31,:).*100,'-r','LineWidth',2)% Hemispheres 50
% 
% title ('Dense Cubes, Hemispheres, and AR1 Cylinders')
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Dense AR1 SWSsReflected.png'])
% 
% %============Sparse aspect ratio 1 shapes=======================================
% figure
% hold on
% plot(WaveLength, percentReflArray(1,:).*100,'-k','LineWidth',2)
% plot(WaveLength, percentReflArray(14,:).*100,':b','LineWidth',2) % Cube 200
% plot(WaveLength, percentReflArray(38,:).*100,'--b','LineWidth',2)% Square cylinder sparse 200 nm
% plot(WaveLength, percentReflArray(15,:).*100,':g','LineWidth',2) % Cube 100
% plot(WaveLength, percentReflArray(39,:).*100,'--g','LineWidth',2)% Square cylinder sparse 100 nm
% plot(WaveLength, percentReflArray(69,:).*100,'-g','LineWidth',2)% Sparse Hemispheres 100
% plot(WaveLength, percentReflArray(16,:).*100,':r','LineWidth',2) % Cube 50
% plot(WaveLength, percentReflArray(40,:).*100,'--r','LineWidth',2)% Square cylinder sparse 50 nm
% plot(WaveLength, percentReflArray(70,:).*100,'-r','LineWidth',2)%  Sparse Hemispheres 50
% 
% title ('Sparse Cubes, Hemispheres, and AR1 Cylinders')
%  ylim([0 5])
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Sparse AR1 SWSsReflected.png'])
% 
% 
% %% Additional processing
% 
% %Analyzes avearage reflectance against weighted average reflectance
% AvgCompare = [avgReflArray; avgWeightedReflArray];