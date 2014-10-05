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


%% Plots data
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
% 
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

%=======================All Cones 2.5 nm mesh=======================================
figure
hold on
plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2)
plot(WaveLength, percentReflArray(156,:).*100,'b-','LineWidth',2) % Cones 200 nm wide 100 nm long 2.5 nm mesh
plot(WaveLength, percentReflArray(157,:).*100,'g-','LineWidth',2) % Cones 200 nm wide 200 nm long
plot(WaveLength, percentReflArray(158,:).*100,'r-','LineWidth',2) % Cones 200 nm wide 400 nm long
plot(WaveLength, percentReflArray(159,:).*100,'c-','LineWidth',2) % Cones 200 nm wide 800 nm long

plot(WaveLength, percentReflArray(160,:).*100,'b--','LineWidth',2)% Cones 100 nm wide 100 nm long
plot(WaveLength, percentReflArray(161,:).*100,'g--','LineWidth',2)% Cones 100 nm wide 200 nm long
plot(WaveLength, percentReflArray(162,:).*100,'r--','LineWidth',2)% Cones 100 nm wide 400 nm long
plot(WaveLength, percentReflArray(163,:).*100,'c--','LineWidth',2)% Cones 100 nm wide 800 nm long

plot(WaveLength, percentReflArray(164,:).*100,'b:','LineWidth',2)% Cones 50 nm wide 100 nm long
plot(WaveLength, percentReflArray(165,:).*100,'g:','LineWidth',2)% Cones 50 nm wide 200 nm long
plot(WaveLength, percentReflArray(166,:).*100,'r:','LineWidth',2)% Cones 50 nm wide 400 nm long
plot(WaveLength, percentReflArray(167,:).*100,'c:','LineWidth',2)% Cones 50 nm wide 800 nm long
title ('Cones 2.5 nm Mesh')
 ylim([0 5])

%--legend('Free Space','Transmitted','Reflected');
xlabel ('Wavelength (nm)');
ylabel ('Percent Reflected');

%Sets the font to 16, bold and saving the figure
set(findall(gcf,'-property','FontSize'),'FontSize',16) 
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Cones_color_boldReflected2_5nmMesh.png'])

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
%=============Dense aspect ratio 1 shapes 5 nm mesh=======================================
figure
hold on

plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2)
plot(WaveLength, percentReflArray(122,:).*100,':b','LineWidth',2) % Checkergboard 200
plot(WaveLength, percentReflArray(116,:).*100,'--b','LineWidth',2)% Dense Cylinders 200 nm
plot(WaveLength, percentReflArray(110,:).*100,'-b','LineWidth',2)% Hemispheres 200
plot(WaveLength, percentReflArray(123,:).*100,':g','LineWidth',2) % Checkergboard 100
plot(WaveLength, percentReflArray(117,:).*100,'--g','LineWidth',2)% Dense Cylinders100 nm
plot(WaveLength, percentReflArray(111,:).*100,'-g','LineWidth',2)% Hemispheres 100
plot(WaveLength, percentReflArray(124,:).*100,':r','LineWidth',2) % Checkergboard 50
plot(WaveLength, percentReflArray(118,:).*100,'--r','LineWidth',2)% Dense Cylinders 50 nm
plot(WaveLength, percentReflArray(112,:).*100,'-r','LineWidth',2)% Hemispheres 50

title ('Dense Cubes, Hemispheres, and AR1 Cylinders 5 nm Mesh')
 ylim([0 5])

%--legend('Free Space','Transmitted','Reflected');
xlabel ('Wavelength (nm)');
ylabel ('Percent Reflected');

%Sets the font to 16, bold and saving the figure
set(findall(gcf,'-property','FontSize'),'FontSize',16) 
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Dense AR1 SWSsReflected 5 nm mesh.png'])
% 
%============Sparse aspect ratio 1 shapes=======================================
figure
hold on
plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2)
% plot(WaveLength, percentReflArray(125,:).*100,':b','LineWidth',2) % Cube 200
% plot(WaveLength, percentReflArray(119,:).*100,'--b','LineWidth',2)% Square cylinder sparse 200 nm
plot(WaveLength, percentReflArray(126,:).*100,':g','LineWidth',2) % Cube 100
plot(WaveLength, percentReflArray(120,:).*100,'--g','LineWidth',2)% Square cylinder sparse 100 nm
plot(WaveLength, percentReflArray(114,:).*100,'-g','LineWidth',2)% Sparse Hemispheres 100
plot(WaveLength, percentReflArray(127,:).*100,':r','LineWidth',2) % Cube 50
plot(WaveLength, percentReflArray(121,:).*100,'--r','LineWidth',2)% Square cylinder sparse 50 nm
plot(WaveLength, percentReflArray(115,:).*100,'-r','LineWidth',2)%  Sparse Hemispheres 50

title ('Sparse Cubes, Hemispheres, and AR1 Cylinders 5 nm mesh')
 ylim([0 5])

%--legend('Free Space','Transmitted','Reflected');
xlabel ('Wavelength (nm)');
ylabel ('Percent Reflected');

%Sets the font to 16, bold and saving the figure
set(findall(gcf,'-property','FontSize'),'FontSize',16) 
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Sparse AR1 SWSsReflected 5 nm mesh.png'])
% 
% 
% %% Additional processing
% 
% %Analyzes avearage reflectance against weighted average reflectance
% AvgCompare = [avgReflArray; avgWeightedReflArray];












%Downsampling divisor
DownSampleTMM = 9;

%Imports the data
[numTMMCompiled,txtTMMCompiled,rawTMMCompiled]=xlsread('Compiled TMM Results.xlsx');




%% Plots the results - FIGURE 1
%Note to self:  The next 40-ish lines will plot the same FDTD lines that I
%would normally make in the AnalyzeDataReflExcel.m file, but there is code
%tagged on the end that will plot TMM data from Luxpop on top of those
%lines in the form of points.  The lines come from data entered in to the
%Excel file called "DataPath.xlsx".  The number in (51,:), etc., refers to
%which line of data in DataPath we would like to plot.  

%=========================Hollow Spheres================================
figure
hold on
plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2) % no ARC
plot(WaveLength, percentReflArray(128,:).*100,'-b','LineWidth',2) % 1 layer hollow spheres 150 nm diameter 5nm shell thickness
plot(WaveLength, percentReflArray(129,:).*100,'-g','LineWidth',2) % 1 layer hollow spheres 100 nm diameter 5nm shell thickness
plot(WaveLength, percentReflArray(130,:).*100,'-r','LineWidth',2) % 1 layer hollow spheres 50 nm diameter 5nm shell thickness
plot(WaveLength, percentReflArray(131,:).*100,'--b','LineWidth',2)% 2 layer2 hollow spheres 150 nm diameter 5nm shell thickness
plot(WaveLength, percentReflArray(132,:).*100,'--g','LineWidth',2)% 2 layer2 hollow spheres 100 nm diameter 5nm shell thickness
plot(WaveLength, percentReflArray(133,:).*100,'--r','LineWidth',2)% 2 layer2 hollow spheres 50 nm diameter 5nm shell thickness
plot(WaveLength, percentReflArray(134,:).*100,'-.b','LineWidth',2)% 4 layer2 hollow spheres 150 nm diameter 5nm shell thickness
plot(WaveLength, percentReflArray(135,:).*100,'-.g','LineWidth',2)% 4 layer2 hollow spheres 100 nm diameter 5nm shell thickness
plot(WaveLength, percentReflArray(136,:).*100,'-.r','LineWidth',2)% 4 layer2 hollow spheres 150 nm diameter 5nm shell thickness
% plot(WaveLength, percentReflArray(137,:).*100,':b','LineWidth',2) % 8 layer2 hollow spheres 150 nm diameter 5nm shell thickness
% plot(WaveLength, percentReflArray(138,:).*100,':g','LineWidth',2) % 8 layer2 hollow spheres 100 nm diameter 5nm shell thickness
% plot(WaveLength, percentReflArray(139,:).*100,':r','LineWidth',2) % 8 layer2 hollow spheres 50 nm diameter 5nm shell thickness
title ('Hollow Spheres')
% legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
 ylim([0 5])

%--legend('Free Space','Transmitted','Reflected');
xlabel ('Wavelength (nm)');
ylabel ('Percent Reflected');

%Sets the font to 16, bold and saving the figure
set(findall(gcf,'-property','FontSize'),'FontSize',16) 
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
print(gcf,'-dpng','-r600',['N:\Kat FDTD Data\Figures\HollowSpheres.png'])

% %=========================375 nm Cylinders 5 nm mesh and thin film================================
figure
hold on
plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2) % no ARC

plot(WaveLength, percentReflArray(53,:).*100,'-g','LineWidth',2) % Cylinder full width 200 nm wide
plot(WaveLength, percentReflArray(54,:).*100,'--g','LineWidth',2) % Cylinder full width 100 nm wide
plot(WaveLength, percentReflArray(55,:).*100,':g','LineWidth',2)% Cylinder full width 50 nm wide
plot(WaveLength, percentReflArray(56,:).*100,'-r','LineWidth',2)% Cylinder 80% width 200 nm wide
plot(WaveLength, percentReflArray(57,:).*100,'--r','LineWidth',2)% Cylinder 80% width 100 nm wide
plot(WaveLength, percentReflArray(58,:).*100,':r','LineWidth',2) % Cylinder 80% width 50 nm wide
% plot(WaveLength, percentReflArray(59,:).*100,'-b','LineWidth',2) % Cylinder half width 200 nm wide
plot(WaveLength, percentReflArray(60,:).*100,'--b','LineWidth',2) % Cylinder half width 100 nm wide
plot(WaveLength, percentReflArray(61,:).*100,':b','LineWidth',2) % Cylinder half width 50 nm wide
plot(WaveLength, percentReflArray(52,:).*100,'-c','LineWidth',2) % Thin film 375 nm n=1.22

%Now we have a nice plot with 972 points that has one x vector with many y
%vectors plotted against it.  We also have about a hundred points over the
%same wavelength ranges that have TMM-calculated reflectance values from
%Luxpop.  These values can be found in the excel file "Compiled TMM
%Results.xlsx", which should be in the same folder as all these other
%MATLAB files.  We need to get this data onto the same axis as the lines
%with the right density of points to see the data without obscuring the
%rest of the lines.  To do this we will interpolate between the points and
%plot the interpolated values at a spacing of "DownSampleTMM".  
%=========================TMM================================

%Removes the last row of data, which is invalid
numTMMCompiled = numTMMCompiled(1:length(numTMMCompiled) - 1,:);

%Interpolates the TMM data
numTMMData(:,12) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,12)',WaveLength);
numTMMData(:,13) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,13)',WaveLength);
numTMMData(:,14) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,14)',WaveLength);

%Sets the first and last rows of the data to the data point next to it.
numTMMData(1,:) = numTMMData(2,:);
numTMMData(length(numTMMData),:) = numTMMData(length(numTMMData)-1,:);

%Gets the length of the wavelength
wvLength = size(WaveLength,2);

%Plots the TMM interpolated points on the same plot as the FDTD lines
plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),12).*100,'go','LineWidth',2);
plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),13).*100,'ro','LineWidth',2);
plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),14).*100,'bo','LineWidth',2);

title ('Cylinders 375 nm long and thin film, n=1.22, 375 nm thick, FDTD and TMM')
% legend(char([strLegend(1) strLegend(17:28)]),'Location','Best');
 ylim([0 5])

%--legend('Free Space','Transmitted','Reflected');
xlabel ('Wavelength (nm)');
ylabel ('Percent Reflected');

%Sets the font to 16, bold and saving the figure
set(findall(gcf,'-property','FontSize'),'FontSize',16) 
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Cylinders375nmFDTDandTMM.png'])

%================== Pyramids FDTD 2.5 nm mesh and TMM 1 nm slices=======================================
figure
hold on
plot(WaveLength, percentReflArray(51,:).*100,'-k','LineWidth',2)
plot(WaveLength, percentReflArray(144,:).*100,'b-','LineWidth',2) % Pyramids 200 nm wide 100 nm long 2.5 nm mesh
plot(WaveLength, percentReflArray(145,:).*100,'g-','LineWidth',2) % Pyramids 200 nm wide 200 nm long
plot(WaveLength, percentReflArray(146,:).*100,'r-','LineWidth',2) % Pyramids 200 nm wide 400 nm long
plot(WaveLength, percentReflArray(147,:).*100,'c-','LineWidth',2) % Pyramids 200 nm wide 800 nm long

plot(WaveLength, percentReflArray(148,:).*100,'b--','LineWidth',2)% Pyramids 100 nm wide 100 nm long
plot(WaveLength, percentReflArray(149,:).*100,'g--','LineWidth',2)% Pyramids 100 nm wide 200 nm long
plot(WaveLength, percentReflArray(150,:).*100,'r--','LineWidth',2)% Pyramids 100 nm wide 400 nm long
plot(WaveLength, percentReflArray(151,:).*100,'c--','LineWidth',2)% Pyramids 100 nm wide 800 nm long

plot(WaveLength, percentReflArray(152,:).*100,'b:','LineWidth',2)% Pyramids 50 nm wide 100 nm long
plot(WaveLength, percentReflArray(153,:).*100,'g:','LineWidth',2)% Pyramids 50 nm wide 200 nm long
plot(WaveLength, percentReflArray(154,:).*100,'r:','LineWidth',2)% Pyramids 50 nm wide 400 nm long
plot(WaveLength, percentReflArray(155,:).*100,'c:','LineWidth',2)% Pyramids 50 nm wide 800 nm long

%=========================TMM================================

%Removes the last row of data, which is invalid
numTMMCompiled = numTMMCompiled(1:length(numTMMCompiled) - 1,:);

%Interpolates the TMM data
numTMMData(:,8) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,8)',WaveLength);
numTMMData(:,9) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,9)',WaveLength);
numTMMData(:,10) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,10)',WaveLength);
numTMMData(:,11) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,11)',WaveLength);

%Sets the first and last rows of the data to the data point next to it.
numTMMData(1,:) = numTMMData(2,:);
numTMMData(length(numTMMData),:) = numTMMData(length(numTMMData)-1,:);

%Gets the length of the wavelength
wvLength = size(WaveLength,2);

%Plots the TMM interpolated points on the same plot as the FDTD lines
plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),8).*100,'bo','LineWidth',2);
plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),9).*100,'go','LineWidth',2);
plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),10).*100,'ro','LineWidth',2);
plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),11).*100,'co','LineWidth',2);


title ('Pyramids FDTD and TMM 5 nm resolution')
 ylim([0 5])

%--legend('Free Space','Transmitted','Reflected');
xlabel ('Wavelength (nm)');
ylabel ('Percent Reflected');

%Sets the font to 16, bold and saving the figure
set(findall(gcf,'-property','FontSize'),'FontSize',16) 
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\Pyramids2_5nmMeshReflectedFDTDandTMM.png'])

% %================== Pyramids 5 nm mesh=======================================
% 
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
% %=========================TMM================================
% 
% %Removes the last row of data, which is invalid
% numTMMCompiled = numTMMCompiled(1:length(numTMMCompiled) - 1,:);
% 
% %Interpolates the TMM data
% numTMMData(:,8) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,8)',WaveLength);
% numTMMData(:,9) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,9)',WaveLength);
% numTMMData(:,10) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,10)',WaveLength);
% numTMMData(:,11) = interp1(numTMMCompiled(:,1)',numTMMCompiled(:,11)',WaveLength);
% 
% %Sets the first and last rows of the data to the data point next to it.
% numTMMData(1,:) = numTMMData(2,:);
% numTMMData(length(numTMMData),:) = numTMMData(length(numTMMData)-1,:);
% 
% %Gets the length of the wavelength
% wvLength = size(WaveLength,2);
% 
% %Plots the TMM interpolated points on the same plot as the FDTD lines
% plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),8).*100,'bo','LineWidth',2);
% plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),9).*100,'go','LineWidth',2);
% plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),10).*100,'ro','LineWidth',2);
% plot(WaveLength(1:DownSampleTMM:wvLength),numTMMData((1:DownSampleTMM:wvLength),11).*100,'co','LineWidth',2);
% 
% 
% title ('Pyramids FDTD 5 nm mesh and TMM 1 nm resolution')
%  ylim([0 5])
% 
% 
% %--legend('Free Space','Transmitted','Reflected');
% xlabel ('Wavelength (nm)');
% ylabel ('Percent Reflected');
% 
% %Sets the font to 16, bold and saving the figure
% set(findall(gcf,'-property','FontSize'),'FontSize',16) 
% set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
% print(gcf,'-dpng','-r600',['M:\Kat FDTD Data\Figures\PyramidsWideTallSweepsReflected5nmMeshTMM1nmres.png'])
% 
