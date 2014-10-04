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
[~,~,RawSettings] = xlsread('DataPath.xlsx','Settings');

%Constructs arrDataPath from Excel spreadsheet
arrDataPath = RawParse(2:(size(RawParse,1)),:);

%Length of data
intDataXLength = cell2mat(RawSettings(1,2));
intDataYLength = cell2mat(RawSettings(2,2));

%Where to store the data
strAnalyzedDataPath = [RawSettings{4,2} '\'];

% Loads the reference data
load (RawSettings{3,2});

%Creates the wavelength vector
c = 299792458;
WaveLength = c./freq.*1e9;

%Loads the AM 1.5 data
load('AM1.5.mat')

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

%Loops through the sets of data to parse
for intDataParseIndex = 1:size(arrDataPath,1)

    %Prints progress of overall parsing
    fprintf('\n%i of %i datasets started.\n\n',intDataParseIndex, size(arrDataPath,1));
    
    %Loads the data
    load([strAnalyzedDataPath char(arrDataPath(intDataParseIndex,2)) '.mat'], 'txElecAvg','refElecAvg');
    
    %Builds the transmitted data array
    txArray(intDataParseIndex,:) = txElecAvg;
    
    %Builds the received data array
    refArray(intDataParseIndex,:) = refElecAvg;
    
    %Calculates the percent transmitted
    percentTxArray(intDataParseIndex,:) = txElecAvg ./ txElecAvgFS;
    
    %Calculates the percent reflected
    percentReflArray(intDataParseIndex,:) = refElecAvg ./ txElecAvgFS;
    
    %Calculates the average percent reflected 
    avgReflArray (intDataParseIndex) = sum(refElecAvg ./ txElecAvgFS) / length(txElecAvgFS);
    
    %===================================================
    % Performs weighted average
    %===================================================
    
    %Calculates the percent reflected
    PercentReflected = refElecAvg ./ txElecAvgFS;

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
clear ('arrDataPath', 'intDataParseIndex', 'strAnalyzedDataPath','refElecAvg','txElecAvg');

%Performs post processing
run AnalyzeDataReflPostProcessing
%run AnalyzeDataReflCompareToTMM