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

%Sorts the data path array by whether we're parsing or not
arrDataPath = sortrows(arrDataPath,6);

%Finds the first instance where we start parsing
for intIndex= 1:size(arrDataPath,1)
    %Checks to find the first instance of 1
    if isequal(arrDataPath(intIndex,6), {1})
        %Figures out the parse size
        intDataParseSize = size(arrDataPath,1) - intIndex + 1;
        
        %Truncates the array
        arrDataPath = arrDataPath(intIndex:size(arrDataPath,1),:);
        
        %Breaks out of the search loop
        break;
    end
end


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

%Performs post processing
run AnalyzeDataReflPostProcessing
%run AnalyzeDataReflCompareToTMM