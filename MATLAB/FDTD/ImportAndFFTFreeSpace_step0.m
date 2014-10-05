%% Initializes the script

%Clears the screen
clc;

%Clears all the variables
clear all;

%Closes all windows
close all;


%% File Path Info 

%Reads the Excel spreadsheet for parse locations and file settings
[~,~,RawParse] = xlsread('DataPath.xlsx','Free Space Raw FDTD Data');

%Constructs arrDataPath from Excel spreadsheet
arrDataPath = RawParse(2:(size(RawParse,1)),:);

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


%% Loads the reflected light data from the hard drive

%Prints starting message
fprintf('\n\nGetting started!\n\nMake sure to close all other unnecessary programs...this will take some memory.\nClick on the command window and press Ctrl+C to terminate the job.\n');

%Loops through the sets of data to parse
for intDataParseIndex = 1:intDataParseSize

    %Prints progress of overall parsing
    fprintf('\n%i of %i datasets started.\n\n',intDataParseIndex, intDataParseSize);
    
    %% Parses XFTDT data

    %Clears the data
    clear reflected*
    clear transmitted*

    %Length of data
    intDataXLength = cell2mat(arrDataPath(intDataParseIndex,9));
    intDataYLength = cell2mat(arrDataPath(intDataParseIndex,10));

    %Where to store the data
    strAnalyzedDataPath = [char(arrDataPath(intDataParseIndex,11)) '\'];

    %"Reflected" file name stub
    strFileReflected = ['\' char(arrDataPath(intDataParseIndex,12))];

    %"Transmitted" file name stub
    strFileTransmitted = ['\' char(arrDataPath(intDataParseIndex,13))];
    
    %Loops using Parallel FOR to load reflected light
    for intIndex = 0:(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)
        %Sets up the file paths
        dataFileExr =  [char(arrDataPath(intDataParseIndex,1)) strFileReflected 'Ex_total_ts' int2str(intIndex) '.bin'];
        dataFileExi = '';
        dataFileEyr =  [ char(arrDataPath(intDataParseIndex,1)) strFileReflected 'Ey_total_ts' int2str(intIndex) '.bin'];
        dataFileEyi = '';
        dataFileHxr =  [char(arrDataPath(intDataParseIndex,1)) strFileReflected 'Hx_total_ts' int2str(intIndex) '.bin'];
        dataFileHxi = '';
        dataFileHyr =  [ char(arrDataPath(intDataParseIndex,1)) strFileReflected 'Hy_total_ts' int2str(intIndex) '.bin'];
        dataFileHyi = '';

        %Retrieves the data from the BIN file
        [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData(dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength,intDataYLength);

        %Concatinates the data matricies
        reflectedExr(:,:,intIndex+1) = dataExr;
        reflectedEyr(:,:,intIndex+1) = dataEyr;
        reflectedHxr(:,:,intIndex+1) = dataHxr;
        reflectedHyr(:,:,intIndex+1) = dataHyr;
        reflectedExi(:,:,intIndex+1) = dataExi;
        reflectedEyi(:,:,intIndex+1) = dataEyi;
        reflectedHxi(:,:,intIndex+1) = dataHxi;
        reflectedHyi(:,:,intIndex+1) = dataHyi;

        %Displays the progress
        if mod(intIndex,30) == 0 
            ReflectedLightPercentDone = round(10000*(intIndex/(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)))/100;
            fprintf('Reflected Percent Done: %3.1f%%\n',ReflectedLightPercentDone);
        end
    end



    %% Loads the transmitted light data from the hard drive

    %Prints half done
    fprintf('\n\nHalf done with current run of data!\n');

    %Loops using Parallel FOR to load transmitted light
    for intIndex = 0:(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)
        %Sets up the file paths
        dataFileExr =  [char(arrDataPath(intDataParseIndex,1)) strFileTransmitted 'Ex_total_ts' int2str(intIndex) '.bin'];
        dataFileExi = '';
        dataFileEyr =  [char(arrDataPath(intDataParseIndex,1)) strFileTransmitted 'Ey_total_ts' int2str(intIndex) '.bin'];
        dataFileEyi = '';
        dataFileHxr =  [char(arrDataPath(intDataParseIndex,1)) strFileTransmitted 'Hx_total_ts' int2str(intIndex) '.bin'];
        dataFileHxi = '';
        dataFileHyr =  [char(arrDataPath(intDataParseIndex,1)) strFileTransmitted 'Hy_total_ts' int2str(intIndex) '.bin'];
        dataFileHyi = '';

        %Retrieves the data from the BIN file
        [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData(dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength,intDataYLength);

        %Concatinates the data matricies
        transmittedExr(:,:,intIndex+1) = dataExr;
        transmittedEyr(:,:,intIndex+1) = dataEyr;
        transmittedHxr(:,:,intIndex+1) = dataHxr;
        transmittedHyr(:,:,intIndex+1) = dataHyr;
        transmittedExi(:,:,intIndex+1) = dataExi;
        transmittedEyi(:,:,intIndex+1) = dataEyi;
        transmittedHxi(:,:,intIndex+1) = dataHxi;
        transmittedHyi(:,:,intIndex+1) = dataHyi;

        %Displays the progress
        if mod(intIndex,30) == 0 
            TransmittedLightPercentDone = round(10000*(intIndex/(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)))/100;
            fprintf('Transmitted Percent Done: %3.1f%%\n',TransmittedLightPercentDone);
        end 
    end

    %Performs FFT on data
    run PerformFFTFreeSpace
    
    %Creates freq
    freq = f(freqStartIndex:freqStopIndex);
    
    
    %=========================================================================
    % Average and Compare Data
    %=========================================================================
    %% Performs averaging

    %Creates the averaged frequency vectors
    txAvgFS = zeros(1,size(tx_FFT,3));
    refAvgFS = zeros(1,size(ref_FFT,3));


    %Loops through the entire data block, averaging every frequency point, for
    %transmitted light
    for intIndex = 1:size(tx_FFT,3)
        %Averages a data slice of frequency
        txAvgFS(intIndex) = sum(sum(tx_FFT(:,:,intIndex)))/(size(tx_FFT,1)*size(tx_FFT,2));
    end


    %Loops through the entire data block, averaging every frequency point, for
    %reflected light
    for intIndex = 1:size(ref_FFT,3)
        %Averages a data slice of frequency
        refAvgFS(intIndex) = sum(sum(ref_FFT(:,:,intIndex)))/(size(ref_FFT,1)*size(ref_FFT,2));
    end

    %Accounts for index of refraction on transmitted light
    txAvgFS =  txAvgFS .* cell2mat((arrDataPath(intDataParseIndex,3)));

 
        
    %Saves datafile
    fprintf('\n\nSaving your data - this will also take some time...\n');
    save([strAnalyzedDataPath char(arrDataPath(intDataParseIndex,2)) '.mat'], 'freq','txAvgFS','refAvgFS')
    
    %Prints half done
    fprintf('\n\nAll done with one run of data!\n'); 
end


