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

%Loops through the sets of data to parse
for intDataParseIndex = 1:intDataParseSize

    %Prints progress of overall parsing
    fprintf('\n%i of %i datasets started.\n\n',intDataParseIndex, intDataParseSize);
    
    %Length of data
    intDataXLength = cell2mat(arrDataPath(intDataParseIndex,9));
    intDataYLength = cell2mat(arrDataPath(intDataParseIndex,10));

    %Where to store the data
    strAnalyzedDataPath = [char(arrDataPath(intDataParseIndex,12)) '\'];

    %"Reflected" file name stub
    strFileReflected = ['\' char(arrDataPath(intDataParseIndex,13))];

    %"Transmitted" file name stub
    strFileTransmitted = ['\' char(arrDataPath(intDataParseIndex,14))];

    %% Parses Reflected Ex 
    
    %Preallocates the reflected data
    reflectedExr = zeros(intDataXLength,intDataYLength,cell2mat(arrDataPath(intDataParseIndex,4)));

    %Loops to load reflected light Ex
    for intIndex = 0:(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)
        %Sets up the file paths
        dataFileExr =  [char(arrDataPath(intDataParseIndex,1)) strFileReflected 'Ex_total_ts' int2str(intIndex) '.bin'];
        dataFileExi = '';
        dataFileEyr = '';
        dataFileEyi = '';
        dataFileHxr = '';
        dataFileHxi = '';
        dataFileHyr = '';
        dataFileHyi = '';

        %Retrieves the data from the BIN file
        [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData(dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength,intDataYLength);

        %Concatinates the data matricies
        reflectedExr(:,:,intIndex+1) = dataExr;

        %Displays the progress
        if mod(intIndex,30) == 0 
            ReflectedLightPercentDone = round(10000*(intIndex/(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)))/100;
            fprintf('%i of %i datasets - reflected Ex Percent Done: %3.1f%%\n',intDataParseIndex, intDataParseSize,ReflectedLightPercentDone); 
        end
    end

    %Performs FFT on data
    fprintf('\n\nPerforming Reflected Ex FFT...\n\n');
    refEx_FFT = PerformFFT(reflectedExr,cell2mat(arrDataPath(intDataParseIndex,5)),(cell2mat(arrDataPath(intDataParseIndex,7))*1e-9),(cell2mat(arrDataPath(intDataParseIndex,8))*1e-9));
    
    %Clears unused variables
    clear reflectedExr;

     %% Parses Reflected Ey
    
    %Preallocates the reflected data
    reflectedEyr = zeros(intDataXLength,intDataYLength,cell2mat(arrDataPath(intDataParseIndex,4)));

    %Loops to load reflected light Ey
    for intIndex = 0:(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)
        %Sets up the file paths
        dataFileExr = '';
        dataFileExi = '';
        dataFileEyr =  [ char(arrDataPath(intDataParseIndex,1)) strFileReflected 'Ey_total_ts' int2str(intIndex) '.bin'];
        dataFileEyi = '';
        dataFileHxr = '';
        dataFileHxi = '';
        dataFileHyr = '';
        dataFileHyi = '';

        %Retrieves the data from the BIN file
        [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData(dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength,intDataYLength);

        %Concatinates the data matricies
        reflectedEyr(:,:,intIndex+1) = dataEyr;

        %Displays the progress
        if mod(intIndex,30) == 0 
            ReflectedLightPercentDone = round(10000*(intIndex/(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)))/100;
            fprintf('%i of %i datasets - reflected Ey Percent Done: %3.1f%%\n',intDataParseIndex, intDataParseSize,ReflectedLightPercentDone); 
        end
    end

    %Performs FFT on data
    fprintf('\n\nPerforming Reflected Ey FFT...\n\n');
    refEy_FFT = PerformFFT(reflectedEyr,cell2mat(arrDataPath(intDataParseIndex,5)),(cell2mat(arrDataPath(intDataParseIndex,7))*1e-9),(cell2mat(arrDataPath(intDataParseIndex,8))*1e-9));
    
    %Clears unused variables
    clear reflectedEyr;
    
    %% Parses Reflected Hx 
    
    %Preallocates the reflected data
    reflectedHxr = zeros(intDataXLength,intDataYLength,cell2mat(arrDataPath(intDataParseIndex,4)));

    %Loops to load reflected light Ex
    for intIndex = 0:(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)
        %Sets up the file paths
        dataFileExr = '';
        dataFileExi = '';
        dataFileEyr = '';
        dataFileEyi = '';
        dataFileHxr = [char(arrDataPath(intDataParseIndex,1)) strFileReflected 'Hx_total_ts' int2str(intIndex) '.bin'];
        dataFileHxi = '';
        dataFileHyr = '';
        dataFileHyi = '';

        %Retrieves the data from the BIN file
        [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData(dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength,intDataYLength);

        %Concatinates the data matricies
        reflectedHxr(:,:,intIndex+1) = dataHxr;

        %Displays the progress
        if mod(intIndex,30) == 0 
            ReflectedLightPercentDone = round(10000*(intIndex/(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)))/100;
            fprintf('%i of %i datasets - reflected Hx Percent Done: %3.1f%%\n',intDataParseIndex, intDataParseSize,ReflectedLightPercentDone); 
        end
    end
   
    %Performs FFT on data
    fprintf('\n\nPerforming Reflected Hx FFT...\n\n');
    refHx_FFT = PerformFFT(reflectedHxr,cell2mat(arrDataPath(intDataParseIndex,5)),(cell2mat(arrDataPath(intDataParseIndex,7))*1e-9),(cell2mat(arrDataPath(intDataParseIndex,8))*1e-9));
    
    %Clears unused variables
    clear reflectedHxr;
    
    %% Parses Reflected Hy
    
    %Preallocates the reflected data
    reflectedHyr = zeros(intDataXLength,intDataYLength,cell2mat(arrDataPath(intDataParseIndex,4)));

    %Loops to load reflected light Ey
    for intIndex = 0:(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)
        %Sets up the file paths
        dataFileExr = '';
        dataFileExi = '';
        dataFileEyr = ''; 
        dataFileEyi = '';
        dataFileHxr = '';
        dataFileHxi = '';
        dataFileHyr = [ char(arrDataPath(intDataParseIndex,1)) strFileReflected 'Hy_total_ts' int2str(intIndex) '.bin'];
        dataFileHyi = '';

        %Retrieves the data from the BIN file
        [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData(dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength,intDataYLength);

        %Concatinates the data matricies
        reflectedHyr(:,:,intIndex+1) = dataHyr;

        %Displays the progress
        if mod(intIndex,30) == 0 
            ReflectedLightPercentDone = round(10000*(intIndex/(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)))/100;
            fprintf('%i of %i datasets - reflected Hy Percent Done: %3.1f%%\n',intDataParseIndex, intDataParseSize,ReflectedLightPercentDone); 
        end
    end
    
    %Performs FFT on data
    fprintf('\n\nPerforming Reflected Hx FFT...\n\n');
    refHy_FFT = PerformFFT(reflectedHyr,cell2mat(arrDataPath(intDataParseIndex,5)),(cell2mat(arrDataPath(intDataParseIndex,7))*1e-9),(cell2mat(arrDataPath(intDataParseIndex,8))*1e-9));
    
    %Clears unused variables
    clear reflectedHyr;
    
    %% Averages Reflected Data
    
    %Performs addition on the light intensity on each electric field direction
    ref_FFT = real(1/2.* refEx_FFT .* conj(refHy_FFT)) - real(1/2.*refEy_FFT .* conj(refHx_FFT));

    %Deletes unused variables
    clear refEx_FFT refEy_FFT refHx_FFT refHy_FFT;
    
    %Creates the averaged frequency vectors
    refAvg = zeros(1,size(ref_FFT,3));

    %Loops through the entire data block, averaging every frequency point, for
    %reflected light
    for intIndex = 1:size(ref_FFT,3)
        %Averages a data slice of frequency
        refAvg(intIndex) = sum(sum(ref_FFT(:,:,intIndex)))/(size(ref_FFT,1)*size(ref_FFT,2));
    end
   
    %% Parses Transmitted Ex
   
    % Preallocates the transmitted data
    transmittedExr = zeros(intDataXLength,intDataYLength,cell2mat(arrDataPath(intDataParseIndex,4)));

    %Loops using Parallel FOR to load transmitted light
    for intIndex = 0:(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)
        %Sets up the file paths
        dataFileExr =  [char(arrDataPath(intDataParseIndex,1)) strFileTransmitted 'Ex_total_ts' int2str(intIndex) '.bin'];
        dataFileExi = '';
        dataFileEyr = '';
        dataFileEyi = '';
        dataFileHxr = '';
        dataFileHxi = '';
        dataFileHyr = '';
        dataFileHyi = '';

        %Retrieves the data from the BIN file
        [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData( dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength,intDataYLength);

        %Concatinates the data matricies
        transmittedExr(:,:,intIndex+1) = dataExr;
        
         %Displays the progress
        if mod(intIndex,30) == 0 
            TransmittedLightPercentDone = round(10000*(intIndex/(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)))/100;
            fprintf('%i of %i datasets - transmitted Ex Percent Done: %3.1f%%\n',intDataParseIndex, intDataParseSize,TransmittedLightPercentDone); 
        end 
    end
    
    %Performs FFT on data
    fprintf('\n\nPerforming Transmitted Ex FFT...\n\n');
    txEx_FFT = PerformFFT(transmittedExr,cell2mat(arrDataPath(intDataParseIndex,5)),(cell2mat(arrDataPath(intDataParseIndex,7))*1e-9),(cell2mat(arrDataPath(intDataParseIndex,8))*1e-9));
    
    %Clears unused variables
    clear transmittedExr;
    
    
    %% Parses Transmitted Ey
   
    % Preallocates the transmitted data
    transmittedEyr = zeros(intDataXLength,intDataYLength,cell2mat(arrDataPath(intDataParseIndex,4)));

    %Loops using Parallel FOR to load transmitted light
    for intIndex = 0:(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)
        %Sets up the file paths
        dataFileExr = '';
        dataFileExi = '';
        dataFileEyr =  [char(arrDataPath(intDataParseIndex,1)) strFileTransmitted 'Ey_total_ts' int2str(intIndex) '.bin'];
        dataFileEyi = '';
        dataFileHxr = '';
        dataFileHxi = '';
        dataFileHyr = '';
        dataFileHyi = '';

        %Retrieves the data from the BIN file
        [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData( dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength,intDataYLength);

        %Concatinates the data matricies
        transmittedEyr(:,:,intIndex+1) = dataEyr;
        
         %Displays the progress
        if mod(intIndex,30) == 0 
            TransmittedLightPercentDone = round(10000*(intIndex/(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)))/100;
            fprintf('%i of %i datasets - transmitted Ey Percent Done: %3.1f%%\n',intDataParseIndex, intDataParseSize,TransmittedLightPercentDone); 
        end 
    end

    %Performs FFT on data
    fprintf('\n\nPerforming Transmitted Ey FFT...\n\n');
    txEy_FFT = PerformFFT(transmittedEyr,cell2mat(arrDataPath(intDataParseIndex,5)),(cell2mat(arrDataPath(intDataParseIndex,7))*1e-9),(cell2mat(arrDataPath(intDataParseIndex,8))*1e-9));
    
    %Clears unused variables
    clear transmittedEyr;
    
    %% Parses Transmitted Hx
   
    % Preallocates the transmitted data
    transmittedHxr = zeros(intDataXLength,intDataYLength,cell2mat(arrDataPath(intDataParseIndex,4)));

    %Loops using Parallel FOR to load transmitted light
    for intIndex = 0:(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)
        %Sets up the file paths
        dataFileExr = ''; 
        dataFileExi = '';
        dataFileEyr = '';
        dataFileEyi = '';
        dataFileHxr = [char(arrDataPath(intDataParseIndex,1)) strFileTransmitted 'Hx_total_ts' int2str(intIndex) '.bin'];
        dataFileHxi = '';
        dataFileHyr = '';
        dataFileHyi = '';

        %Retrieves the data from the BIN file
        [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData( dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength,intDataYLength);

        %Concatinates the data matricies
        transmittedHxr(:,:,intIndex+1) = dataHxr;
        
         %Displays the progress
        if mod(intIndex,30) == 0 
            TransmittedLightPercentDone = round(10000*(intIndex/(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)))/100;
            fprintf('%i of %i datasets - transmitted Hx Percent Done: %3.1f%%\n',intDataParseIndex, intDataParseSize,TransmittedLightPercentDone); 
        end 
    end

    %Performs FFT on data
    fprintf('\n\nPerforming Transmitted Hx FFT...\n\n');
    txHx_FFT = PerformFFT(transmittedHxr,cell2mat(arrDataPath(intDataParseIndex,5)),(cell2mat(arrDataPath(intDataParseIndex,7))*1e-9),(cell2mat(arrDataPath(intDataParseIndex,8))*1e-9));
    
    %Clears unused variables
    clear transmittedHxr;
    
    
    %% Parses Transmitted Hy
   
    % Preallocates the transmitted data
    transmittedHyr = zeros(intDataXLength,intDataYLength,cell2mat(arrDataPath(intDataParseIndex,4)));

    %Loops using Parallel FOR to load transmitted light
    for intIndex = 0:(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)
        %Sets up the file paths
        dataFileExr = '';
        dataFileExi = '';
        dataFileEyr = ''; 
        dataFileEyi = '';
        dataFileHxr = '';
        dataFileHxi = '';
        dataFileHyr = [char(arrDataPath(intDataParseIndex,1)) strFileTransmitted 'Hy_total_ts' int2str(intIndex) '.bin'];
        dataFileHyi = '';

        %Retrieves the data from the BIN file
        [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData( dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength,intDataYLength);

        %Concatinates the data matricies
        transmittedHyr(:,:,intIndex+1) = dataHyr;
        
         %Displays the progress
        if mod(intIndex,30) == 0 
            TransmittedLightPercentDone = round(10000*(intIndex/(cell2mat(arrDataPath(intDataParseIndex,4)) - 1)))/100;
            fprintf('%i of %i datasets - transmitted Hy Percent Done: %3.1f%%\n',intDataParseIndex, intDataParseSize,TransmittedLightPercentDone); 
        end 
    end

    %Performs FFT on data
    fprintf('\n\nPerforming Transmitted Hy FFT...\n\n');
    txHy_FFT = PerformFFT(transmittedHyr,cell2mat(arrDataPath(intDataParseIndex,5)),(cell2mat(arrDataPath(intDataParseIndex,7))*1e-9),(cell2mat(arrDataPath(intDataParseIndex,8))*1e-9));
    
    %Clears unused variables
    clear transmittedHyr;
    
    %% Averages transmitted data
    
    % Performs addition on the light intensity on each electric field direction
    tx_FFT =  real(1/2.*txEy_FFT .* conj(txHx_FFT)) - real(1/2.* txEx_FFT .* conj(txHy_FFT));
    
    %Deletes unused variables
    clear txEx_FFT txEy_FFT txHx_FFT txHy_FFT;
 
    %Creates the averaged frequency vectors
    txAvg = zeros(1,size(tx_FFT,3));

    %Loops through the entire data block, averaging every frequency point, for
    %transmitted light
    for intIndex = 1:size(tx_FFT,3)
        %Averages a data slice of frequency
        txAvg(intIndex) = sum(sum(tx_FFT(:,:,intIndex)))/(size(tx_FFT,1)*size(tx_FFT,2));
    end
    
    %Reads permittivity
    n =  cell2mat((arrDataPath(intDataParseIndex,3)));
    
    %% Runs post processing
    run LoadAndAnalyzeDataPostProcessing
        
    %% Saves datafile
    save([strAnalyzedDataPath char(arrDataPath(intDataParseIndex,2)) '.mat'], 'freq','txAvg','refAvg*' ,'n')
end


