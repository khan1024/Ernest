%==================================
%Performs FFT
%==================================

%Saves datafile
fprintf('\n\nPerforms FFT - this will also take some time...\n');

%Defines the speed of light
c = 299792458;

%Sets the frequencies we're interested in
freqStart = c/(cell2mat(arrDataPath(intDataParseIndex,8))*1e-9); % frequency in Hz
freqStop = c/(cell2mat(arrDataPath(intDataParseIndex,7))*1e-9);

%Sets the matrix to transform
matTransform = transmittedExr;

%Figures out the length of the sample
L = size(matTransform,3);

%Time step
dT = cell2mat(arrDataPath(intDataParseIndex,5));


%Sampling Rate
fs = 1/dT;

 %Next power of 2 from length of y, or power of 2 that is larger than y
NFFT = 2^16;

%Calculates the frequency axis
f = fs/2*linspace(0,1,NFFT/2+1);


%Calculates the starting and stopping indicies for the frequencies we're
%interested in
freqStartIndex = find(f>freqStart);
freqStartIndex = freqStartIndex(1);
freqStopIndex = find(f <=freqStop);
freqStopIndex = freqStopIndex(length(freqStopIndex));

%Initializes the 3D FFT matrix
txEx_FFT = zeros(size(matTransform,1),size(matTransform,2),(freqStopIndex-freqStartIndex+1));
txEy_FFT = zeros(size(matTransform,1),size(matTransform,2),(freqStopIndex-freqStartIndex+1));
txHx_FFT = zeros(size(matTransform,1),size(matTransform,2),(freqStopIndex-freqStartIndex+1));
txHy_FFT = zeros(size(matTransform,1),size(matTransform,2),(freqStopIndex-freqStartIndex+1));

refEx_FFT = zeros(size(matTransform,1),size(matTransform,2),(freqStopIndex-freqStartIndex+1));
refEy_FFT = zeros(size(matTransform,1),size(matTransform,2),(freqStopIndex-freqStartIndex+1));
refHx_FFT = zeros(size(matTransform,1),size(matTransform,2),(freqStopIndex-freqStartIndex+1));
refHy_FFT = zeros(size(matTransform,1),size(matTransform,2),(freqStopIndex-freqStartIndex+1));


%% ======================================
%PERFORMS FFT: transmittedExr
%========================================


%Performs 2D FFT
for intIndex = 1:size(transmittedExr,1)
    %Retrieves a slice of data from the matrix
    DataSlice = reshape(transmittedExr(intIndex,:,:),[size(transmittedExr,2),size(transmittedExr,3)]);

    %Performs FFT on the data slice
    DataSliceFFT= fft(DataSlice,NFFT,2);
    
    %Chops down the FFT to the right size?
    DataSliceFFT = DataSliceFFT(:,1:NFFT/2+1);

    %Trims down the data slice to just the data we're interested in and
    %stores the data into the 3D FFT matrix
    txEx_FFT(intIndex,:,:) = DataSliceFFT(:,freqStartIndex:freqStopIndex);
    
    %Removes the data slice FFT
    clear DataSliceFFT;

end 


%% ======================================
%PERFORMS FFT: transmittedEyr
%========================================

%Performs 2D FFT
for intIndex = 1:size(transmittedEyr,1)
    %Retrieves a slice of data from the matrix
    DataSlice = reshape(transmittedEyr(intIndex,:,:),[size(transmittedEyr,2),size(transmittedEyr,3)]);

    %Performs FFT on the data slice
    DataSliceFFT= fft(DataSlice,NFFT,2);
    
    %Chops down the FFT to the right size?
    DataSliceFFT = DataSliceFFT(:,1:NFFT/2+1);

    %Trims down the data slice to just the data we're interested in and
    %stores the data into the 3D FFT matrix
    txEy_FFT(intIndex,:,:) = DataSliceFFT(:,freqStartIndex:freqStopIndex);
    
    %Removes the data slice FFT
    clear DataSliceFFT;

end 


%% ======================================
%PERFORMS FFT: transmittedHxr
%========================================

%Performs 2D FFT
for intIndex = 1:size(transmittedHxr,1)
    %Retrieves a slice of data from the matrix
    DataSlice = reshape(transmittedHxr(intIndex,:,:),[size(transmittedHxr,2),size(transmittedHxr,3)]);

    %Performs FFT on the data slice
    DataSliceFFT= fft(DataSlice,NFFT,2);
    
    %Chops down the FFT to the right size?
    DataSliceFFT = DataSliceFFT(:,1:NFFT/2+1);

    %Trims down the data slice to just the data we're interested in and
    %stores the data into the 3D FFT matrix
    txHx_FFT(intIndex,:,:) = DataSliceFFT(:,freqStartIndex:freqStopIndex);
    
    %Removes the data slice FFT
    clear DataSliceFFT;

end 

%% ======================================
%PERFORMS FFT: transmittedEyr
%========================================

%Performs 2D FFT
for intIndex = 1:size(transmittedHyr,1)
    %Retrieves a slice of data from the matrix
    DataSlice = reshape(transmittedHyr(intIndex,:,:),[size(transmittedHyr,2),size(transmittedHyr,3)]);

    %Performs FFT on the data slice
    DataSliceFFT= fft(DataSlice,NFFT,2);
    
    %Chops down the FFT to the right size?
    DataSliceFFT = DataSliceFFT(:,1:NFFT/2+1);

    %Trims down the data slice to just the data we're interested in and
    %stores the data into the 3D FFT matrix
    txHy_FFT(intIndex,:,:) = DataSliceFFT(:,freqStartIndex:freqStopIndex);
    
    %Removes the data slice FFT
    clear DataSliceFFT;

end 


%% ======================================
%PERFORMS FFT: reflectedExr
%========================================

%Performs 2D FFT
for intIndex = 1:size(reflectedExr,1)
    %Retrieves a slice of data from the matrix
    DataSlice = reshape(reflectedExr(intIndex,:,:),[size(reflectedExr,2),size(reflectedExr,3)]);

    %Performs FFT on the data slice
    DataSliceFFT= fft(DataSlice,NFFT,2);
    
    %Chops down the FFT to the right size?
    DataSliceFFT = DataSliceFFT(:,1:NFFT/2+1);

    %Trims down the data slice to just the data we're interested in and
    %stores the data into the 3D FFT matrix
    refEx_FFT(intIndex,:,:) = DataSliceFFT(:,freqStartIndex:freqStopIndex);
    
    %Removes the data slice FFT
    clear DataSliceFFT;

end 


%% ======================================
%PERFORMS FFT: reflectedEyr
%========================================

%Performs 2D FFT
for intIndex = 1:size(reflectedEyr,1)
    %Retrieves a slice of data from the matrix
    DataSlice = reshape(reflectedEyr(intIndex,:,:),[size(reflectedEyr,2),size(reflectedEyr,3)]);

    %Performs FFT on the data slice
    DataSliceFFT= fft(DataSlice,NFFT,2);
    
    %Chops down the FFT to the right size?
    DataSliceFFT = DataSliceFFT(:,1:NFFT/2+1);

    %Trims down the data slice to just the data we're interested in and
    %stores the data into the 3D FFT matrix
    refEy_FFT(intIndex,:,:) = DataSliceFFT(:,freqStartIndex:freqStopIndex);
    
    %Removes the data slice FFT
    clear DataSliceFFT;

end 

%% ======================================
%PERFORMS FFT: reflectedHxr
%========================================

%Performs 2D FFT
for intIndex = 1:size(reflectedHxr,1)
    %Retrieves a slice of data from the matrix
    DataSlice = reshape(reflectedHxr(intIndex,:,:),[size(reflectedHxr,2),size(reflectedHxr,3)]);

    %Performs FFT on the data slice
    DataSliceFFT= fft(DataSlice,NFFT,2);
    
    %Chops down the FFT to the right size?
    DataSliceFFT = DataSliceFFT(:,1:NFFT/2+1);

    %Trims down the data slice to just the data we're interested in and
    %stores the data into the 3D FFT matrix
    refHx_FFT(intIndex,:,:) = DataSliceFFT(:,freqStartIndex:freqStopIndex);
    
    %Removes the data slice FFT
    clear DataSliceFFT;

end 

%% ======================================
%PERFORMS FFT: reflectedEyr
%========================================

%Performs 2D FFT
for intIndex = 1:size(reflectedHyr,1)
    %Retrieves a slice of data from the matrix
    DataSlice = reshape(reflectedHyr(intIndex,:,:),[size(reflectedHyr,2),size(reflectedHyr,3)]);

    %Performs FFT on the data slice
    DataSliceFFT= fft(DataSlice,NFFT,2);
    
    %Chops down the FFT to the right size?
    DataSliceFFT = DataSliceFFT(:,1:NFFT/2+1);

    %Trims down the data slice to just the data we're interested in and
    %stores the data into the 3D FFT matrix
    refHy_FFT(intIndex,:,:) = DataSliceFFT(:,freqStartIndex:freqStopIndex);
    
    %Removes the data slice FFT
    clear DataSliceFFT;

end 

%Performs addition on the light intensity on each electric field direction
if (intEzCoeff == 0)
    %Light traveling in -Z direction
    tx_FFT =  real(1/2.*txEy_FFT .* conj(txHx_FFT)) - real(1/2.* txEx_FFT .* conj(txHy_FFT));
elseif (intExCoeff == 0)
    %Light traveling in +X direction
    tx_FFT =  -1*real(1/2.*txEz_FFT .* conj(txHy_FFT)) + real(1/2.* txEy_FFT .* conj(txHz_FFT));
else
    %Prints error message
    beep;
    fprintf('\n\nYour coefficients could be messed up! Check them, or write a new case for light direction\n\n');
end

%Performs addition on the light intensity on each electric field direction
if (intEzCoeff == 0)
    %Light traveling in -Z direction
    ref_FFT = real(1/2.* refEx_FFT .* conj(refHy_FFT)) - real(1/2.*refEy_FFT .* conj(refHx_FFT));
elseif (intExCoeff == 0)
    %Light traveling in +X direction
    ref_FFT = -1*real(1/2.* refEz_FFT .* conj(refHy_FFT)) + real(1/2.*refEy_FFT .* conj(refHz_FFT));
else
    %Prints error message
    beep;
    fprintf('\n\nYour coefficients could be messed up! Check them, or write a new case for light direction\n\n');
end

