function [fftResult] = PerformFFT(matTransform, dT, lambdaMin, lambdaMax)
%Defines the speed of light
c = 299792458;

%Sets the frequencies we're interested in
freqStart = c/lambdaMax; % frequency in Hz
freqStop = c/lambdaMin;

%Figures out the length of the sample
L = size(matTransform,3);

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
%fftResult = zeros(size(matTransform,1),size(matTransform,2),(freqStopIndex-freqStartIndex+1));
%fftResult = zeros(size(matTransform,1),size(matTransform,2),size(matTransform,3));


%% Performs 2D FFT
for intIndex = 1:size(matTransform,1)
    %Retrieves a slice of data from the matrix
    DataSlice = reshape(matTransform(intIndex,:,:),[size(matTransform,2),size(matTransform,3)]);

    %Performs FFT on the data slice
    DataSliceFFT= fft(DataSlice,NFFT,2);
    
    %Chops down the FFT to the right size
    DataSliceFFT = DataSliceFFT(:,1:NFFT/2+1);

    %Trims down the data slice to just the data we're interested in and
    %stores the data into the 3D FFT matrix
    fftResult(intIndex,:,:) = DataSliceFFT(:,freqStartIndex:freqStopIndex);
    
    %Removes the data slice FFT
    clear DataSliceFFT;
end