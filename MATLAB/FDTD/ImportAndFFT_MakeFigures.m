    

%% =======================
%Compares to reference data
%=========================

% Loads the reference data
load (char(arrDataPath(intDataParseIndex,11)));

%Defines the speed of light
c = 299792458;

%Sets the frequencies we're interested in
freqStart = c/(cell2mat(arrDataPath(intDataParseIndex,8))*1e-9); % frequency in Hz
freqStop = c/(cell2mat(arrDataPath(intDataParseIndex,7))*1e-9);

%Interpolates the transmitted light to have the same scale as the
%transmitted Free Space
X = freqStart:(freqStop-freqStart)/(length(txAvg)-1):freqStop;                 %Old X axis spacing
Xq = freqStart:(freqStop-freqStart)/(length(txAvgFS)-1):freqStop;              %New X axis spacing
interpRatioTx = (length(txAvgFS)/length(txAvg))^2;                                     %Calculates the interpolation ratio 
txAvg = interp1(X,txAvg,Xq)./interpRatioTx;


%Interpolates the reflected light
X = freqStart:(freqStop-freqStart)/(length(refAvg)-1):freqStop;                 %Old X axis spacing
Xq = freqStart:(freqStop-freqStart)/(length(refAvgFS)-1):freqStop;              %New X axis spacing
interpRatioRef = (length(refAvgFS)/length(refAvg))^2;                                   %Calculates the interpolation ratio 
refAvg = interp1(X,refAvg,Xq)./interpRatioRef;                                  


% Creates the ratio of transmitted light vs. the reference transmitted light
txRatio = sum(txAvg ./ txAvgFS) / length(txAvgFS);

%Performs error checking by making sure that reflected + transmitted = reference
errFFTPct = 100*sum(1-((txAvg +refAvg)./ txAvgFS)) / length(txAvgFS);
errFFT = ((txAvg +refAvg)./ txAvgFS);
txAvg100 = txAvg ./ errFFT;
refAvg100 = refAvg ./ errFFT;


%% ===========================
%Generate Figures
%=============================

%% Figure: Transmitted, Reflected and Free Space Comparison
figure
hold on
plot(freq/1e9, txAvgFS,'b')
plot(freq/1e9, txAvg,'r')
plot(freq/1e9,refAvg,'g')
legend('Free Space','Transmitted','Reflected');
xlabel ('Frequency (GHz)');
ylabel ('Intensity');
title(sprintf('Experiment: %s\nPercent Transmitted: %s%% \nPercent Error: %s%%',char(arrDataPath(intDataParseIndex,2)),num2str(round(txRatio*1000)/10),num2str(abs(round(errFFTPct*100)/100))));

%Sets the font to 16, bold and saving the figure
set(findall(gcf,'-property','FontSize'),'FontSize',16) 
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
print(gcf,'-dpng','-r600',['N:\Kat FDTD Data\Figures\' char(arrDataPath(intDataParseIndex,2)) ', Intensity.png'])

%% Figure: Percent light transmitted
percenttransmitted = txAvg ./ txAvgFS;
percentreflected = refAvg ./ txAvgFS;

% for i = 1:length(percenttransmitted)
%     reversepercenttransmitted(i) =  percenttransmitted(length(percenttransmitted) - i+1);
% end
c = 299792458;
wavelength = c./freq.*1e9; 

figure
hold on
plot(wavelength, percenttransmitted.*100,'b')
% ylim([80 100])

%--legend('Free Space','Transmitted','Reflected');
xlabel ('Wavelength (nm)');
ylabel ('Percent Transmitted');
% title(sprintf('Experiment: %s\nPercent Transmitted: %s%% \nPercent Error: %s%%',char(arrDataPath(intDataParseIndex,2)),num2str(round(txRatio*1000)/10),num2str(abs(round(errFFTPct*100)/100))));
title(sprintf('Experiment: %s\n',char(arrDataPath(intDataParseIndex,2))))
%Sets the font to 16, bold and saving the figure
set(findall(gcf,'-property','FontSize'),'FontSize',16) 
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
print(gcf,'-dpng','-r600',['N:\Kat FDTD Data\Figures\' char(arrDataPath(intDataParseIndex,2)) ', Percent Transmitted.png'])

figure
hold on
plot(wavelength, percentreflected.*100,'b')
% ylim([80 100])

%--legend('Free Space','Transmitted','Reflected');
xlabel ('Wavelength (nm)');
ylabel ('Percent Reflected');
% title(sprintf('Experiment: %s\nPercent Transmitted: %s%% \nPercent Error: %s%%',char(arrDataPath(intDataParseIndex,2)),num2str(round(txRatio*1000)/10),num2str(abs(round(errFFTPct*100)/100))));
title(sprintf('Experiment: %s\n',char(arrDataPath(intDataParseIndex,2))))
%Sets the font to 16, bold and saving the figure
set(findall(gcf,'-property','FontSize'),'FontSize',16) 
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold') 
print(gcf,'-dpng','-r600',['N:\Kat FDTD Data\Figures\' char(arrDataPath(intDataParseIndex,2)) ', Percent Reflected.png'])

%% Closes all figures
% close all;