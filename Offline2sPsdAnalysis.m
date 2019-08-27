function Offline2sPsdAnalysis
%%                Version 1.0 Created by Brooks A. Gross
%                                08.28.2014
%                            SLEEP AND MEMORY LAB
%                           UNIVERSITY OF MICHIGAN
%
%% DESCRIPTION:
%   'Offline2sPsdAnalysis.m' generates power spectral density for each 2 
%   second epoch.  Results are saved to a .xlsx file. The power estimates 
%   for the defined bandwidths are written to one spreadsheet, and the
%   spectral power is written to another spreadsheet.  EEG data in the .NCS
%   file (CSC) CANNOT have any gaps in time.  The user must split the .NCS
%   to remove gaps and run individual sections separately.

%% DEFINE CONSTANTS
    
% EEG bandwidths:
D_lo = 0.1; % Specify low end of Delta band
D_hi = 4; % Specify high end of Delta band
T_lo = 4.1; % Specify low end of Theta band
T_hi = 10;   % Specify high end of Theta band
S_lo = 10.1; % Specify low end of Sigma band
S_hi = 15;  % Specify high end of Sigma band
B_lo = 15.1; 
B_hi = 20;

% EEG filter set:
EEG_Fc = 30; % The EEG low-pass cut-off frequency. Default is 30 Hz.
EEG_highpass_enable = 0; % Set to '0' to turn off high-pass filter, '1' to turn on.
EEG_HP_Fc = 1; % The EEG high-pass cut-off frequency. Default is 1 Hz.
EEG_Notch_enable = 0; % Set to '0' to turn off notch filter, '1' to turn on.

% -------------------------------------------------------------------------
% SELECT DATA FILES:

% Select EEG file
working_dir=pwd;
current_dir='C:\SleepData\DataFiles';
cd(current_dir);
[fileName, pathName] = uigetfile({'*.ncs','Neuralynx CSC File (*.ncs)'},...
    'Select the EEG data file');
cscFile= fullfile(pathName, fileName);
cd(working_dir);
fileLabel = strrep(fileName, '.ncs', '');
clear fileName pathName

% Load the CSC file:
[Timestamps, SampleFrequencies, Samples] = Nlx2MatCSC(cscFile, [1 0 1 0 1], 0, 1, [] );
SampFreq1 = SampleFrequencies(1);
clear SampleFrequencies
[binSize,numBins]=size(Samples);
newM = binSize*numBins;
Samples = reshape(Samples, newM, 1);
interTimestamps = zeros(newM, 1);

%Interpolate time stamps:
idx = 1;
for i = 1:numBins
  if i < numBins
    t1 = Timestamps(i);
    t2 = Timestamps(i+1);
    interval = (t2-t1)/binSize;
    trange =([t1 : interval : t2-interval]);
    interTimestamps(idx:idx+binSize-1,1) = trange;
  else
    t1 = Timestamps(i);
    t2 = t1+interval*binSize;
    trange =([t1 :interval : t2-interval]);
    interTimestamps(idx:idx+binSize-1,1) = trange;
  end
  idx = idx + binSize;
end
clear trange Timestamps
interTimestamps=interTimestamps/1000000; %Convert from usec to seconds.


%Calculate down-sampling factor:
DS = (1:1:32);
DSampSF = SampFreq1./DS;
indSampfactor = find(DSampSF >= 1000);
reducedSamplingRate = DSampSF(indSampfactor(end));
sampfactor = DS(indSampfactor(end));
msgbox({['Orginal Sampling Rate:  ' num2str(SampFreq1) 'Hz'];...
    ['Down-Sampled Sampling Rate:  ' round(num2str(reducedSamplingRate)) 'Hz'];...
    ['Sampling Factor:  ' num2str(sampfactor) '']});
    
physInput = 2;  %Needed to select proper error box in HeaderADBit.
ADBit2uV = HeaderADBit(cscFile, physInput);    %Calls a function to extract the AD Bit Value.
Samples = Samples * ADBit2uV;   %Convert EEG amplitude of signal from AD Bits to microvolts.
%  Low pass filter for EEG signals
[Blow,Alow]=ellip(7,1,60, EEG_Fc/(SampFreq1/2));           % Default setting implements low pass filter with 30hz cutoff
filtered_samples=filter(Blow,Alow,Samples);
clear Samples

%  OPTIONAL highpass filter for EEG signals
if EEG_highpass_enable>0
    [EEG_Bhi,EEG_Ahi] = ellip(7,1,60, EEG_HP_Fc/(SampFreq1/2),'high');   % Default is OFF
    filtered_samples = filter(EEG_Bhi,EEG_Ahi, filtered_samples);
end
%  OPTIONAL 60Hz Notch filter for EEG signals
if EEG_Notch_enable > 0
    wo = 60/(SampFreq1/2);
    [B_EEG_Notch,A_EEG_Notch] =  iirnotch(wo, wo/35);   % Default is OFF
    filtered_samples = filter(B_EEG_Notch,A_EEG_Notch, filtered_samples);
end
interTimestamps = interTimestamps(1:sampfactor:end);
Samples = filtered_samples(1:sampfactor:end);
clear physInput ADBit2uV filtered_samples

dt = 1/reducedSamplingRate; % Define the sampling interval.
fNyquist = reducedSamplingRate/2; % Determine the Nyquist frequency.

reducedSamplingRate = floor(reducedSamplingRate);
epochSize = 2*reducedSamplingRate;
num2sEpochs = floor(size(interTimestamps,1)/epochSize);
epochDuration = interTimestamps(epochSize+1) - interTimestamps(1);

df = 1/epochDuration; % Determine the frequency resolution.

faxis = (0:df:fNyquist); % Construct the frequency axis.
lastFreqIndex = find(faxis <= EEG_Fc, 1, 'last');
faxis = faxis(1:lastFreqIndex);
epochPowerSpectrum = zeros(num2sEpochs, lastFreqIndex);

warning('off', 'signal:spectrum:obsoleteFunction');

epochTime = zeros(num2sEpochs,1);
deltaPower = zeros(num2sEpochs,1);
thetaPower = zeros(num2sEpochs,1);
sigmaPower = zeros(num2sEpochs,1);
betaPower = zeros(num2sEpochs,1);
st_power = zeros(num2sEpochs,1);
dt_ratio = zeros(num2sEpochs,1);

for i = 1:num2sEpochs
    startPoint = (i-1)*epochSize + 1;
    endPoint = startPoint+ epochSize - 1;
    eeg2sData = Samples(startPoint:endPoint);
    epochTime(i) = interTimestamps(startPoint);        
    % This is calculating EEG power in frequency domain
    fft_in=double(eeg2sData);
%     windowsize =length(fft_in);
%     df = epochSize;
%     if df < windowsize
%         df = windowsize;
%     end

    % Compute the power spectrum of the Hann tapered data:
    
    xh = hann(length(eeg2sData)).*eeg2sData; 
    Sxx = 2*dt^2/epochDuration * fft(xh).*conj(fft(xh)); % Compute the power spectrum of Hann tapered data.
%     Sxx=10*log10(Sxx); % ... convert to a decibel scale.
    Sxx = Sxx(1:length(eeg2sData)/2+1); % ... ignore negative frequencies.		    
    Sxx = real(Sxx(1:lastFreqIndex));
    epochPowerSpectrum(i, :) = Sxx';
    % SPECTRUM is obsolete. Replaced with above 4 lines that does the equivalent.
%     [Pxx2,F2]=spectrum(fft_in,df,0,ones(windowsize,1),reducedSamplingRate);
    % ******  [P,F] = SPECTRUM(X,NFFT,NOVERLAP,WINDOW,reducedSamplingRate)

    %For the EEG signal
    index_delta=[];index_theta=[];index_sigma=[]; index_beta=[];
    index_delta=find(faxis(1)+D_lo< faxis & faxis < faxis(1)+D_hi);      % Default delta band 0.4 -4 Hz
    index_theta=find(faxis(1)+T_lo< faxis & faxis < faxis(1)+T_hi);    % Default theta band 5-9 Hz
    index_sigma=find(faxis(1)+S_lo< faxis & faxis < faxis(1)+S_hi);     % Default sigma band 10-14 Hz
    index_beta =find(faxis(1)+B_lo< faxis & faxis < faxis(1)+B_hi);     % Default Beta band 15-20 Hz
%     spectPower(i,:) = Sxx(:,1)';
    deltaPower(i)=sum(Sxx(index_delta))/epochSize *2;
    thetaPower(i)=sum(Sxx(index_theta))/epochSize *2;
    sigmaPower(i)=sum(Sxx(index_sigma))/epochSize *2;
    betaPower(i)=sum(Sxx(index_beta))/epochSize *2;
    clear index_delta index_theta index_sigma index_beta
    st_power(i)=abs(sigmaPower(i).*thetaPower(i));   % Used to indicate waking
    dt_ratio(i)=abs(deltaPower(i)./thetaPower(i));   
end      
clear EMG_SAMPLES EMG_TIMESTAMPS Samples interTimestamps  


%Write header for PSD file before entering loop to calculate PSd values for
%all 2s epochs
c = clock;
dt = datestr(c,'mmddyy-HHMM');
resultExcelFileName = ['Psd2sEpochsBasedOn' fileLabel '_' dt '.xlsx'];
warning off MATLAB:xlswrite:AddSheet
sheetName = 'Files used';
xlswrite(resultExcelFileName,cellstr('EEG'), sheetName, 'A2');
xlswrite(resultExcelFileName,cellstr(cscFile), sheetName, 'B2');

%Write PSD information to a new sheet in the Excel file.
sheetName = 'PSDs-2sEpochs';
xlswrite(resultExcelFileName,cellstr('Timestamp'), sheetName, 'A1');
xlswrite(resultExcelFileName,cellstr('DeltaPower'), sheetName, 'B1');
xlswrite(resultExcelFileName,cellstr('ThetaPower'), sheetName, 'C1');
xlswrite(resultExcelFileName,cellstr('SigmaPower'), sheetName, 'D1');
xlswrite(resultExcelFileName,cellstr('BetaPower'), sheetName, 'E1');
xlswrite(resultExcelFileName,cellstr('SxT_Power'), sheetName, 'F1');
xlswrite(resultExcelFileName,cellstr('DT_Ratio'), sheetName, 'G1');
xlswrite(resultExcelFileName,[epochTime deltaPower thetaPower sigmaPower betaPower st_power dt_ratio],sheetName, 'A2');

%Write full spectrum to a new sheet in the Excel file.
sheetName = 'spectrumVsFreq';
xlswrite(resultExcelFileName,cellstr('Timestamp'), sheetName, 'A1');
xlswrite(resultExcelFileName,faxis, sheetName, 'B1');
xlswrite(resultExcelFileName,epochTime,sheetName, 'A2');
xlswrite(resultExcelFileName,epochPowerSpectrum,sheetName, 'B2');
