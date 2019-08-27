%% Select EEG file
working_dir=pwd;
current_dir='C:\SleepData\DataFiles';
cd(current_dir);
[fileName, pathName] = uigetfile({'*.ncs','Neuralynx CSC File (*.ncs)'},...
    'Select the EEG data file');
cscFile= fullfile(pathName, fileName);
cd(working_dir);
fileLabel = strrep(fileName, '.ncs', '');
clear fileName pathName

%% Load the CSC file:
[Timestamps, SampleFrequencies, Samples] = Nlx2MatCSC(cscFile, [1 0 1 0 1], 0, 1, [] );
cscSamplingRate = SampleFrequencies(1);
clear SampleFrequencies
[binSize,numBins]=size(Samples);
newM = binSize*numBins;
Samples = reshape(Samples, newM, 1);
interTimestamps = zeros(newM, 1);

%% Interpolate time stamps:
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
clear trange Timestamps idx

EEG_timestamps=interTimestamps/1000000; %Convert from usec to seconds.
clear interTimestamps

events = events/1000000; %convert event time stamps from usec to seconds.
numBins = size(events,1); 

powerSpectraResults = cell(numBins,1);
dt = 1/cscSamplingRate; % Define the sampling interval.
fNyquist = cscSamplingRate/2; % Determine the Nyquist frequency.
maxFrequency = 30; %Define the max frequency you want to look at.

%% Define frequency bands:
delta = [0.4 4];
theta = [5 9];
sigma = [10 14]; 
beta = [15 20]; 
freqBands = [delta; theta; sigma; beta];
numBands = size(freqBands,1);
bandPowerResults = zeros(numBins, numBands);
df = 0.2; % Minimum df = 1/duration
faxis = (0.4:df:30); % Defines the frequencies to calculate PSDs for in 'pwelch'

%% Calculate power spectra and band power for each user defined segment of data:
for i = 1:numBins
    diff1 = EEG_timestamps - events(i,1);
    [~,startIdx]=min(abs(diff1));
    startTime = EEG_timestamps(startIdx);

    diff2 = EEG_timestamps - events(i,2);
    [~,endIdx]=min(abs(diff2));
    endTime = EEG_timestamps(endIdx);

    x = Samples(startIdx:endIdx);

%     % Compute the power spectrum of the Hann tapered data:    
%     duration = endTime - startTime;
%     df = 1/duration;
%     faxis = (0:df:fNyquist); % Construct the frequency axis.
%     lastFreqIndex = find(faxis <= maxFrequency, 1, 'last');
%     faxis = faxis(1:lastFreqIndex);
%     xh = hann(length(x)).*x; 
%     Sxx = 2*dt^2/duration * fft(xh).*conj(fft(xh)); % Compute the power spectrum of Hann tapered data.
%     Sxx=real(Sxx); %Ignores imaginary component.
%     %Should not convert to dB power before calculating band power and/or taking averages:
%     %Sxx=10*log10(real(Sxx)); % ... convert to a decibel scale. Ignores imaginary component.
%     %Sxx = Sxx(1:length(x)/2+1); % ... ignore negative frequencies. Only needed if not limiting frequency range like next line:		
%     faxis = faxis(:); %Force into column vector
%     Sxx = Sxx(1:lastFreqIndex);
%     powerSpectraResults{i,1} = [faxis Sxx];
    
    % Calculate the power spectral density (PSD) in power/Hz:
    window = cscSamplingRate;
    overlap = floor(window/2);
    [pxx,f] = pwelch(x,window,overlap,faxis,cscSamplingRate);
    powerSpectraResults{i,1} = [f' pxx'];
    
    % Calculate the band power:
    width = diff(f');
    lastRectWidth = 0;  % Don't include last point of PSD data.
    width = [width; lastRectWidth];
    for j = 1:numBands
        idx1 = find(f<=freqBands(j,1), 1, 'last' );
        idx2 = find(f>=freqBands(j,2), 1, 'first');
        bandPowerResults(i,j) = width(idx1:idx2)'*pxx(idx1:idx2)';
        clear idx1 idx2
    end
end