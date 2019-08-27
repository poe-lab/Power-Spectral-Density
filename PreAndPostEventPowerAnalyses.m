function PreAndPostEventPowerAnalyses
%%                Version 1.0 Created by Brooks A. Gross
%                                11.10.2014
%                            SLEEP AND MEMORY LAB
%                           UNIVERSITY OF MICHIGAN
%
%% DESCRIPTION:
%   'PreAndPostEventPowerAnalyses.m' generates power spectral density for each 
%   user-defined pre and post event epoch.  Results are saved to a .xlsx file. The power estimates 
%   for the defined bandwidths are written to one spreadsheet, and the
%   spectral power is written to another spreadsheet.  EEG data in the .NCS
%   file (CSC) CANNOT have any gaps in time.  The user must split the .NCS
%   to remove gaps and run individual sections separately.
%--04.15.2015: Added a time shift to the events to compensate for the delay
%               between the TTL pulse writing the event to the Event file.

%% DEFINE CONSTANTS
    
% EEG bandwidths:
delta = [0.4 4];
theta = [5 9];
sigma = [10 14]; 
beta = [15 20]; 
freqBands = [delta; theta; sigma; beta];


% EEG filter set:
EEG_Fc = 30; % The EEG low-pass cut-off frequency. Default is 30 Hz.
% EEG_highpass_enable = 0; % Set to '0' to turn off high-pass filter, '1' to turn on.
EEG_HP_Fc = 0.4; % The EEG high-pass cut-off frequency. Default is 1 Hz.
% EEG_Notch_enable = 0; % Set to '0' to turn off notch filter, '1' to turn on.

% -------------------------------------------------------------------------
% SELECT DATA FILES:
%Request user input to select event file:
working_dir = pwd;
current_dir = 'C:\SleepData\Datafiles'; % This is the default directory that opens.
cd(current_dir);
[filename, pathname] = uigetfile('*.nev', 'Select a Cheetah event file'); %This waits for user input and limits selection to .nev files.
% Check for whether or not a file was selected
if isempty(filename) || isempty(pathname)
uiwait(errordlg('You need to select an event file. Please try again',...
'ERROR','modal'));
cd(working_dir);
else
cd(working_dir);
NevFile= fullfile(pathname, filename);
end
%load event file
ExtractHeader = 0;  % 0 for no and 1 for yes
ExtractMode = 1;  %Extract all data points
[eventTimeStamps, EventStrings] =Nlx2MatEV( NevFile, [1 0 0 0 1], ExtractHeader, ExtractMode, []);

% Remove all events that do not pertain to air puffing:
k = strfind(EventStrings, '10 sec');

numEvents=size(k,1);
for i = 1:numEvents
    if isempty(k{i,1})
        eventTimeStamps(i)=0;
    end
end
m = eventTimeStamps>0;
EventStrings = EventStrings(m);
eventTimeStamps = eventTimeStamps(m)';
eventTimeStamps = eventTimeStamps/1000000;
clear k numEvents EventStrings

numEvents = size(eventTimeStamps,1);
lapsedTime = eventTimeStamps(1) + 5; %Five second blackout period for air puffs
realAirpuffs = zeros(numEvents,1);
realAirpuffs(1) = 1;
for i = 2:numEvents
    if eventTimeStamps(i) < lapsedTime
        realAirpuffs(i) = 0;
    else
        realAirpuffs(i) = 1;
        lapsedTime = eventTimeStamps(i) + 5;
    end
end
eventTimeStamps = eventTimeStamps(logical(realAirpuffs));
numEvents = size(eventTimeStamps,1);
%% Label Events with States Based off of a Sleep Scored File:
% Call the sleep scored file name:
current_dir='C:\SleepData';
cd(current_dir);
[scoredFile, scoredPath] = uigetfile({'*.xls','Excel 1997-2003 File (*.xls)'},'Select the Sleep Scored File');
if isequal(scoredFile,0) || isequal(scoredPath,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    sleepScoredFile= fullfile(scoredPath, scoredFile);

end

%Load sleep scored file:
try
    [numData, stringData] = xlsread(sleepScoredFile);
catch %#ok<*CTCH>
    uiwait(errordlg('Check if the file is saved in Microsoft Excel format.',...
     'ERROR','modal'));
end

%Detect if states are in number or 2-letter format:
if isequal(size(numData,2),3)
    
    scoredStates = numData(:,2:3);
    clear numData stringData
else
    scoredStates = numData(:,2);
    clear numData
    stringData = stringData(3:end,3);
    [stateNumber] = stateLetter2NumberConverter(stringData);
    scoredStates = [scoredStates stateNumber];
    clear stateNumber stringData
end

%label events with states:
labeledEventStates = zeros(numEvents,1);
stateTargetInterval=[];
lengthScoredStates =size(scoredStates,1);
lengthScoredSubStates =[];
for i = 1:8
    isoCount = 1;
    lengthScoredSubStates = lengthScoredStates;
    scoredSubStates = scoredStates;
    for j = 1:lengthScoredSubStates
        if isequal(scoredStates(j,2),i)
            scoredSubStates(j,2) = 1;
        else
            scoredSubStates(j,2) = 0;   
        end  
    end
    firstIsoInd = find(scoredSubStates(:,2)==1, 1);
    if isempty(firstIsoInd)
        scoredSubStates = [];
    elseif firstIsoInd > 1 %First of target stage detected is not at index = 1.
        scoredSubStates(1:firstIsoInd-1,:) = [];
    end
    lengthScoredSubStates =size(scoredSubStates,1); %Recalculate the length of the array due to removal of initial rows.
    if lengthScoredSubStates < 2
        stateTargetInterval = [];
    else
        stateTargetInterval(isoCount,1) = scoredSubStates(1,1); 
        %The following FOR loop generates isolated intervals based on user-selected states:
        for j = 2:lengthScoredSubStates   
            if isequal(scoredSubStates(j,2),1)
                if isequal(scoredSubStates(j-1,2),0)
                    stateTargetInterval(isoCount,1) = scoredSubStates(j,1); %Looking at time
                end
                if isequal(j, lengthScoredSubStates)
                    stateTargetInterval(isoCount,2) = scoredSubStates(j,1); %Looking at time
                end
            elseif isequal(scoredSubStates(j,2),0) && isequal(scoredSubStates(j-1,2),1)
                stateTargetInterval(isoCount,2) = scoredSubStates(j,1); %Looking at time
                isoCount = isoCount + 1;
            end
        end
    end

    [lengthIsoArray, ~] = size(stateTargetInterval);
    if isequal(0,lengthIsoArray)
    else
        for m = 1:lengthIsoArray % Extract all of the sub-intervals for states containing events.
           subIntervalIndx = find(eventTimeStamps >= stateTargetInterval(m,1) & eventTimeStamps <= stateTargetInterval(m,2));
           if isempty(subIntervalIndx)
           else
               labeledEventStates(subIntervalIndx) = i;
           end
        end
    end
    clear stateTargetInterval subIntervalIndx
end



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

% Load the CSC file:
[Timestamps, SampleFrequencies, Samples] = Nlx2MatCSC(cscFile, [1 0 1 0 1], 0, 1, [] );
cscSamplingRate = SampleFrequencies(1);
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
DSampSF = cscSamplingRate./DS;
indSampfactor = find(DSampSF >= 1000);
reducedSamplingRate = DSampSF(indSampfactor(end));
sampfactor = DS(indSampfactor(end));
msgbox({['Orginal Sampling Rate:  ' num2str(cscSamplingRate) 'Hz'];...
    ['Down-Sampled Sampling Rate:  ' round(num2str(reducedSamplingRate)) 'Hz'];...
    ['Sampling Factor:  ' num2str(sampfactor) '']});
    
physInput = 2;  %Needed to select proper error box in HeaderADBit.
ADBit2uV = HeaderADBit(cscFile, physInput);    %Calls a function to extract the AD Bit Value.
Samples = Samples * ADBit2uV;   %Convert EEG amplitude of signal from AD Bits to microvolts.
%Design bandpass filter:
[z, p, k] = ellip(7,1,60, [EEG_HP_Fc EEG_Fc]/(cscSamplingRate/2),'bandpass');
[sos, g] = zp2sos(z,p,k);
%Apply bandpass filter to EEG data:
filtered_samples = filtfilt(sos, g, Samples);
clear Samples

% %  OPTIONAL highpass filter for EEG signals
% if EEG_highpass_enable>0
%     [EEG_Bhi,EEG_Ahi] = ellip(7,1,60, EEG_HP_Fc/(cscSamplingRate/2),'high');   % Default is OFF
%     filtered_samples = filtfilt(sos,g, filtered_samples);
% end
% %  OPTIONAL 60Hz Notch filter for EEG signals
% if EEG_Notch_enable > 0
%     wo = 60/(cscSamplingRate/2);
%     [B_EEG_Notch,A_EEG_Notch] =  iirnotch(wo, wo/35);   % Default is OFF
%     filtered_samples = filtfilt(B_EEG_Notch,A_EEG_Notch, filtered_samples);
% end
interTimestamps = interTimestamps(1:sampfactor:end);
Samples = filtered_samples(1:sampfactor:end);
clear physInput ADBit2uV filtered_samples

dt = 1/reducedSamplingRate; % Define the sampling interval.
fNyquist = reducedSamplingRate/2; % Determine the Nyquist frequency.

reducedSamplingRate = floor(reducedSamplingRate);
epochSize = 2*reducedSamplingRate;
epochDuration = interTimestamps(epochSize+1) - interTimestamps(1);

df = 1/epochDuration; % Determine the frequency resolution.

freqVector = (0:df:fNyquist); % Construct the frequency axis.
lastFreqIndex = find(freqVector <= EEG_Fc, 1, 'last');
freqVector = freqVector(1:lastFreqIndex);
freqVector = freqVector(:); %Force into column vector


% deltaPower = zeros(numEvents,2);
% thetaPower = zeros(numEvents,2);
% sigmaPower = zeros(numEvents,2);
% betaPower = zeros(numEvents,2);

epochPowerSpectrum = zeros(lastFreqIndex,numEvents, 2);

for i = 1:numEvents
    shiftTS = floor(epochSize/20); %Shifts the time of the event to account for delay.
    %find time stamp closest to event:
    lowerLfpTimstampIndex = find(ge(eventTimeStamps(i), interTimestamps), 1, 'last') - shiftTS;
    upperLfpTimstampIndex = find(le(eventTimeStamps(i), interTimestamps), 1, 'first') - shiftTS;
    
    startPoint = [(lowerLfpTimstampIndex - epochSize +1), upperLfpTimstampIndex];
    endPoint = [lowerLfpTimstampIndex, (upperLfpTimstampIndex + epochSize - 1)];
    for j = 1:2
        eeg2sData = Samples(startPoint(j):endPoint(j));
      
    % Compute the power spectrum of the Hann tapered data:
        xh = hann(length(eeg2sData)).*eeg2sData; 
        Sxx = 2*dt^2/epochDuration * fft(xh).*conj(fft(xh)); % Compute the power spectrum of Hann tapered data.
    %     Sxx=10*log10(Sxx); % ... convert to a decibel scale.
        Sxx = Sxx(1:length(eeg2sData)/2+1); % ... ignore negative frequencies.		    
        Sxx = real(Sxx(1:lastFreqIndex));
        epochPowerSpectrum(:, i, j) = Sxx;
    end
end
numBands = size(freqBands,1);
allEvntsBandPwr = zeros(numBands, numEvents,2);
width = diff(freqVector);
lastRectWidth = 0;  % Don't include last point of PSD data.
width = [width; lastRectWidth];
for i = 1:numBands
    idx1 = find(freqVector<=freqBands(i,1), 1, 'last' );
    idx2 = find(freqVector>=freqBands(i,2), 1, 'first');
    for j = 1:2
        allEvntsBandPwr(i,:,j) = width(idx1:idx2)'*epochPowerSpectrum(idx1:idx2,:,j);
    end
    clear idx1 idx2
end
             
clear Samples interTimestamps  

%Calculate average spectrum for sleep:
sleepIdx = ismember(labeledEventStates,2);
sizeSleepEvents = sum(sleepIdx);

preEventSpectrumSleepOnly = epochPowerSpectrum(:,sleepIdx,1);
meanSleepPreEventSpectrum = mean(preEventSpectrumSleepOnly,2);
stdDevSleepPreEventSpectrum = std(preEventSpectrumSleepOnly,0,2);
semDevSleepPreEventSpectrum = stdDevSleepPreEventSpectrum/sqrt(sizeSleepEvents);
clear preEventSpectrumSleepOnly
postEventSpectrumSleepOnly = epochPowerSpectrum(:,sleepIdx,2);
meanSleepPostEventSpectrum = mean(postEventSpectrumSleepOnly,2);
stdDevSleepPostEventSpectrum = std(postEventSpectrumSleepOnly,0,2);
semDevSleepPostEventSpectrum = stdDevSleepPostEventSpectrum/sqrt(sizeSleepEvents);
clear postEventSpectrumSleepOnly

% bandPowerSleep = allEvntsBandPwr(:, sleepIdx, :);
% 
% meanBandPowerSleep = mean(bandPowerSleep,1);
% stdBandPowerSleep = std(bandPowerSleep,1);
% semBandPowerSleep = stdBandPowerSleep/sqrt(sizeSleepEvents);


%Write header for PSD file before entering loop to calculate PSd values for
%all 2s epochs
c = clock;
dt = datestr(c,'mmddyy-HHMM');
resultExcelFileName = ['EventPSD_' fileLabel '_' dt 'v02.xlsx'];
warning off MATLAB:xlswrite:AddSheet
sheetName = 'Files used';
xlswrite(resultExcelFileName,cellstr('EEG'), sheetName, 'A2');
xlswrite(resultExcelFileName,cellstr(cscFile), sheetName, 'B2');

%Write PSD information to a new sheet in the Excel file.
sheetName = 'BandPower-PrePost';
%Generalize on next version:
xlswrite(resultExcelFileName,cellstr('DeltaPower'), sheetName, 'C1');
xlswrite(resultExcelFileName,cellstr('ThetaPower'), sheetName, 'E1');
xlswrite(resultExcelFileName,cellstr('SigmaPower'), sheetName, 'G1');
xlswrite(resultExcelFileName,cellstr('BetaPower'), sheetName, 'I1');

columnHeaders = {'Timestamp','State', 'Pre', 'Post','Pre', 'Post','Pre', 'Post','Pre', 'Post'};
xlswrite(resultExcelFileName,columnHeaders, sheetName, 'A2');
% xlswrite(resultExcelFileName,cellstr('Timestamp'), sheetName, 'A2');
% xlswrite(resultExcelFileName,cellstr('Pre'), sheetName, 'B2');
% xlswrite(resultExcelFileName,cellstr('Post'), sheetName, 'C2');
% xlswrite(resultExcelFileName,cellstr('Pre'), sheetName, 'D2');
% xlswrite(resultExcelFileName,cellstr('Post'), sheetName, 'E2');
% xlswrite(resultExcelFileName,cellstr('Pre'), sheetName, 'F2');
% xlswrite(resultExcelFileName,cellstr('Post'), sheetName, 'G2');
% xlswrite(resultExcelFileName,cellstr('Pre'), sheetName, 'H2');
% xlswrite(resultExcelFileName,cellstr('Post'), sheetName, 'I2');
prePostBndPwr = [squeeze(allEvntsBandPwr(1,:,:)), squeeze(allEvntsBandPwr(2,:,:)),...
    squeeze(allEvntsBandPwr(3,:,:)), squeeze(allEvntsBandPwr(4,:,:))];
xlswrite(resultExcelFileName,[eventTimeStamps labeledEventStates prePostBndPwr],sheetName, 'A3');

%Calculate band power averages for QS for entire record:
sleepPrePostBndPwr = prePostBndPwr(sleepIdx,:);
avgSleepPrePostBndPwr = mean(sleepPrePostBndPwr);
stdSleepPrePostBndPwr = std(sleepPrePostBndPwr);
semSleepPrePostBndPwr = stdSleepPrePostBndPwr/sqrt(sizeSleepEvents);

%Write band power averages for QS for entire record:
sheetName = 'QS-BndPwrAvg';
%Generalize on next version:
xlswrite(resultExcelFileName,cellstr('DeltaPower'), sheetName, 'B1');
xlswrite(resultExcelFileName,cellstr('ThetaPower'), sheetName, 'D1');
xlswrite(resultExcelFileName,cellstr('SigmaPower'), sheetName, 'F1');
xlswrite(resultExcelFileName,cellstr('BetaPower'), sheetName, 'H1');
columnHeaders = {'Pre', 'Post','Pre', 'Post','Pre', 'Post','Pre', 'Post'};
xlswrite(resultExcelFileName,columnHeaders, sheetName, 'B2');
xlswrite(resultExcelFileName,cellstr('Avg'), sheetName, 'A3');
xlswrite(resultExcelFileName,avgSleepPrePostBndPwr,sheetName, 'B3');
xlswrite(resultExcelFileName,cellstr('Std'), sheetName, 'A4');
xlswrite(resultExcelFileName,stdSleepPrePostBndPwr,sheetName, 'B4');
xlswrite(resultExcelFileName,cellstr('SEM'), sheetName, 'A5');
xlswrite(resultExcelFileName,semSleepPrePostBndPwr,sheetName, 'B5');
xlswrite(resultExcelFileName,cellstr('SampleSize'), sheetName, 'A6');
xlswrite(resultExcelFileName,sizeSleepEvents,sheetName, 'B6');

%Write full spectrum Pre-Event to a new sheet in the Excel file.
sheetName = 'PreEventSpectrum';
xlswrite(resultExcelFileName,cellstr('Timestamp'), sheetName, 'A1');
xlswrite(resultExcelFileName,cellstr('State'), sheetName, 'B1');
xlswrite(resultExcelFileName,freqVector', sheetName, 'C1');
xlswrite(resultExcelFileName,eventTimeStamps,sheetName, 'A2');
xlswrite(resultExcelFileName,labeledEventStates,sheetName, 'B2');
prePwrSpectra = squeeze(epochPowerSpectrum(:,:,1))';
xlswrite(resultExcelFileName,prePwrSpectra,sheetName, 'C2');
clear prePwrSpectra
%Write full spectrum Post-Event to a new sheet in the Excel file.
sheetName = 'PostEventSpectrum';
xlswrite(resultExcelFileName,cellstr('Timestamp'), sheetName, 'A1');
xlswrite(resultExcelFileName,cellstr('State'), sheetName, 'B1');
xlswrite(resultExcelFileName,freqVector', sheetName, 'C1');
xlswrite(resultExcelFileName,eventTimeStamps,sheetName, 'A2');
xlswrite(resultExcelFileName,labeledEventStates,sheetName, 'B2');
postPwrSpectra = squeeze(epochPowerSpectrum(:,:,2))';
xlswrite(resultExcelFileName,postPwrSpectra,sheetName, 'C2');
clear postPwrSpectra
% xlswrite(resultExcelFileName,cellstr('SleepMean'), sheetName, 'A3');
% xlswrite(resultExcelFileName,meanBandPowerSleep,sheetName, 'B3');
% xlswrite(resultExcelFileName,cellstr('SleepStd'), sheetName, 'A4');
% xlswrite(resultExcelFileName,stdBandPowerSleep,sheetName, 'B4');
% xlswrite(resultExcelFileName,cellstr('SleepSem'), sheetName, 'A5');
% xlswrite(resultExcelFileName,semBandPowerSleep,sheetName, 'B5');
% xlswrite(resultExcelFileName,[eventTimeStamps deltaPower thetaPower sigmaPower betaPower],sheetName, 'A6');

%Sleep Averages:
sheetName = 'SleepSpectrumAvg';
xlswrite(resultExcelFileName,cellstr('PreEventAvg'), sheetName, 'A2');
xlswrite(resultExcelFileName,meanSleepPreEventSpectrum',sheetName, 'B2');
xlswrite(resultExcelFileName,cellstr('PreEventStd'), sheetName, 'A3');
xlswrite(resultExcelFileName,stdDevSleepPreEventSpectrum',sheetName, 'B3');
xlswrite(resultExcelFileName,cellstr('PreEventSem'), sheetName, 'A4');
xlswrite(resultExcelFileName,semDevSleepPreEventSpectrum',sheetName, 'B4');
xlswrite(resultExcelFileName,cellstr('SampleSize'), sheetName, 'A5');
xlswrite(resultExcelFileName,sizeSleepEvents,sheetName, 'B5');
xlswrite(resultExcelFileName,cellstr('PostEventAvg'), sheetName, 'A6');
xlswrite(resultExcelFileName,meanSleepPostEventSpectrum',sheetName, 'B6');
xlswrite(resultExcelFileName,cellstr('PostEventStd'), sheetName, 'A7');
xlswrite(resultExcelFileName,stdDevSleepPostEventSpectrum',sheetName, 'B7');
xlswrite(resultExcelFileName,cellstr('PostEventSem'), sheetName, 'A8');
xlswrite(resultExcelFileName,semDevSleepPostEventSpectrum',sheetName, 'B8');

