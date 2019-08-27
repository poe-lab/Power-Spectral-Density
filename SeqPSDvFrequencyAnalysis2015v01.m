function varargout = SeqPSDvFrequencyAnalysis2015v01(varargin)
% SEQPSDVFREQUENCYANALYSIS2015V01 Application M-file for SeqPSDvFrequencyAnalysis2015v01.fig
%    FIG = SEQPSDVFREQUENCYANALYSIS2015V01 launch SeqPSDvFrequencyAnalysis2015v01 GUI.
%    SEQPSDVFREQUENCYANALYSIS2015V01('callback_unhooked', ...) invoke the named callback.
% Last Modified by GUIDE v2.5 13-Apr-2015 09:23:30
% 
%04022015: changed to convert only REAL part of Sxx to decibels so that
% warning no longer appears. -BAG
%04072015: Removed conversion to dB of power since addition of dB values is
%not the same as addition of power values (i.e., log A +log B = log (A x B)
if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse',varargin{:});
    
    % Generate a structure of handles to pass to callbacks, and store it.  
    handles = guihandles(fig);
    guidata(fig, handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr); %#ok<LERR>
    end
    
end
%##########################################################################
function LoadFiles_pushbutton_Callback(hObject, eventdata, handles)
global EEG_Fc EEG_HP_Fc
EEG_Fc = 60; set(handles.EEG_cutoff, 'String', EEG_Fc);
EEG_HP_Fc = 0.4; set(handles.EEG_HP_cutoff, 'String', EEG_HP_Fc);
% Call the time stamp file name:
working_dir=pwd;
current_dir='C:\SleepData';
cd(current_dir);
[filename, pathname] = uigetfile('*.xls', 'Pick the timestamp file for these datafiles');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    timestampfile= fullfile(pathname, filename);
    set(handles.tstampsfile,'string',filename);
    set(handles.tstampsfile,'Tooltipstring',timestampfile);
end

% Call the Cortical CSC file name:
current_dir='C:\SleepData';
cd(current_dir);
[filename, pathname] = uigetfile({'*.ncs','Neuralynx CSC File (*.ncs)'},'Select the EEG');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    cortexCscFile= fullfile(pathname, filename);
    set(handles.cortexCSC,'string',filename);
    set(handles.cortexCSC,'Tooltipstring',cortexCscFile);
end

% Call the sleep scored file name:
current_dir='C:\SleepData';
cd(current_dir);
[filename, pathname] = uigetfile({'*.xls','Excel 1997-2003 File (*.xls)'},'Select the Stage Scored File');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    cortexScoredFile = fullfile(pathname, filename);
    set(handles.cortexScored,'string',filename);
    set(handles.cortexScored,'Tooltipstring',cortexScoredFile);
end


function runSpectralAnalysis_Callback(hObject, eventdata, handles)
global EEG_Fc EEG_HP_Fc
handles=guihandles(SeqPSDvFrequencyAnalysis2015v01);
% This file uses PERSISTENT variables, and hence we clear them on every new use

% Get UI file names
sleepScoredFile=get(handles.cortexScored,'TooltipString');
CscFilename=get(handles.cortexCSC,'TooltipString');
timeStampFile=get(handles.tstampsfile,'TooltipString');

% Read in the timestampsfilebutton file:
try
    timestampSections = xlsread(timeStampFile);
catch
    uiwait(errordlg('Check if the file is saved in Microsoft Excel format.',...
        'ERROR','modal'));
end

%% LOAD SLEEP SCORED DATA
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
nEpochPoints=length(scoredStates);
%% LOAD EEG DATA
[EEG_samples, EEG_timestamps, reducedSamplingRate] = loadEEGData(CscFilename, timestampSections, EEG_Fc, EEG_HP_Fc);

% Find the EPOCHSIZE of 10 seconds:
index=find((EEG_timestamps(1)+9.999 < EEG_timestamps) & (EEG_timestamps < EEG_timestamps(1)+10.001));
if (isempty(index)) == 1
    index=find((EEG_timestamps(1)+9.99 < EEG_timestamps) & (EEG_timestamps < EEG_timestamps(1)+10.01));
end
diff= EEG_timestamps(index(1):index(end)) - (EEG_timestamps(1)+10);
[minimum,ind]=min(abs(diff)); %#ok<ASGLU>
try
    EPOCHSIZE=index(ind);
catch %#ok<CTCH>
    fprintf('There is an error in calculating the EPOCHSIZE of 10sec in read_n_extract_datafiles\n');
end
epochDuration = EEG_timestamps(EPOCHSIZE+1) - EEG_timestamps(1);
dt = 1/reducedSamplingRate; % Define the sampling interval.
df = 1/epochDuration; % Determine the frequency resolution.
fNyquist = reducedSamplingRate/2; % Determine the Nyquist frequency.

%%
faxis = (0:df:fNyquist); % Construct the frequency axis.
lastFreqIndex = find(faxis <= EEG_Fc, 1, 'last');
faxis = faxis(1:lastFreqIndex);
scoredEpochPowerSpectrum = zeros(nEpochPoints, lastFreqIndex);

for i = 1:nEpochPoints
    stPointIndex = find((EEG_timestamps > (scoredStates(i,1) - 0.1 )) & (EEG_timestamps <scoredStates(i,1) + 0.1));
    if isempty(stPointIndex)
        
    else
        diff2 = EEG_timestamps(stPointIndex) - scoredStates(i,1);
        [~,ind2]=min(abs(diff2));
        st_pt = stPointIndex(ind2);
        end_pt=st_pt+EPOCHSIZE-1;
        if end_pt <= length(EEG_timestamps)
        % Compute the power spectrum of the Hann tapered data:
            x = EEG_samples(st_pt:end_pt);
            xh = hann(length(x)).*x; 
            Sxx = 2*dt^2/epochDuration * fft(xh).*conj(fft(xh)); % Compute the power spectrum of Hann tapered data.
            Sxx=real(Sxx); %Ignores imaginary component.
            %Should not convert to dB power before calculating band power and/or taking
%averages:
%             Sxx=10*log10(real(Sxx)); % ... convert to a decibel scale. Ignores imaginary component.
%             Sxx = Sxx(1:length(x)/2+1); % ... ignore negative frequencies. Only needed if not limiting frequency range like next line:		
            scoredEpochPowerSpectrum(i, :) = Sxx(1:lastFreqIndex);
        end
    end
end

meanSpectrumState = zeros(8, lastFreqIndex);
stdDevSpectrumState = zeros(8, lastFreqIndex);
sampleSizeSpectrumState = zeros(8,1);

for i = 1:8
    if i == 7
    else
        stateIndex = find(scoredStates(:,2) == i); %#ok<*AGROW>
        if isempty(stateIndex)==0
            subsetEpochPowerSpectrum = scoredEpochPowerSpectrum(stateIndex, :);    % Separate spectrums for epochs by scored state
            meanSpectrumState(i,:) = mean(subsetEpochPowerSpectrum); % Average spectrums for each state.
            stdDevSpectrumState(i,:) = std(subsetEpochPowerSpectrum);
            sampleSizeSpectrumState(i) = size(subsetEpochPowerSpectrum,1);
            % Generate figures for average spectrums of each state.
            figure(i)
            plot(faxis,meanSpectrumState(i,:)); % Plot power versus frequency.
            %set(gca,'XLim',[0 EEG_Fc])
%             ylim([-80 10]) % Set range of y-values in plot.
            xlabel('Frequency [Hz]');  ylabel('Power [uV^2/Hz]') % Label the axes.
            switch i
                case 1
                    title('ACTIVE WAKE: Averaged Spectrum over All Epochs')
                case 2
                    title('QUIET SLEEP: Averaged Spectrum over All Epochs')
                case 3
                    title('REM: Averaged Spectrum over All Epochs')
                case 4
                    title('QUIET WAKE: Averaged Spectrum over All Epochs')
                case 5
                    title('UNHOOKED: Averaged Spectrum over All Epochs')
                case 6
                    title('TRANSITION to REM: Averaged Spectrum over All Epochs')
                case 7
                    title('NOT SCORED: Averaged Spectrum over All Epochs')
                case 8
                    title('INTERMEDIATE WAKE: Averaged Spectrum over All Epochs')
            end
        end
    end 
end
% Write all Results to an Excel file:

%Request user input to name the file:
prompt = {'File name for spectral analysis results:'};
def = {'SubjectNumber'};
dlgTitle = 'Save Results';
lineNo = 1;
answer = inputdlg(prompt,dlgTitle,lineNo,def);
filenameUserInput = char(answer(1,:));
dateString = date;
timeSaved = clock;
timeSavedString = [num2str(timeSaved(4)) '-' num2str(timeSaved(5))];
resultsFilename = strcat('C:\SleepData\', filenameUserInput, dateString, '_', timeSavedString, '.xlsx');
clear timeSaved timeSavedString filenameUserInput

%Request user input to name the file:
prompt = {'Enter user name:'};
def = {'Your name'};
dlgTitle = 'User Name';
lineNo = 1;
answer = inputdlg(prompt,dlgTitle,lineNo,def);
username = char(answer(1,:));
warning off MATLAB:xlswrite:AddSheet

%Write all file names used for analyses to the 'headerInfo' sheet:
headerArray = {'Date_of_Analysis', dateString; 'User_Name', username;...
    'CSC_File', CscFilename; 'Scored_File', sleepScoredFile};
sheetName = 'HeaderInfo';
xlswrite(resultsFilename,headerArray, sheetName);
clear headerArray username prompt def lineNo answer tetrodeFileName dateString

%Write frequency settings used for analyses to the 'headerInfo' sheet:
headerArray = {'Low_Pass_Hz', EEG_Fc};
xlswrite(resultsFilename,headerArray, sheetName, 'A7');
clear headerArray sheetName

%Write column headers for results to the 'PowerSpectra' sheet:
sheetName = 'PowerSpectra';
columnHeaders = {'Power_Spectra(uV^2/Hz)'};
xlswrite(resultsFilename,columnHeaders, sheetName, 'A1');
clear columnHeaders
columnHeaders = {'Frequency(Hz)'};
xlswrite(resultsFilename,columnHeaders, sheetName, 'C1');
clear columnHeaders
columnHeaders = {'Time_(s)','Scored_State'};
xlswrite(resultsFilename,columnHeaders, sheetName, 'A2');
clear columnHeaders
xlswrite(resultsFilename,faxis, sheetName, 'C2');
%Write results to the 'PowerSpectra' sheet:
xlswrite(resultsFilename,[scoredStates scoredEpochPowerSpectrum], sheetName, 'A3');

%Write column headers for results to the 'StatsByState' sheet:
sheetName = 'StatsByState';
columnHeaders = {'MEAN_POWER_OF_State'; 'State'};
columnHeader2 = {'AW'; 'QS'; 'RE'; 'QW'; 'UH'; 'TR'; 'NS'; 'IW'};
xlswrite(resultsFilename,columnHeaders, sheetName, 'A1');
clear columnHeaders
xlswrite(resultsFilename,faxis, sheetName, 'B2');
xlswrite(resultsFilename,columnHeader2, sheetName, 'A3');

%Write means to the 'StatsByState' sheet:
xlswrite(resultsFilename, meanSpectrumState, sheetName, 'B3');

columnHeaders = {'STD_DEV'; 'State'};
xlswrite(resultsFilename,columnHeaders, sheetName, 'A15');
clear columnHeaders
xlswrite(resultsFilename,faxis, sheetName, 'B16');
xlswrite(resultsFilename,columnHeader2, sheetName, 'A17');

%Write standard deviation results to the 'StatsByState' sheet:
xlswrite(resultsFilename, stdDevSpectrumState, sheetName, 'B17');

columnHeaders = {'SAMPLE_SIZE'; 'State'};
xlswrite(resultsFilename,columnHeaders, sheetName, 'A30');
clear columnHeaders
columnHeaders = {'SampleSize'};
xlswrite(resultsFilename,columnHeaders, sheetName, 'B31');
clear columnHeaders
xlswrite(resultsFilename,columnHeader2, sheetName, 'A32');

%Write standard deviation results to the 'StatsByState' sheet:
xlswrite(resultsFilename, sampleSizeSpectrumState, sheetName, 'B32');

function [EEG_samples, EEG_timestamps, Fs] = loadEEGData(CscFilename, timestampSections, EEG_Fc, EEG_HP_Fc)
% This function loads in the CSC amplitude and timestamps for the selected
% CSC file.  It also automatically down-samples.

waithandle= waitbar(0.2,'Loading the EEG data');pause(0.2);
lowertimestamp = timestampSections(1,1);
uppertimestamp = timestampSections(end,2);
[Timestamps,cscSamplingRate,Samples]=Nlx2MatCSC(CscFilename,[1 0 1 0 1],0,4,[lowertimestamp uppertimestamp]);
cscSamplingRate = cscSamplingRate(1);
unfilt_samples=double(Samples(:)');
exactLow = timestampSections(1,3);
exactHi = length(unfilt_samples);
clear Samples
close(waithandle);

% Precise time stamps should be calculated here:
waithandle = waitbar(0.6,'Extracting EEG Timestamps...'); pause(0.2);
[EEG_TIMESTAMPS_temp,EEG_SAMPLES_temp] = generate_timestamps_from_Ncsfiles(Timestamps,unfilt_samples,exactLow, exactHi,[]);
clear Timestamps unfilt_samples
close(waithandle);

% Automated down-sampling to down-smple to 1kHz if greater or output an error if below 250Hz:
if cscSamplingRate > 1050
    DS = (1:1:32);
    DSampSF = cscSamplingRate./DS;
    indSampfactor = find(DSampSF >= 1000);
    Fs = DSampSF(indSampfactor(end));
    sampFactor = DS(indSampfactor(end));
    msgbox({['Recording Sampling Rate:  ' num2str(cscSamplingRate) 'Hz'];...
        ['Down-Sampled Sampling Rate:  ' num2str(round(Fs)) 'Hz'];...
        ['Sampling Factor:  ' num2str(sampFactor) '']});
elseif cscSamplingRate < 250
    msgbox({['Recording Sampling Rate:  ' num2str(cscSamplingRate) 'Hz'];...
        'This sampling rate is too low to run the phase analyses.';...
        'Phase-O-Matic will now close.'});
    exit
else  % In practice, we should be sampling at at least 1kHz, but older files may have been recorded at slower frequencies.
    Fs = cscSamplingRate;
    sampFactor = 1;
    msgbox({['Recording Sampling Rate:  ' num2str(cscSamplingRate) 'Hz'];...
        'This sampling rate is acceptable to run the phase analyses.'});
end

% Convert amplitude to uV
physInput = 2;  %Needed to select proper error box in HeaderADBit.
ADBit2uV = HeaderADBit(CscFilename, physInput);    %Extract the AD Bit Value.
EEG_SAMPLES_temp = EEG_SAMPLES_temp * ADBit2uV;   %Convert amplitude of signal from AD Bits to microvolts.

%Design bandpass filter:
[z, p, k] = ellip(7,1,60, [EEG_HP_Fc EEG_Fc]/(cscSamplingRate/2),'bandpass');
[sos, g] = zp2sos(z,p,k);
%Apply bandpass filter to EEG data:
EEG_SAMPLES_temp = filtfilt(sos,g, EEG_SAMPLES_temp);
% %  Low pass filter for EEG signals
% [Blow,Alow] = ellip(7,1,60, EEG_Fc/(cscSamplingRate/2));           % Default setting implements low pass filter with 30hz cutoff
% EEG_SAMPLES_temp = filter(Blow,Alow,EEG_SAMPLES_temp);

waithandle = waitbar(0.3,'Downsampling EEG Data...'); pause(0.2);
EEG_TIMESTAMPS_temp = EEG_TIMESTAMPS_temp(1:sampFactor:end);
EEG_SAMPLES_temp = EEG_SAMPLES_temp(1:sampFactor:end);
close(waithandle);

%waithandle = waitbar(0.7,'Deshifting EEG Data...'); pause(0.2);
% EEG_timestamps = [];
% EEG_samples = [];
EEG_samples = EEG_SAMPLES_temp;%( floor( ( length(coeffs)/2 )/sampFactor ):(length(EEG_SAMPLES_temp) - floor( (length(coeffs)/2)/sampFactor) ));
EEG_timestamps = EEG_TIMESTAMPS_temp;%(1:length(EEG_samples));
% close(waithandle);
clear EEG_SAMPLES_temp EEG_TIMESTAMPS_temp

%##########################################################################
%                       User Input for EEG Filters
%##########################################################################
function EEG_cutoff_Callback(hObject, eventdata, handles)
global EEG_Fc
EEG_cutoff = str2double(get(hObject,'String'));   % returns contents of EEG_cutoff as a double
if isnan(EEG_cutoff)
    errordlg('Input must be a number','Error');
end
if EEG_cutoff > 125
    errordlg('Cutoff frequency must be < 125 Hz','Error');
end
EEG_Fc = EEG_cutoff;

function EEG_cutoff_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function cortexScored_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to cortexScored (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cortexScored as text
%        str2double(get(hObject,'String')) returns contents of cortexScored as a double


% --- Executes during object creation, after setting all properties.
function cortexScored_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cortexScored (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EEG_HP_cutoff_Callback(hObject, eventdata, handles)
global EEG_HP_Fc
EEG_HP_cutoff = str2double(get(hObject,'String'));   % returns contents of EEG_HP_cutoff as a double
if isnan(EEG_HP_cutoff)
    errordlg('Input must be a number','Error');
end
if EEG_HP_cutoff < 0.001
    errordlg('Cutoff frequency must be > 0.001 Hz','Error');
end
EEG_HP_Fc = EEG_HP_cutoff;


% --- Executes during object creation, after setting all properties.
function EEG_HP_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EEG_HP_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
