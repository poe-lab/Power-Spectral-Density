function calculateBandPower_01282019
%% This program calculates the band power for the defined frequency ranges.
% Change the bandwidth values below to suit your research needs.

% ***WARNING: This program is not smart enough to know what frequency range
% was used to calculate the power spectra. Therefore, it will give 
% erroneous values if the frequency range is missing!!!!!!!!!!!!!!!!!!!!!!!

%Define frequency bands to calculate power within:
delta = [0.4 4.99]; % DEFAULT: [0.4 4]
theta = [5 9.99]; % DEFAULT: [5 9]
sigma = [10 14.99]; % DEFAULT: [10 14]
beta = [15 20]; % DEFAULT: [15 20]
lowGamma = [30 50];  % DEFAULT: [30 50]
highGamma = [80 100]; % DEFAULT: [80 100]

%% Include frequency bands to analyze in the following array:
freqBands = [delta; theta; sigma; beta; lowGamma; highGamma];

%% Select folder and get list of Excel files to fix:
working_dir=pwd;

[fileSet, filePath] = uigetfile('*.xlsx','Select power spectra file(s)','MultiSelect', 'on');

numFiles = size(fileSet,2); % # of files selected
if iscell(fileSet)    
else
    numFiles = 1;
end


%% Run band power analyses for each file:
for p = 1:numFiles
    %% Load the power spectra file:
    if iscell(fileSet)
        spectraFileName = fileSet{1,p};
    else
        spectraFileName = fileSet;
    end
    spectraFile = fullfile(filePath, spectraFileName);
    sheet = 'PowerSpectra'; %name of sheet to run function on
    num = xlsread(spectraFile,sheet);
    timeStamps = num(2:end,1); %from second row till the end in column 1
    states = num(2:end,2); %from second row till the end in column 2
    freqVector = num(1,3:end);
    epochPowerSpectra = num(2:end,3:end);
    clear num

    epochPowerSpectra = epochPowerSpectra'; %Force into column vector
    freqVector = freqVector(:); %Force into column vector
    numEpochs = size(epochPowerSpectra,2);
    numBands = size(freqBands,1);
    allEpochsBandpwr = zeros(numBands, numEpochs);
    widthTemp = diff(freqVector);
    lastRectWidth = 0;  % Don't include last point of PSD data.
    width = [widthTemp; lastRectWidth];
    targFreq = ones(numBands,1); % Marks which bands to keep for output
    
    %% Calculate band power of each band in the spectral data:
    for i = 1:numBands
        idx1 = find(freqVector<=freqBands(i,1), 1, 'last' );
        % Makes sure start of freq band is in spectra data
        if isempty(idx1) 
            % Remove the frequency band:
            targFeq(i) = 0; %#ok<*AGROW,NASGU>
        else % If the start freq of band is in spectra data
            idx2 = find(freqVector>=freqBands(i,2), 1, 'first');
            % Check if end of freq band is in spectra data
            if isempty(idx2)
                % If end of freq band is not in the data, replace it with
                % last frequency value of spectral data.
                idx2 = length(freqVector);
                freqBands(i,2) = freqVector(end);
            end
            allEpochsBandpwr(i,:) = width(idx1:idx2)'*epochPowerSpectra(idx1:idx2,:);
        end
        
        clear idx1 idx2
    end
    clear epochPowerSpectra width freqVector
    targLogic = targFreq == 1;
    freqBands = freqBands(targLogic, :);
    allEpochsBandpwr = allEpochsBandpwr(targLogic, :);
    allEpochsBandpwr = allEpochsBandpwr';
    
    %% Write column headers for results to the 'bandPower' sheet:
    sheetName = 'bandPower';
    columnHeaders = {'Power(uV^2)'};
    xlswrite(spectraFile,columnHeaders, sheetName, 'A1');
    clear columnHeaders
    columnHeaders = {'Bandwidths(Hz)'};
    xlswrite(spectraFile,columnHeaders, sheetName, 'B1');
    clear columnHeaders
    columnHeaders = {'Time_(s)','Scored_State'};
    xlswrite(spectraFile,columnHeaders, sheetName, 'A2');
    clear columnHeaders
    xlswrite(spectraFile,freqBands', sheetName, 'C1');
    %% Write results to the 'PowerSpectra' sheet:
    xlswrite(spectraFile,[timeStamps states allEpochsBandpwr], sheetName, 'A3');
    clear timeStamps states allEpochsBandpwr
end
cd(working_dir);
msgbox('Band power calculations complete.','Pop-up');

