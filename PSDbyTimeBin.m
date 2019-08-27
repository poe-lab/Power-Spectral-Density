working_dir=pwd;
%Define frequency bands to calculate power within:
delta = [0.4 4];
theta = [5 9];
sigma = [10 14]; 
beta = [15 20]; 
freqBands = [delta; theta; sigma; beta];


% Select folder and get list of Excel files to fix:
fileType = '*.xlsx';
[dataFolder, fileName, numberOfDataFiles] = batchLoadFiles(fileType);
for m = 1:numberOfDataFiles
    spectraFile= fullfile(dataFolder, fileName(m,:));
    sheet = 'PowerSpectra';
    % load the dB version of the power spectra:
    % current_dir='C:\SleepData';
    % cd(current_dir);
    % [filename, pathname] = uigetfile({'*.xlsx','Excel File (*.xlsx)'},'Select the old spectral results file in dB');
    % if isequal(filename,0) || isequal(pathname,0)
    %     uiwait(errordlg('You need to select a file. Please press the button again',...
    %         'ERROR','modal'));
    %     cd(working_dir);
    % else
    %     cd(working_dir);
    %     spectraFile = fullfile(pathname, filename);
    % end
    num = xlsread(spectraFile,sheet);
    timeStamps = num(2:end,1);
    states = num(2:end,2);
    freqVector = num(1,3:end);
    dBpower = num(2:end,3:end);
    clear num
    %Note that the PSD values are in dB, so may need to convert back to normal
    %before taking averages. Each row is an epoch's PSD.
    epochPowerSpectra = db2pow(dBpower);
    clear dBpower
    %Write column headers for results to the 'fixedPower' sheet:
    sheetName = 'fixedPower';
    columnHeaders = {'Power_Spectra(uV^2/Hz)'};
    xlswrite(spectraFile,columnHeaders, sheetName, 'A1');
    clear columnHeaders
    columnHeaders = {'Frequency(Hz)'};
    xlswrite(spectraFile,columnHeaders, sheetName, 'C1');
    clear columnHeaders
    columnHeaders = {'Time_(s)','Scored_State'};
    xlswrite(spectraFile,columnHeaders, sheetName, 'A2');
    clear columnHeaders
    xlswrite(spectraFile,freqVector, sheetName, 'C2');
    %Write results to the 'PowerSpectra' sheet:
    xlswrite(spectraFile,[timeStamps states epochPowerSpectra], sheetName, 'A3');
    epochPowerSpectra = epochPowerSpectra'; %Force into column vector
    freqVector = freqVector(:); %Force into column vector
    numEpochs = size(epochPowerSpectra,2);
    numBands = size(freqBands,1);
    allEpochsBandpwr = zeros(numBands, numEpochs);
    width = diff(freqVector);
    lastRectWidth = 0;  % Don't include last point of PSD data.
    width = [width; lastRectWidth];
    for i = 1:numBands
        idx1 = find(freqVector<=freqBands(i,1), 1, 'last' );
        idx2 = find(freqVector>=freqBands(i,2), 1, 'first');
        allEpochsBandpwr(i,:) = width(idx1:idx2)'*epochPowerSpectra(idx1:idx2,:);
        clear idx1 idx2
    end
    clear epochPowerSpectra width freqVector
    allEpochsBandpwr = allEpochsBandpwr';
    %Write column headers for results to the 'bandPower' sheet:
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
    %Write results to the 'PowerSpectra' sheet:
    xlswrite(spectraFile,[timeStamps states allEpochsBandpwr], sheetName, 'A3');
    clear timeStamps states allEpochsBandpwr
end
cd(working_dir);
msgbox('Fix complete.','Pop-up');

%Now need to find averages per time bin for each state

