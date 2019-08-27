function normalizedData = normalizePower(sleepStates, power)
%find all unhooked & not scored epochs to exclude from calculation of mean
row2Exclude = find(sleepStates==5|sleepStates==7);
incSpectra = power;
incSpectra(row2Exclude,:) = [];
meanSpectra = mean(incSpectra);
normalizedData = 0*power;
numBands = size(normalizedData,2);
for i=1:numBands
normalizedData(:,i) = power(:,i)./meanSpectra(i);
end