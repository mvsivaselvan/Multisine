function [t, wmsin] = WindowedMsinGen(numFFTpoints, nFreq, ...
                                      nWindows, nRampSamples, scale, ...
                                      sampleRate, csvFileName)

% numFFTpoints = number of samples in the msine signals
% nFreq = number of frequency components (so bandwidth =
%         nFreq*(1/SampleRate))
% nWindows = number of repetitions of the msin
% nRampSamples = number of samples over which to ramp up and down
% scale = scale factor for msin signal
% sampleRate = sampling frequency (Hz)
% csvFileName = CSV filename to write data in to use with Veristand
%               Stimulus Profile

msin = msingen(numFFTpoints, nFreq);

wmsin = repmat(msin(:), nWindows, 1);

ramp = ones(size(wmsin));

ramp(1:nRampSamples) = (1-cos((0:(nRampSamples-1))/(nRampSamples-1)*pi))/2;
ramp(end:-1:(end-nRampSamples+1)) = ramp(1:nRampSamples);

wmsin = (wmsin.*ramp)*scale;
t = (0:length(wmsin)-1)'/sampleRate;

tab = table(t*1000, wmsin, 'VariableNames', {'timestamp', 'Multisine'});
    % t*1000 to convert to milli-s
writetable(tab, csvFileName);
