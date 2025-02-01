% PIlots, Cyclostationarity, Statistical QAM, DC Offset, Guard bands, PAPR
% (Limit with correlation is -15-17 dB)
% Conclusion: Pilot Addition is definitely better as well. Higher peaks but more spread out.
clear;

%% Initial Parameters
numBits = 1000000;


%% Bit Stream
bitStream = randi([0 1],numBits,1);

%% Waveform Generation
cfgnonHT = wlanNonHTConfig('MCS', 3); 
waveform = wlanWaveformGenerator(bitStream, cfgnonHT);
nfft = 64;
cp = 16;
ofdmSymbolLength = nfft + cp;
pilotIndices = [-21, -7, 7, 21] + 33; % [1,1,1-1] for every OFDM symbol

txSignal = [zeros(80000,1); waveform; zeros(50000,1)];
txsignal_woNoise = txSignal;
sigpower = pow2db(mean(abs(waveform).^2));
txSignal = awgn(txSignal,-18,sigpower); % 2 dB gain only with Pilots

% Baseline xcorr
basePreamble = waveform(1:4*ofdmSymbolLength);
[r, idx] = xcorr(txSignal, basePreamble);

%% Start xcorr
rcorr = zeros(length(txSignal)-length(waveform),1);
timecorrArray = zeros(length(txSignal)-length(waveform),1);
freqCorrArray = zeros(length(txSignal)-length(waveform),1);
for curIdx = 1:length(txSignal)-length(waveform)
    curwindow = txSignal(curIdx:curIdx+length(waveform)-1);
    timeCorr = 0;
    freqCorr = 0;
    timeCorr = timeCorr + sum(curwindow(1:length(basePreamble)).*conj(basePreamble));
    curcurIdx = length(basePreamble)+1;
    symbolNum = 0;
    while curcurIdx < length(waveform)
        currSymbol = curwindow(curcurIdx:curcurIdx+ofdmSymbolLength-1);
        currSymbolfft = fftshift(fft(currSymbol(cp+1:end)));
        % currSymbolfft(abs(currSymbolfft)<0.1) = 0;
        pilots = currSymbolfft(pilotIndices);
        actualPilots = wlan.internal.nonHTPilots(1,symbolNum);
        freqCorr = freqCorr + sum(pilots.*conj(actualPilots));
        curcurIdx = curcurIdx + ofdmSymbolLength;
        symbolNum = symbolNum + 1;
    end
    timecorrArray(curIdx) = timeCorr;
    freqCorrArray(curIdx) = freqCorr;
end
rcorr = abs(timecorrArray) + abs(freqCorrArray);

subplot(2,2,1);
plot(idx, abs(r));
subplot(2,2,2);
plot(abs(timecorrArray));
subplot(2,2,3);
plot(abs(freqCorrArray));
subplot(2,2,4);
plot(abs(rcorr));

