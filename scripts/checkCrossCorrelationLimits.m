clc;
clear;

rng(1,"twister");

%% Initial Parameters
numBits = 1000000;


%% Bit Stream
bitStream = randi([0 1],numBits,1);

%% Waveform Generation
cfgnonHT = wlanNonHTConfig('MCS', 3); 
waveform = wlanWaveformGenerator(bitStream, cfgnonHT);

preambleSize_stfLtf = 320;
preamble = waveform(1:preambleSize_stfLtf);

snrRange = -10:-1:-20;
startPads = 80000;
endPads = 50000;

txSignal = [zeros(startPads, 1); waveform; zeros(endPads, 1)];

for idx = 1:length(snrRange)
    sigpower = pow2db(mean(abs(waveform).^2));
    rxsignal = awgn(txSignal, snrRange(idx), sigpower);
    [r, indices] = xcorr(rxsignal, preamble);
    subplot(2,1,2); plot(abs(rxsignal));
    subplot(2,1,1); plot(indices, abs(r));
    title(['snr = ', num2str(snrRange(idx))]);
    close all;
end