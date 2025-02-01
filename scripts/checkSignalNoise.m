clc;
clear;

period = 10e-3;
fs = 5e6;
periodSignalSamples = period * fs;
packetLength = 20000;
start_sample = 10000;

onePeriodSignal = zeros(periodSignalSamples, 1);
onePeriodSignal(start_sample + 1:start_sample + packetLength) = (1+1i)/sqrt(2);
repetitions = 10;

signal = repmat(onePeriodSignal, repetitions, 1);

signalWithNoise = awgn(signal, -30, 'measured');
noise = signalWithNoise - signal;
noise_sq = noise.^2;
figure;
plot(abs(signalWithNoise))

signalPower = (abs(signalWithNoise)).^2;

lowerWindowsize = 48000;
% lowerLog = log10(lowerWindowsize/2);
higherWindowsize = periodSignalSamples+3000;
% higherLog = log10(higherWindowsize/2);

sliding_window_sd_signNoise = [];
sliding_window_sd_Noise = [];

for actualWindowSize = lowerWindowsize:1000:higherWindowsize
    actualWindowSize
%     actualWindowSize = floor(2*10^logPeriods);
    windowSum = zeros(length(signalWithNoise)-actualWindowSize+1,1);
    windowSum_Noise = windowSum;
    for i = 1:length(windowSum)
        windowSum(i) = mean(signalPower(i:i+actualWindowSize-1));
        windowSum_Noise(i) = mean(noise_sq(i:i+actualWindowSize-1));
    end
    sliding_window_sd_signNoise = [sliding_window_sd_signNoise; sqrt(var(windowSum))];
    sliding_window_sd_Noise  = [sliding_window_sd_Noise; sqrt(var(windowSum_Noise))];   
%     figure;
%     plot(windowSum)
%     title(num2str(actualWindowSize))
end

figure;
hold on;
plot(lowerWindowsize:1000:higherWindowsize, sliding_window_sd_signNoise)
plot(lowerWindowsize:1000:higherWindowsize, sliding_window_sd_Noise)
legend('Signal with noise', 'Noise');
title('SD variation with window size')




