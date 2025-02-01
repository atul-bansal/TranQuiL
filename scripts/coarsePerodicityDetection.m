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

figure;
hold on;
sliding_window_sd_snrVariation = [];
for snrvalues = 10:-5:-30
    signalWithNoise = awgn(signal, snrvalues, 'measured');
    signalPower = (abs(signalWithNoise)).^2;
    
    lowerWindowsize = periodSignalSamples-6000;
    % lowerLog = log10(lowerWindowsize/2);
    higherWindowsize = periodSignalSamples+6000;
    % higherLog = log10(higherWindowsize/2);
    
    sliding_window_sd_signNoise = [];
    for actualWindowSize = lowerWindowsize:2000:higherWindowsize
        actualWindowSize
    %     actualWindowSize = floor(2*10^logPeriods);
        windowSum = zeros(length(signalWithNoise)-actualWindowSize+1,1);
        for i = 1:length(windowSum)
            windowSum(i) = mean(signalPower(i:i+actualWindowSize-1));
        end
        sliding_window_sd_signNoise = [sliding_window_sd_signNoise; sqrt(var(windowSum))];
    %     figure;
    %     plot(windowSum)
    %     title(num2str(actualWindowSize))
    end
    sliding_window_sd_snrVariation = [sliding_window_sd_snrVariation sliding_window_sd_signNoise];
    plot(lowerWindowsize:2000:higherWindowsize, sliding_window_sd_signNoise)
end

%% Separate Processing For better visualization
sliding_window_sd_snrVariation = sliding_window_sd_snrVariation(:,1:7);
sliding_window_sd_snrVariation(:,end) = sliding_window_sd_snrVariation(:,end) - 0.1;
sliding_window_sd_snrVariation(:,end-1) = sliding_window_sd_snrVariation(:,end-1) - 0.05;
snrValues = 10:-5:-20;

figure; 
plot(lowerWindowsize:2000:higherWindowsize, sliding_window_sd_snrVariation)
