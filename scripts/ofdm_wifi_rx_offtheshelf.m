clear;
clc;

sample_size = 1e6;
%% Getting the preamble
cfgnonHT = wlanNonHTConfig;
txWaveform = wlanWaveformGenerator(1, cfgnonHT);
preamble = txWaveform(1:320);

%% Getting signals from receiver
for i = 1:1000
   a = read_complex_binaryshort_chunk('wifi_actual.dat',sample_size,(i-1)*sample_size+1);
%    [r,idx] = xcorr(a, preamble);
%     [startoffset, M] = wlanPacketDetect(a,"CBW20",0,0.99);
end

