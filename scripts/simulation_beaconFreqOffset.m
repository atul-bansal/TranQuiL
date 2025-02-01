clear;

% Creating beacon frame
useSDR = false;
saveToFile = false;

ssid = "TEST_BEACON";
beaconInterval = 100;
band = 5;
chNum = 52;

frameBodyConfig = wlanMACManagementConfig( ...
    BeaconInterval=beaconInterval, ...
    SSID=ssid,...
    Timestamp=randi(2^53-1));

beaconFrameConfig = wlanMACFrameConfig(FrameType="Beacon", ...
    ManagementConfig=frameBodyConfig);

[mpduBits,mpduLength] = wlanMACFrame(beaconFrameConfig, OutputFormat="bits");
% [mpduChars,mpduLength] = wlanMACFrame(beaconFrameConfig);

% Creating beacon packet
fc = wlanChannelFrequency(chNum,band);

nonHTcfg = wlanNonHTConfig;       % Create packet configuration
% nonHTcfg.MCS = 6;                 % MCS = 6: (Modulation: 64QAM Rate: 2/3)
nonHTcfg.MCS = 0;                 % MCS = 0: (Modulation: BPSK Rate: 1/2)
nonHTcfg.NumTransmitAntennas = 1; % Number of transmit antenna
nonHTcfg.PSDULength = mpduLength; 
txWaveform = wlanWaveformGenerator(mpduBits,nonHTcfg);
packetSize = length(txWaveform);

beaconPreamble = txWaveform(1:944); % 944 number by - 400 Physical PReamble samples + Some MAC header samples

%% TimeStamp + Sequence Control + Source MAC address Variation
timestamp = randi(2^53-1);
seqValue = randi(4095);
sourceMAC = '60B76E4363FC';
txWaveform = [];
% deltaf = 125e3;
snrdB = -34;

for timeValue = 1:15
frameBodyConfig = wlanMACManagementConfig( ...
    BeaconInterval=beaconInterval, ...
    SSID=ssid,...
    Timestamp = timestamp+timeValue);

dsElementID = 3;
dsInformation = dec2hex(chNum,2);
frameBodyConfig = frameBodyConfig.addIE(dsElementID,dsInformation);

beaconFrameConfig = wlanMACFrameConfig(FrameType="Beacon", ...
    ManagementConfig=frameBodyConfig, SequenceNumber=seqValue+timeValue, Address2 = sourceMAC, Address3 = sourceMAC,...
    Address4 = sourceMAC);

[mpduBits,mpduLength] = wlanMACFrame(beaconFrameConfig, OutputFormat="bits");

fc = wlanChannelFrequency(chNum,band);

nonHTcfg = wlanNonHTConfig;       % Create packet configuration
% nonHTcfg.MCS = 6;               % MCS = 6: (Modulation: 64QAM Rate: 2/3)
nonHTcfg.MCS = 0;                 % MCS = 0: (Modulation: BPSK Rate: 1/2)
nonHTcfg.NumTransmitAntennas = 1; % Number of transmit antenna
nonHTcfg.PSDULength = mpduLength; 
tbtt = beaconInterval*1024e-6;
txWaveform = [txWaveform; wlanWaveformGenerator(mpduBits,nonHTcfg,Idletime=tbtt)];
end
% txWaveform = txWaveform .* exp(2i*pi*deltaf*(1:length(txWaveform)).'*(1/20e6))* exp(1i*2*pi*rand);
rxWaveform = awgn(txWaveform, snrdB, 'measured');

figure
r = xcorr(rxWaveform, beaconPreamble);
plot(abs(r));
title('Timestamp + Sequence Number + Source MAC variation + Frequency Offset');

%% Using Periodicity by adding packets
r_new = r(length(rxWaveform)-mod(length(rxWaveform),10000):end);
periodicSamples = 20e6*102.4e-3+packetSize;
r_new = r_new(1:floor(length(r_new)/periodicSamples)*periodicSamples);
r_new_matrix = reshape(r_new, periodicSamples, []);

dtwMatrix = Inf*ones(size(r_new_matrix,2));
indices = cell(size(r_new_matrix,2), size(r_new_matrix,2));

for i = 1:size(r_new_matrix,2)
    for j = i+1:size(r_new_matrix,2)
%         figure;
        [dtwMatrix(i,j), ix, iy] = dtw(downsample(r_new_matrix(:,i),64), downsample(r_new_matrix(:,j),64));
        indices{i,j} = {ix, iy};
%         title(['i = ', num2str(i),' j = ',num2str(j)]);
%         pause
%         close all;
    end
end

%% Further Processing
r_new_downsampled_abs = [];
for i = 1:size(r_new_matrix,2)
    r_new_downsampled_abs(:,i) = downsample(r_new_matrix(:,i),64);
end
[R, C] = ndgrid(1:size(dtwMatrix,1), 1:size(dtwMatrix,2));
[val, idx] = sort(dtwMatrix(:));
nonInfnum = size(r_new_matrix,2)*(size(r_new_matrix,2)-1)/2;
idx = idx(1:nonInfnum);
lengthWrap = zeros(size(r_new_matrix,2),1);
X1 = r_new_downsampled_abs(:,R(idx(1)));
X2 = r_new_downsampled_abs(:,C(idx(1)));
cellEntry = indices(R(idx(1)), C(idx(1)));
ix1 = cellEntry{1}(1);
ix2 = cellEntry{1}(2);
ix1 = ix1{1};
ix2 = ix2{1};
dtwX1 = X1(ix1);
dtwX2 = X2(ix2);
combination = dtwX1 + dtwX2;
original = dtwX2;
figure;
hold on;
plot(abs(dtwX1));
plot(abs(combination));
xlim([1, 300]);
pause
close all
globalLength = length(ix1);
lengthWrap(R(idx(1))) = globalLength;
lengthWrap(C(idx(1))) = globalLength;
i = 2;

while ~all(lengthWrap)
    if R(idx(i)) >= C(idx(i))
        i = i + 1;
        if i > length(idx)
            i = 1;
        end
        continue;
    end
    if lengthWrap(R(idx(i))) > 0 && lengthWrap(C(idx(i))) > 0
        i = i + 1;
        if i > length(idx)
            i = 1;
        end
        continue;
    elseif lengthWrap(R(idx(i))) > 0 
        X = r_new_downsampled_abs(:,C(idx(i)));
        cellEntry = indices(R(idx(i)), C(idx(i)));
        ix2 = cellEntry{1}(2);
        ix2 = ix2{1};
        dtwX = X(ix2);
        combination = combination + resample(dtwX, globalLength, length(dtwX));
        lengthWrap(C(idx(i))) = globalLength;
        figure;
        hold on;
        plot(abs(original));
        plot(abs(combination));
        xlim([1, 300]);
        pause
        close all
    elseif lengthWrap(C(idx(i))) > 0 
        X = r_new_downsampled_abs(:,R(idx(i)));
        cellEntry = indices(R(idx(i)), C(idx(i)));
        ix2 = cellEntry{1}(1);
        ix2 = ix2{1};
        dtwX = X(ix2);
        combination = combination + resample(dtwX, globalLength, length(dtwX));
        lengthWrap(R(idx(i))) = globalLength;
        figure;
        hold on;
        plot(abs(original));
        plot(abs(combination));
        xlim([1, 300]);
        pause
        close all
    else
        X1 = r_new_downsampled_abs(:,R(idx(i)));
        X2 = r_new_downsampled_abs(:,C(idx(i)));
        cellEntry = indices(R(idx(i)), C(idx(i)));
        ix1 = cellEntry{1}(1);
        ix2 = cellEntry{1}(2);
        ix1 = ix1{1};
        ix2 = ix2{1};
        dtwX1 = X1(ix1);
        dtwX2 = X2(ix2);
        dtwX = dtwX1+dtwX2;
        combination = combination + resample(dtwX, globalLength, length(dtwX));
        lengthWrap(R(idx(i))) = globalLength;
        lengthWrap(C(idx(i))) = globalLength;
        figure;
        hold on;
        plot(abs(original));
        plot(abs(combination));
        xlim([1, 300]);
        pause
        close all
    end
    i = i + 1;
    if i > length(idx)
        i = 1;
    end
end

%% Old code
% dtwMatrix = Inf*ones(size(r_new_matrix,2));
% indices = cell(size(r_new_matrix,2), size(r_new_matrix,2));
% 
% for i = 1:size(r_new_matrix,2)
%     for j = i+1:size(r_new_matrix,2)
% %         figure;
%         [dtwMatrix(i,j), ix, iy] = dtw(downsample(abs(r_new_matrix(:,i)),40), downsample(abs(r_new_matrix(:,j)),40));
%         indices{i,j} = {ix, iy};
% %         title(['i = ', num2str(i),' j = ',num2str(j)]);
% %         pause
% %         close all;
%     end
% end
% [R, C] = ndgrid(1:size(dtwMatrix,1), 1:size(dtwMatrix,2));
% [val, idx] = sort(dtwMatrix(:));
% 
% for i = 1:length(idx)
%     
% end


