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

beaconPreamble = txWaveform(1:944); % 944 number by - 400 Physical PReamble samples + Some MAC header samples
periodicbeaconPreamble = [beaconPreamble; zeros(length(txWaveform)-length(beaconPreamble),1)]; 
periodicbeaconPreamble15 = repmat([periodicbeaconPreamble; zeros(20e6*102.4e-3+80,1)],15,1);

%% Timestamp Variation
timestamp = randi(2^53-1);
txWaveform = [];

for timeValue = 1:15
frameBodyConfig = wlanMACManagementConfig( ...
    BeaconInterval=beaconInterval, ...
    SSID=ssid,...
    Timestamp = timestamp+timeValue);

dsElementID = 3;
dsInformation = dec2hex(chNum,2);
frameBodyConfig = frameBodyConfig.addIE(dsElementID,dsInformation);

beaconFrameConfig = wlanMACFrameConfig(FrameType="Beacon", ...
    ManagementConfig=frameBodyConfig);
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
figure
plot(abs(xcorr(txWaveform, beaconPreamble)));
title('TimeStamp Variation');

%% TimeStamp + Sequence Control field variation
timestamp = randi(2^53-1);
seqValue = randi(4095);
txWaveform = [];

for timeValue = 1:15
frameBodyConfig = wlanMACManagementConfig( ...
    BeaconInterval=beaconInterval, ...
    SSID=ssid,...
    Timestamp = timestamp+timeValue);

dsElementID = 3;
dsInformation = dec2hex(chNum,2);
frameBodyConfig = frameBodyConfig.addIE(dsElementID,dsInformation);

beaconFrameConfig = wlanMACFrameConfig(FrameType="Beacon", ...
    ManagementConfig=frameBodyConfig, SequenceNumber=seqValue+timeValue);
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
figure
plot(abs(xcorr(txWaveform, beaconPreamble)));
title('Timestamp + Sequence Number variation');

%% TimeStamp + Sequence Control + Source MAC address Variation
timestamp = randi(2^53-1);
seqValue = randi(4095);
sourceMAC = '60B76E4363FC';
txWaveform = [];
deltaf = 125e3;
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
txWaveform = txWaveform .* exp(2i*pi*deltaf*(1:length(txWaveform)).'*(1/20e6))* exp(1i*2*pi*rand);
rxWaveform = awgn(txWaveform, snrdB, 'measured');

figure
plot(abs(xcorr(rxWaveform, beaconPreamble)));
title('Timestamp + Sequence Number + Source MAC variation');

%% Introducing Periodicity - got 13 dB gain in Simulation
timestamp = randi(2^53-1);
seqValue = randi(4095);
sourceMAC = '60B76E4363FC';
txWaveform = [];
snrdB = -56;

for timeValue = 1:30
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
rxWaveform = awgn(txWaveform, snrdB, 'measured');

figure;
plot(abs(xcorr(rxWaveform, periodicbeaconPreamble15)));
% plot(abs(xcorr(txWaveform, periodicbeaconPreamble15)));
title('Timestamp + Sequence Number + Source MAC variation +  Periodicity');

%% Introducing Frequency Offsets
timestamp = randi(2^53-1);
seqValue = randi(4095);
sourceMAC = '60B76E4363FC';
txWaveform = [];
snrdB = -25;
deltaf = 125e3;

for timeValue = 1:30
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
t = ((1:length(txWaveform)).*(1/20e6)).';
txWaveform = txWaveform .* exp(-2i*pi*deltaf*t) * exp(1i*2*pi*rand);
rxWaveform = awgn(txWaveform, snrdB, 'measured');

figure;
plot(abs(xcorr(rxWaveform, periodicbeaconPreamble15)));
% plot(abs(xcorr(txWaveform, periodicbeaconPreamble15)));
title('Timestamp + Sequence Number + Source MAC variation +  Periodicity + Frequency Offsets');