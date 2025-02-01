clear;
clc;
format long
%% Transmitter
message = 'Hello! My name is Atul Bansal. I am a 3rd year PhD student and currently working under Professor Swarun Kumar.';
nfft = 64;
preamble = [1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1,1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1];
preamble_symbols = qammod(preamble.',4,'InputType','bit');
len = length(message);
bits = reshape(dec2bin(double(message).', 8).', 1, [])-'0';
bits = [bits, zeros(1, mod(-length(bits), 2*nfft))];
nsym = length(bits)/(2*nfft);
output = [];

output = bits;
output = qammod(output.',4,'InputType','bit');
output = [preamble_symbols; output];

nsym = length(output)/nfft;
for ii = 1:nsym
% Collect the iith symbol
symbol = output((ii-1)*nfft+1:ii*nfft);
% Run an IFFT on the symbol
output((ii-1)*nfft+1:ii*nfft) = ifft(symbol);
end
output = [output; zeros(length(output),1)];
% write_complex_binary(output, 'wifi_tx_ofdm.dat');

%% Receiver
% a = read_complex_binary('ofdm_movement_AoAch_0_binary_2');
% b = read_complex_binary('ofdm_movement_AoAch_1_binary_2');
a = read_complex_binary('ofdm_movement_AoA_1_1.dat');
b = read_complex_binary('ofdm_movement_AoA_1_2.dat');
a = a(1e5+1:end);
b = b(1e5+1:end);


%%
len = length(a);
relchannelMatrix = [];
relPhase = [];

for cnt = 1:24
aoffsetcorrect = a((cnt-1)*2e6+1:cnt*2e6);
boffsetcorrect = b((cnt-1)*2e6+1:cnt*2e6);
% aoffsetcorrect = a(1:2e6);
% boffsetcorrect = b(1:2e6);

[r1, idx1] = xcorr(aoffsetcorrect, ifft(preamble_symbols));
[r2, idx2] = xcorr(boffsetcorrect, ifft(preamble_symbols));

[~,a1]=findpeaks(abs(r1),'MinPeakHeight',prctile(abs(r1),85),'MinPeakDistance',1e3);
a1 = a1 - length(aoffsetcorrect);
a1(a1<0) = [];
for i = 1:100:length(a1)
    for j = 4
        if(a1(i)+j*64 < length(aoffsetcorrect))
            relchannelMatrix = [relchannelMatrix angle(fftshift(fft(aoffsetcorrect(a1(i)+(j-1)*64+1:a1(i)+j*64))./fft(boffsetcorrect(a1(i)+(j-1)*64+1:a1(i)+j*64))))];
        end
    end
end    
end

k = length(relchannelMatrix(32,:));
PHASE = zeros(nfft,1);
for j = 1:nfft
p1 = polyfit(1:55, unwrap(relchannelMatrix(j,1:55)),1);
p2 = polyfit(k-55+1:k, unwrap(relchannelMatrix(j,k-55+1:k)),1);

slopediff = p2(1) - p1(1);
slopediff

phaseJump = p2(2) - p1(2);
phaseJump = mod(phaseJump, 2*pi);
PHASE(j) = phaseJump;
end  

% frequencies = linspace(2461e6, 2463e6, nfft);
% h = exp(1i*PHASE);
% frequencies = frequencies.';
% dBS = 0.16;
% BS_station_distance =  ceil(dBS) + 1;
% BS_station_distance_ns = BS_station_distance/0.3;
% alpha_start = 100;
% alpha_end = 3000;
% 
% [optimal_p, tau_in_ns, val, locs, optimal_alpha] = inverseNDFT_robust_alphaReturn(h, frequencies, BS_station_distance_ns,alpha_start, alpha_end);

% I think Bartlett makes more sense here as ToF will give me 0 ToF. Because
% of so close base stations, cant differentiate between 2 paths separated
% by less than c/BW = 15 m. To be continued. 
