clear;
clc;

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
a = read_complex_binary('ofdm_movement_AoAch_0_binary_2');
b = read_complex_binary('ofdm_movement_AoAch_1_binary_2');
% a = read_complex_binary('ofdm_movement_AoA_1_1.dat');
% b = read_complex_binary('ofdm_movement_AoA_1_2.dat');
a = a(1e5+1:end);
b = b(1e5+1:end);
%%
len = length(a);
relPhase = [];
figure;
hold on;

for cnt = 1:37
acnt = a((cnt-1)*2e6+1:cnt*2e6);
bcnt = b((cnt-1)*2e6+1:cnt*2e6);
[r1, idx1] = xcorr(acnt, ifft(preamble_symbols));
[r2, idx2] = xcorr(bcnt, ifft(preamble_symbols));

[~,a1]=findpeaks(abs(r1),'MinPeakHeight',prctile(abs(r1),85),'MinPeakDistance',1e3);
[~,a2]=findpeaks(abs(r2),'MinPeakHeight',prctile(abs(r2),85),'MinPeakDistance',1e3);
a1 = a1 - length(acnt);
a1(a1<0) = [];
for i = 1:500:length(a1)
    for j = 4
        if(a1(i)+j*64 < length(acnt))
            relPhase = [relPhase; angle(fftshift(fft(acnt(a1(i)+(j-1)*64+1:a1(i)+j*64))./fft(bcnt(a1(i)+(j-1)*64+1:a1(i)+j*64))))];
        end
    end
    plot(relPhase,'bx');
    pause(3);
end
end    

% a5 = a(1.337e6+68:1.337e6+1000);
% plot(real(a5))
% b5 = b(1.337e6+68:1.337e6+1000);
% plot(real(b5))
% plot(angle(fft(a5(1:64))./fft(b5(1:64))), 'bx')
% plot(angle(fft(a5(64+1:64+64))./fft(b5(64+1:64+64))), 'bx')
% plot(angle(fft(a5(2*64+1:2*64+64))./fft(b5(2*64+1:2*64+64))), 'bx')
% plot(angle(fft(a5(8*64+1:8*64+64))./fft(b5(8*64+1:8*64+64))), 'bx')
% 
% 
% a6 = a(4.1945e6+82:4.1945e6+1000+82);
% plot(real(a6))
% b6 = b(4.1945e6+82:4.1945e6+1000+82);
% plot(real(b6))
% plot(angle(fft(a6(64+1:64+64))./fft(b6(64+1:64+64))), 'bx'); ylim([-pi pi])
% plot(angle(fft(a6(2*64+1:2*64+64))./fft(b6(2*64+1:2*64+64))), 'bx'); ylim([-pi pi])
% plot(angle(fft(a6(8*64+1:8*64+64))./fft(b6(8*64+1:8*64+64))), 'bx'); ylim([-pi pi])
% xlabel('Number of Subcarriers');
% ylabel('h1/h2 in radians');
% 
