clear;
clc;

BW = 1.25e5;
SF = 10;
chirp_size=1024;
chirpsize=chirp_size;

extra_sampling_factor = 1;
Fs = BW;
symbol_length=chirp_size;
symbol_length_upsampled = extra_sampling_factor*chirp_size;
freq_shift_per_sample =  Fs/symbol_length; % How each frequency bin maps to a difference in frequency
Ts = 1/freq_shift_per_sample; % Symbol Duration
% fq = linspace(-BW/2,BW/2-freq_shift_per_sample,symbol_length); % The X-Axis
reset_freq = -BW/2; % The initial frequency of the base chirp
final_freq = (BW/2)-freq_shift_per_sample; % The final frequency
[up,down] = my_create_chirpspecial1(extra_sampling_factor*Fs,Ts,reset_freq,final_freq,chirp_size);
upfft = fft(up);
downfft = fft(down);

% BW = 250 kHz
% up250fftshift = [zeros(length(up)/2,1);fftshift(upfft);zeros(length(up)/2,1)];
% up250 = ifft(fftshift(up250fftshift));
% downfft = fft(down);
% down250fftshift = [zeros(length(down)/2,1);fftshift(downfft);zeros(length(down)/2,1)];
% down250 = ifft(fftshift(down250fftshift));
% corrSymbol = repmat(up250,2,1);
% packet = [zeros(10000,1);repmat(up250,10,1); zeros(10000,1)];

% BW = 1 MHz
% up1Mfftshift = [zeros(7*length(up)/2,1);fftshift(upfft);zeros(7*length(up)/2,1)];
% up1M = ifft(fftshift(up1Mfftshift));
% down1Mfftshift = [zeros(7*length(down)/2,1);fftshift(downfft);zeros(7*length(down)/2,1)];
% down1M = ifft(fftshift(down1Mfftshift));
% corrSymbol = repmat(up1M,2,1);
% packet = [zeros(10000,1);repmat(up1M,10,1); zeros(10000,1)];

% BW = 6 MHz
up6Mfftshift = [zeros(47*length(up)/2,1);fftshift(upfft);zeros(47*length(up)/2,1)];
up6M = ifft(fftshift(up6Mfftshift));
down6Mfftshift = [zeros(47*length(down)/2,1);fftshift(downfft);zeros(47*length(down)/2,1)];
down6M = ifft(fftshift(down6Mfftshift));

packet = [zeros(10000,1);repmat(up6M,10,1); zeros(10000,1)];
corrSymbol = repmat(up6M,2,1);
% OFDM signal

message = randi([0 1], 10240,1);
msg_symbols = qammod(message,4,'InputType','bit');
num_subcarriers = 64;
blocks = length(msg_symbols)/num_subcarriers;
msg = reshape(msg_symbols,num_subcarriers,[]);
ofdm_msg = zeros(num_subcarriers,blocks);
for i = 1:blocks
    ofdm_msg(:,i) = ifft(msg(:,i));
end
time_ofdm_msg = reshape(ofdm_msg, length(msg_symbols),1);
time_ofdm_msg = [zeros(10000,1); time_ofdm_msg; zeros(10000,1)];
pow_ofdm = rms(time_ofdm_msg)^2;

% Noise signal
noise = sqrt(pow_ofdm)*randn(length(time_ofdm_msg),1);


% Plots
[r1,lag1] = xcorr(packet, corrSymbol);
[r2, lag2] = xcorr(time_ofdm_msg, corrSymbol);
[r3, lag3] = xcorr(noise, corrSymbol);

figure;
hold on;
% plot(lag1, abs(r1));

plot(lag3, abs(r3));
plot(lag2, abs(r2));
legend('Noise','OFDM');
title('6 MHz OFDM');