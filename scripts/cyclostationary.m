clc;
clear;

rng(1,"twister");

%% Initial Parameters
M = 4;
numBits = 1000;
nfft = 64;
bw = 20e6;

%% Bit Stream
bitStream = randi([0 1],numBits,1);


%% 4-QAM Modulation
symbolStream = qammod(bitStream, M, InputType="bit",UnitAveragePower=true);
% scatterplot(symbolStream)
numZerosAppend = ceil(length(symbolStream)/nfft)*nfft - length(symbolStream);

% append zeros
symbolStream = [symbolStream; zeros(numZerosAppend, 1)];

%% OFDM Modulaiton
freqDomainSymbols = reshape(symbolStream, nfft, []);
ofdmSymbol = ofdmmod(freqDomainSymbols, nfft, 16);


%% Create a time domain waveform
txSignal = [zeros(800000, 1); repmat(ofdmSymbol, 100, 1); zeros(1200000,1)];

%% Add AWGN noise
% rxsignal = awgn(txSignal, -30, "measured");
rxsignal = txSignal;

%% Demodulation
% rxSymParallel = ofdmdemod(rxsignal, nfft, 16);
% rxSymSerial = reshape(rxSymParallel, [], 1);
% rxSymSerial = rxSymSerial(1:length(rxSymSerial)-numZerosAppend);
% scatterplot(rxSymSerial)

%% QAM demod
% rxbitStream = qamdemod(rxSymSerial, M, OutputType="bit",UnitAveragePower=true);
% ber = sum(abs(rxbitStream-bitStream))/length(bitStream);


%% Cyclostationarity
% ofdmSymbol = rxsignal;
% totalSamples = length(ofdmSymbol); % N
% ratioSamplesNfft = totalSamples/nfft;
% subcarrierSpacing = 1;
% cyclicspacing = 1/ratioSamplesNfft; % delta alpha
% 
% cyclicfRange = -nfft/4:cyclicspacing:nfft/4-cyclicspacing;
% freqRange = -nfft/2:nfft/2-1;
% 
% cyclicFunction = zeros(length(cyclicfRange), length(freqRange));
% 
% for i = 1:length(cyclicfRange)
%     cyclicFunction(i, :) = fftshift(calculateS(ofdmSymbol, cyclicfRange(i), nfft));
% end
% 
% 
% function S = calculateS(ofdmSymbol, cyclicf, nfft)
%     S = 0;
%     for i = 1:length(ofdmSymbol)-nfft+1
%         S = S + innerProduct(ofdmSymbol, i, cyclicf, nfft);
%     end
%     S = S/(length(ofdmSymbol)-nfft);
% end
% 
% function value = innerProduct(ofdmSymbol, i, cyclicf, nfft)
%     expValues = exp(-1i*2*pi*cyclicf*(i:i+nfft-1)/(2*nfft)).';
%     fft1 = fft(ofdmSymbol(i:i+nfft-1).*expValues);
%     fft2 = fft(ofdmSymbol(i:i+nfft-1).*conj(expValues));
%     value = fft1.*conj(fft2);
% end

%% Online Cyclostationarity
ofdmSymbol = rxsignal;
Nw = nfft;			% window length
Nv = fix(2/3*Nw);	% block overlap
n = length(ofdmSymbol);
da = 1/n;           % cyclic frequency resolution
a1 = -640;            % first cyclic freq. bin to scan (i.e. cyclic freq. a1*da)
a2 = 640;           % last cyclic freq. bin to scan (i.e. cyclic freq. a2*da)

S = zeros(nfft,a2-a1+1);

for k = a1:a2
    S(:, k+641) = calculateCoh(ofdmSymbol,k/n,nfft,Nv,Nw);
end
mesh(a1/n:1/n:a2/n, 1:64, abs(S))

function coh = calculateCoh(ofdmSymbol, alpha, nfft, Nv, Nw)
    Window = hanning(Nw);
    Window = Window(:);
    nwind = length(Window); % length of window
    n = length(ofdmSymbol);   % Number of data points
    K = fix((n-Nv)/(nwind-Nv));	% Number of windows

    index = 1:nwind;
    f = (0:nfft-1)/nfft;
    t = (0:n-1)';
    CPS = 0;

    x1 = ofdmSymbol.*exp(-1i*pi*alpha*t);
    x2 = ofdmSymbol.*exp(1i*pi*alpha*t);

    for i=1:K
        xw = Window.*x1(index);
        yw = Window.*x2(index);
        Yw1 = fftshift(fft(yw,nfft));		% Yw(f+a/2) or Yw(f)
        Xw2 = fftshift(fft(xw,nfft));		% Xw(f-a/2) or Xw(f-a)
        CPS = Yw1.*conj(Xw2) + CPS;
        index = index + (nwind - Nv);
    end
    % normalize
    KMU = K*norm(Window)^2;	% Normalizing scale factor ==> asymptotically unbiased
    CPS = CPS/KMU;
    coh = CPS;
end





