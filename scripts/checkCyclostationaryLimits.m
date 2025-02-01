clear;
% Starting to go bad at SNR = -7 dB

%% Initial Parameters
numBits = 1000000;


%% Bit Stream
bitStream = randi([0 1],numBits,1);

%% Waveform Generation
cfgnonHT = wlanNonHTConfig('MCS', 3); 
waveform = wlanWaveformGenerator(bitStream, cfgnonHT);
nfft = 64;
cp = 16;
ofdmSymbolLength = nfft + cp;

%% Limits of Noise
sigpower = pow2db(mean(abs(waveform).^2));
txSignal = [zeros(10000,1); waveform; zeros(10000,1)];
snrRange = 50:-1:0;
cyclicResolution = 0.1;
alpha = -1+cyclicResolution/nfft:cyclicResolution/nfft:(1-cyclicResolution/nfft);

windowsize = ofdmSymbolLength;
overlap = fix(0.5*windowsize);
% overlap = 0;


%% Extra thing to be deleted later
txSignal = waveform;
for idx = 1:length(snrRange)
    rxsignal = awgn(txSignal, snrRange(idx), sigpower);
    curwindow = txSignal;
    S = zeros(nfft,length(alpha));
    for alphaIdx = 1:length(alpha)
        S(:,alphaIdx) = calculateCoh(curwindow, alpha(alphaIdx), nfft, overlap, windowsize, cp);
    end
    S(:,ceil(length(alpha)/2)) = 0;
    mesh(alpha*nfft, -31:32, abs(S))
    [pks,locs] = findpeaks(abs(S(:)), 'MinPeakDistance', 100, 'MinPeakHeight', 35, 'NPeaks', 20);
    [r,c] = ind2sub(size(S), locs);
    checkPhase = S(locs); % Phases are conjugate of each other
end



%% Actual Code

for idx = 1:length(snrRange)
    rxsignal = awgn(txSignal, snrRange(idx), sigpower);
    cyclicCorrArray = zeros(length(rxsignal)-length(waveform),1);
    for curIdx = 1:length(rxsignal)-length(waveform)
        curwindow = rxsignal(curIdx:curIdx+length(waveform)-1);
        S = zeros(nfft,length(alpha));
        for alphaIdx = 1:length(alpha)
            S(:,alphaIdx) = calculateCoh(curwindow, alpha(alphaIdx), nfft, overlap, windowsize, cp);
        end
        S(:,ceil(length(alpha)/2)) = 0;
        mesh(alpha*nfft, -31:32, abs(S))
        cyclicCorrArray(curIdx) = max(abs(S(:)));
    end
    plot(abs(cyclicCorrArray))
end


function coh = calculateCoh(waveform, alpha, nfft, overlap, windowsize, cp)
    n = length(waveform);   % Number of data points
    K = fix((n-overlap)/(windowsize-overlap));	% Number of windows

    index = 1:windowsize;
    % f = (0:nfft-1)/nfft;
    t = (0:n-1)';
    CPS = 0;

    x1 = waveform.*exp(-1i*pi*alpha*t);
    x2 = waveform.*exp(1i*pi*alpha*t);

    for i=1:K
        xw = x1(index);
        yw = x2(index);
        xw = xw(cp+1:end);
        yw = yw(cp+1:end);
        Yw1 = fftshift(fft(yw,nfft));		% Yw(f+a/2) or Yw(f)
        Xw2 = fftshift(fft(xw,nfft));		% Xw(f-a/2) or Xw(f-a)
        CPS = Yw1.*conj(Xw2) + CPS;
        index = index + (windowsize - overlap);
    end
    % normalize
    KMU = K;	% Normalizing scale factor ==> asymptotically unbiased
    CPS = CPS/KMU;
    coh = CPS;
end
