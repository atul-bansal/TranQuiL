clear;

% The code with all the effects combined work till -18 dB of Signal to
% Noise Ratio as compared to -10 dB with only Cross Correlation. This
% provides about 8 dB of gain. 

%% Initial Parameters
numBits = 1000000;
rng(1,"twister");

%% Bit Stream
bitStream = randi([0 1],numBits,1);

%% Waveform Generation
cfgnonHT = wlanNonHTConfig('MCS', 3); 
waveform = wlanWaveformGenerator(bitStream, cfgnonHT);
nfft = 64;
cp = 16;
fs = 20e6;
ofdmSymbolLength = nfft + cp;
cbw = 'CBW16';
basePreamble = waveform(1:320);

%% Freq selective parameters
thresh = 1e-6;

%% Cyclostationarity Parameters
% cyclicResolution = 0.1;
% alpha = -1+cyclicResolution/nfft:cyclicResolution/nfft:(1-cyclicResolution/nfft);
% cycloWindowSize = ofdmSymbolLength;
% cycloOverlap = fix(0.5*cycloWindowSize);
load("CycloStationaryIndices.mat");


%% Overall Window Size
coarseWindow = length(waveform);
coarseOverlap = fix(0.9*coarseWindow);

%% Channel Model
sigpower = pow2db(mean(abs(waveform).^2));
tgahChan = wlanTGahChannel('SampleRate',fs,'ChannelBandwidth',cbw, ...
    'LargeScaleFadingEffect','Pathloss and shadowing', ...
    'DelayProfile','Model-D');
sigPwr_actual = 10^((sigpower-tgahChan.info.Pathloss)/10);
waveformChannel = tgahChan(waveform);

for snrValue = -10:-4:-25
    chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)',...
        'SNR',snrValue,'SignalPower', sigPwr_actual);
    signalReceived = [zeros(10000,1); waveformChannel; zeros(20000,1)];
    rxWaveform = chNoise(signalReceived); % Add AWGN noise
    K_coarse = fix((length(rxWaveform)-coarseOverlap)/(coarseWindow-coarseOverlap));
    cycloCorrArray = zeros(K_coarse, 1);
    pilotCorrArray = zeros(K_coarse, 1);
    indexCorrArray = zeros(K_coarse, 1);
    startIdxWindow = 1;
    for K = 1:K_coarse
        thisWaveform = rxWaveform(startIdxWindow:startIdxWindow+length(waveform)-1);
        decisionFreq = checkFreqSelectivity(thisWaveform, thresh);
        if decisionFreq
            filteredSignal = filterGuardBands(thisWaveform);
            cycloCorrArray(K) = checkCycloStationarity(filteredSignal,r,c);
            pilotCorrArray(K) = checkPilot(filteredSignal);
            indexCorrArray(K) = startIdxWindow;
        else
            cycloCorrArray(K) = 0;
            pilotCorrArray(K) = 0;
            indexCorrArray(K) = startIdxWindow;
        end
        startIdxWindow = startIdxWindow + coarseWindow - coarseOverlap;
    end

    % Baseline Calculation
    [rcorr, crossIdx] = xcorr(rxWaveform, basePreamble);
    crossCorrArray = rcorr(length(rxWaveform)+indexCorrArray);
    combinedCorrArray = 0.3*abs(crossCorrArray)/max(abs(crossCorrArray)) + 0.5*abs(cycloCorrArray)/max(abs(cycloCorrArray)) + 0.2*abs(pilotCorrArray)/max(abs(pilotCorrArray));

    subplot(4,1,1); plot(indexCorrArray, abs(cycloCorrArray)); title(["Cyclostationarity at SNR = ", num2str(snrValue)]);
    subplot(4,1,2); plot(indexCorrArray, abs(pilotCorrArray)); title(["Pilot Correlation at SNR = ", num2str(snrValue)]);
    % subplot(4,1,3); plot(indexCorrArray, abs(crossCorrArray)); title(["Sampled Cross Correlation at SNR = ", num2str(snrValue)]);
    subplot(4,1,4); plot(crossIdx, abs(rcorr)); title(["Full Cross Correlation at SNR = ", num2str(snrValue)]);
    subplot(4,1,3); plot(indexCorrArray, combinedCorrArray); title(["Combined Correlation at SNR = ", num2str(snrValue)]);
    close all
end

function cycloCorVal = checkCycloStationarity(rxwaveform,r,c)
    cyclicResolution = 0.1;
    nfft = 64;
    cp = 16;
    ofdmSymbolLength = nfft + cp;
    alpha = -1+cyclicResolution/nfft:cyclicResolution/nfft:(1-cyclicResolution/nfft);
    cycloWindowSize = ofdmSymbolLength;
    cycloOverlap = fix(0.5*cycloWindowSize);
    S = zeros(nfft,length(alpha));
    for alphaIdx = 1:length(alpha)
        S(:,alphaIdx) = calculateCoh(rxwaveform, alpha(alphaIdx), nfft, cycloOverlap,cycloWindowSize, cp);
    end
    S(:,ceil(length(alpha)/2)) = 0;
    cycloCorVal = 0;
    for idx = 1:length(r)/2
        cycloCorVal = cycloCorVal + S(r(idx), c(idx))*conj(S(r(length(r)-idx+1), c(length(c)-idx+1))); 
    end
end

function pilotCorVal = checkPilot(rxwaveform)
    pilotIndices = [-21, -7, 7, 21] + 33; % [1,1,1-1] for every OFDM symbol
    nfft = 64;
    cp = 16;
    ofdmSymbolLength = nfft + cp;
    curIdx = 4*ofdmSymbolLength+1;
    symbolNum = 0;
    pilotCorVal = 0;
    while curIdx < (length(rxwaveform)-ofdmSymbolLength)
        currSymbol = rxwaveform(curIdx:curIdx+ofdmSymbolLength-1);
        currSymbolfft = fftshift(fft(currSymbol(cp+1:end)));
        pilots = currSymbolfft(pilotIndices);
        actualPilots = wlan.internal.nonHTPilots(1,symbolNum);
        pilotCorVal = pilotCorVal + sum(pilots.*conj(actualPilots));
        curIdx = curIdx + ofdmSymbolLength;
        symbolNum = symbolNum + 1;
    end
end

function decision = checkFreqSelectivity(rxWaveform, thresh)
    % Getting Spectrogram
    [s, ~, ~] = spectro_OFDM(rxWaveform);

    %  ESD calculation
    esd_arr = mean(abs(s).^2, 2); % Mean of Abs value across each row and getting it wrt frequency
    esd_arr_sorted = sort(esd_arr);
    diff_esd_arr = diff(esd_arr_sorted);
    [pks, ~] = findpeaks(diff_esd_arr, 'NPeaks',5, 'SortStr','descend');
    if sqrt(var(pks)) > thresh
        decision = true;
    else
        decision = false;
    end
end

function filteredSignal = filterGuardBands(waveform)
    N = length(waveform);
    occupiedBw = 16.75e6;
    totalBw = 20e6;
    occupiedSamples = ceil(occupiedBw*length(waveform)/totalBw);
    emptySamplesOneSide = floor((length(waveform) - occupiedSamples)/2);
    filterFreqDomain = [zeros(emptySamplesOneSide, 1); ones(N - 2*emptySamplesOneSide, 1); zeros(emptySamplesOneSide,1)];
    filteredSignal = ifft(fftshift(fftshift(fft(waveform)).*filterFreqDomain));
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