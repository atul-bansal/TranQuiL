clc;
clear;

roadLengths = [23 11 9 4 5 2 10];
defaultRoad = 18;

distancesPoints = zeros(15, 7);
offsetsIdxDist = zeros(15,2);

%% Baseline Distance Calculation
offsetsIdxDist(1,:) = [1 1];
offsetsIdxDist(2,:) = [1 2];
offsetsIdxDist(3,:) = [7 0.5];
offsetsIdxDist(4,:) = [4 3.7];
offsetsIdxDist(5,:) = [2 10.4];
offsetsIdxDist(6,:) = [1 2.6];
offsetsIdxDist(7,:) = [1 3.5];
offsetsIdxDist(8,:) = [1 2.5];
offsetsIdxDist(9,:) = [3 1.6];
offsetsIdxDist(10,:) = [3 0.6];
offsetsIdxDist(11,:) = [2 3];
offsetsIdxDist(12,:) = [2 1.5];
offsetsIdxDist(13,:) = [6 0];
offsetsIdxDist(14,:) = [NaN NaN];
offsetsIdxDist(15,:) = [4 1.5];
defaultOffset = 17;

distances = [1.8 1.36 2.86 2.2 0.5 1 0.85 1.12 1.66 2.1 1.88 2.82 1.78 0.5 2.83];

for point = 1:15
    if ~isnan(offsetsIdxDist(point,1))
        idxwithout = setdiff(1:7, offsetsIdxDist(point,1));
        distancesPoints(point, 1) = defaultRoad+roadLengths(offsetsIdxDist(point,1))-offsetsIdxDist(point,2);
        for dist = 1:6
            distancesPoints(point, dist+1) = distancesPoints(point, dist)+roadLengths(idxwithout(dist)); 
        end
    else
        distancesPoints(point,:) = defaultRoad - defaultOffset;
    end
end

meanDistancesPoints = mean(distancesPoints,2);
ts = [];
CI = [];
for i = 1:15
    sdDistancesPoints(i) = std(distancesPoints(i,:))/sqrt(length(distancesPoints(i,:)));
    ts_temp = tinv([0.025 0.975], length(distancesPoints(i,:))-1);
    ts = [ts; ts_temp];
    CI = [CI; mean(distancesPoints(i,:))+ts_temp*sdDistancesPoints(i)];
end

%% Improved Detection Calculation - Copy results from Results_Matrix_Greenbank_experiments
distancesPointsSamples = zeros(15,14);
numSamplePoints = zeros(15,1);
numSamplePoints(1) = 1; distancesPointsSamples(1,1) = 17.1;
numSamplePoints(2) = 1; distancesPointsSamples(2,1) = 16.6;
numSamplePoints(3) = 14; distancesPointsSamples(3,1:numSamplePoints(3)) = [1.74+18 0.4+18 5+0.4+18 2+0.4+18 2+1.74+18 5+1.74+18 9+1.74+18 9+0.4+18 11+0.4+18 11+1.74+18 23+1.74+18 23+0.4+18 2+5+0.4+18 2+5+1.74+18];
numSamplePoints(4) = 12; distancesPointsSamples(4,1:numSamplePoints(4)) = [0.5+18 0.5+2+18 0.5+5+18 0.5+9+18 0.5+11+18 0.5+23+18 0.5+2+5+18 0.5+2+9+18 0.5+2+11+18 0.5+5+9+18 0.5+2+23+18 0.5+2+5+9+18];
numSamplePoints(5) = 1; distancesPointsSamples(5,1:numSamplePoints(5)) = 0.5;
numSamplePoints(6) = 1; distancesPointsSamples(6, 1:numSamplePoints(6)) = 1.26;
numSamplePoints(7) = 2; distancesPointsSamples(7, 1:numSamplePoints(7)) = [0.3 18+0.3];
numSamplePoints(8) = 14; distancesPointsSamples(8,1:numSamplePoints(8)) = [18+1.13 18+11+1.13 18+2+1.13 18+4+1.13 18+5+1.13 18+2+1+1.13 18+4+1+1.13 18+5+1.13 18+2+4+1.13 18+2+5+1.13 18+4+5+1.13 18+2+4+1+1.13 18+2+5+1+1.13 18+4+5+1+1.13];
numSamplePoints(9) = 12; distancesPointsSamples(9,1:numSamplePoints(9)) = [18+11+0.3 18+11+23+0.1 18+2+11+0.3 18+11+4+0.3 18+11+5+0.3 18+11+2+4+0.3 18+11+2+5+0.3 18+4+5+11+0.3 18+23+0.1 18+2+23+0.1 18+4+23+0.1 18+23+5+0.1];
numSamplePoints(10) = 12; distancesPointsSamples(10, 1:numSamplePoints(10)) = [18+8.7 18+21.2 18+11+8.7 18+11+21.2 18+4+8.7 18+4+21.2 18+5+8.7 18+5+21.2 18+2+21.2 18+2+8.7 18+10+21.2 18+10+8.7];
numSamplePoints(11) = 12; distancesPointsSamples(11, 1:numSamplePoints(11)) = [18+6 18+9.5 18+23+6 18+23+9.5 18+4+6 18+4+9.5 18+5+6 18+5+9.5 18+2+6 18+2+9.5 18+10+6 18+10+9.5];
numSamplePoints(12) = 7; distancesPointsSamples(12, 1:numSamplePoints(12)) = [18+7 18+23+7 18+9+7 18+4+7 18+2+7 18+10+7 18+2+4+7];
numSamplePoints(13) = 1; distancesPointsSamples(13, 1:numSamplePoints(13)) = 16;
numSamplePoints(14) = 1; distancesPointsSamples(14, 1:numSamplePoints(14)) = 0;
numSamplePoints(15) = 10; distancesPointsSamples(15, 1:numSamplePoints(15)) = [18+0.3 18+2+0.3 18+5+0.3 18+10+0.3 18+9+0.3 18+11+0.3 18+23+0.3 18+2+5+0.3 18+2+9+0.3 18+2+10+0.3];


for i = 1:15
    meanDistancesPoints_Improved(i) = mean(distancesPointsSamples(i,1:numSamplePoints(i)));
end
meanDistancesPoints_Improved = transpose(meanDistancesPoints_Improved);

ts_improved = [];
CI_improved = [];
for i = 1:15
    sdDistancesPoints_improved(i) = std(distancesPointsSamples(i,1:numSamplePoints(i)))/sqrt(numSamplePoints(i));
    ts_temp_improved = tinv([0.025 0.975], numSamplePoints(i)-1);
    ts_improved = [ts_improved; ts_temp_improved];
    CI_improved = [CI_improved; mean(distancesPointsSamples(i,1:numSamplePoints(i)))+ts_temp*sdDistancesPoints_improved(i)];
end



%% Calculate Mean and SD of both baseline and improved
[sortedDist, sortIdx] = sort(distances);
meanDistancesPointsSorted = meanDistancesPoints(sortIdx);
sdDistancesPointsSorted = sdDistancesPoints(sortIdx);
CI_sorted = CI(sortIdx,:);

meanDistancesPoints_Improved_sorted = meanDistancesPoints_Improved(sortIdx);
sdDistancesPoints_Improved_sorted = sdDistancesPoints_improved(sortIdx);
CI_improved_sorted = CI_improved(sortIdx,:);
convertFactor = 60/25;

%% Plot Both baseline and Imrpoved
% figure; hold on;
% errorbar(sortedDist*1.6, convertFactor*meanDistancesPointsSorted, convertFactor*(meanDistancesPointsSorted-max(CI_sorted(:,1),0)), convertFactor*(CI_sorted(:,2)-meanDistancesPointsSorted));
% errorbar(sortedDist*1.6, convertFactor*meanDistancesPoints_Improved_sorted, convertFactor*(meanDistancesPoints_Improved_sorted-max(CI_improved_sorted(:,1),0)), convertFactor*(CI_improved_sorted(:,2)-meanDistancesPoints_Improved_sorted));


%% Creating Bar Plots
[Y, E] = discretize(sortedDist,7); % E gives the edges, Y gives the indices where element lies
binEdges = E;
binCounts_original = zeros(length(E)-1,1);
binCount_Sd_original = zeros(length(E)-1,1);
binCounts_improved = zeros(length(E)-1,1);
binCount_Sd_improved = zeros(length(E)-1,1);
for binNo = 1:length(binCounts_original)
    valuesBin_original = meanDistancesPointsSorted(Y == binNo);
    valuesBin_improved = meanDistancesPoints_Improved_sorted(Y == binNo);
    valuesBinSD_original = sdDistancesPointsSorted(Y == binNo);
    valuesBinSD_improved = sdDistancesPoints_Improved_sorted(Y == binNo);
    binCounts_original(binNo) = mean(valuesBin_original);
    binCounts_improved(binNo) = mean(valuesBin_improved);
    binCount_Sd_original(binNo) = mean(valuesBinSD_original);
    binCount_Sd_improved(binNo) = mean(valuesBinSD_improved);
end

errorbarvariationPts = E(1:end-1) + (E(2) - E(1))/2;
histogram('BinEdges',E*1.6,'BinCounts',binCounts_original*convertFactor);
hold on;
histogram('BinEdges',E*1.6,'BinCounts',binCounts_improved*convertFactor);
er1 = errorbar(errorbarvariationPts*1.6, binCounts_original*convertFactor, binCount_Sd_original*convertFactor);
er2 = errorbar(errorbarvariationPts*1.6, binCounts_improved*convertFactor, binCount_Sd_improved*convertFactor);
er1.Color = [0 0 1];
er2.Color = [1 0 0];
er1.LineStyle = 'none';
er2.LineStyle = 'none';
hold off;
