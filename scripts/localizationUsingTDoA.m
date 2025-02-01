clc;
clear;
tic
%% Load TDoA values and Ground Truth Coordinates
load("Results_TDoA_array.mat");
load("Gnd_truth_coordinates_cartesian_greenbank.mat");

% GPS coordinates would be like - 5 Base Station locations in the below order and then 15 Transmitter locations - so in total
% 20 locations
baseStations = {'school','church','northfork','buffalo','hannah'};
tdoaBaseStationOrder = [1 2; 1 3; 1 4; 3 4; 5 4];

%% Preprocessing the TDoA Results
measured_tdoa_values = zeros(30,5);
for i = 1:30
    measured_tdoa_values(i,1) = Results_TDoA(i,2);
    measured_tdoa_values(i,2) = Results_TDoA(i,5);
    measured_tdoa_values(i,3) = Results_TDoA(i,8);
    measured_tdoa_values(i,4) = Results_TDoA(i,11);
    measured_tdoa_values(i,5) = Results_TDoA(i,14);
end

%% Localization Starts
x = (floor(min(xEast))-10):10:(ceil(max(xEast))+10);
y = (ceil(max(yNorth))+10):-10:(floor(min(yNorth))-10);
[X,Y] = meshgrid(x,y);


strip_length = 5;
error_values = zeros(30,1);
for idx = 1:30
    intensity = zeros(length(y), length(x));
    for k = 1:size(measured_tdoa_values,2)
        if isnan(measured_tdoa_values(idx,k))
            continue;
        end
        for i=1:length(x)
            i
            for j=1:length(y)
                if (( (pdist([X(j,i),Y(j,i);xEast(tdoaBaseStationOrder(k,1)),yNorth(tdoaBaseStationOrder(k,1))],'euclidean')-pdist([X(j,i),Y(j,i);xEast(tdoaBaseStationOrder(k,2)),yNorth(tdoaBaseStationOrder(k,2))],'euclidean')) <= (measured_tdoa_values(idx,k) + 5)) ...
                        && ( (pdist([X(j,i),Y(j,i);xEast(tdoaBaseStationOrder(k,1)),yNorth(tdoaBaseStationOrder(k,1))],'euclidean')-pdist([X(j,i),Y(j,i);xEast(tdoaBaseStationOrder(k,2)),yNorth(tdoaBaseStationOrder(k,2))],'euclidean')) >= (measured_tdoa_values(idx,k) - 5)))
                    intensity(j,i) = intensity(j,i) + 1;
                end
            end
        end
    end
    maximum = max(max(intensity));
    [m,n] = find(intensity == maximum);
    min_error = Inf;
    for t=1:length(m)
        candidate_location = [X(m(t),n(t)) Y(m(t),n(t))]; 
        candidate_error = pdist([candidate_location;xEast(ceil(idx/2)+5),yNorth(ceil(idx/2)+5)],'euclidean');
        if candidate_error < min_error
            min_error = candidate_error;
            calculated_location = candidate_location;
        end
    end
    error_values(idx) = min_error;
    hmh = heatmap(intensity, 'GridVisible','off','FontSize',20);
    hmh.ColorLimits = [0, max(max(intensity))];
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    close all;
end
toc