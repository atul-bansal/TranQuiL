clc;
clear;
rng('shuffle');


%% Intialization Parameters
num_tx = 3;  % Number of transmitters in the space
num_sensors = 10;   % Number of spectrum sensors 
tx_power = 20; % Transmission Power in dBm
noise_floor = -80:5:-30; % Noise floor in dBm
freq_range = 501:1:1000; % Frequency Range

x = 1:1000;
y = 1:1000;
[X,Y] = meshgrid(x,y);
total_isDetected = zeros(length(noise_floor),1);

for noise = 1:length(noise_floor)
noise
for trials = 1:20
trials
%% Creating Ground truth
x_rand = randperm(length(x), num_tx+num_sensors);
y_rand = randperm(length(y), num_tx+num_sensors);
x_rand = x_rand.';
y_rand = y_rand.';
sensors_coord = [x_rand(1:num_sensors) y_rand(1:num_sensors)];
tx_coord = [x_rand(num_sensors+1: num_sensors+num_tx) y_rand(num_sensors+1: num_sensors+num_tx)];
tx_freq = zeros(num_tx,1);
for i = 1:num_tx
   tx_freq(i) = freq_range(randi(length(freq_range)));
end

gnd_truth_power = ones(length(y), length(x), length(freq_range));
gnd_truth_power = noise_floor(noise)*gnd_truth_power - 1000;
for tx = 1:num_tx
    r = find(freq_range == tx_freq(tx));
    for i = 1:length(y)
        for j = 1:length(x)
            d = sqrt((X(i,j) - tx_coord(tx,1))^2 + (Y(i,j) - tx_coord(tx,2))^2);
            if d == 0
                gnd_truth_power(i,j,r) = tx_power;
            else
                gnd_truth_power(i,j,r) = tx_power - 20*log10(d) - 20*log10(tx_freq(tx)*1e6) - 20*log10(4*pi/3e8);
            end
        end
    end
end


%% Actual Sensing Task
center_freq = 510:20:990;
sensed_power = noise_floor(noise)*ones(length(center_freq), num_sensors) - 1000;
isdetect = zeros(length(center_freq),1);
for iter = 1:length(center_freq)
    for num = 1:num_sensors
        loc = sensors_coord(num,:);
        sensed_spectrum = gnd_truth_power(loc(2), loc(1), (iter-1)*20 + 1:iter*20);
        if isempty(sensed_spectrum(sensed_spectrum > noise_floor(noise))) ~= 1
            r = find(sensed_spectrum > noise_floor(noise));
            sensed_power(iter,num) = mean(sensed_spectrum(r));
            isdetect(iter) = 1;
            break;
        elseif isempty(sensed_spectrum(sensed_spectrum > noise_floor(noise) - 10)) ~= 1 
            r = find(sensed_spectrum > noise_floor(noise)-10);
            sensed_power(iter,num) = mean(sensed_spectrum(r));
            isdetect(iter) = 0.5;
        else
            sensed_power(iter,num) = max(sensed_spectrum);
        end
    end
    if isdetect(iter) == 0.5
        if 10*log10(sum(sqrt(10.^(0.1*sensed_power(iter,:))))^2) > noise_floor(noise)
            isdetect(iter) = 1;
        else
            isdetect(iter) = 0;
        end
    end
end

detected_tx = sum(isdetect);
total_isDetected(noise) = (total_isDetected(noise)*(trials-1) + detected_tx)/trials;
end
end