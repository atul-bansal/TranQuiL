% GPS coordinates would be like - 5 Base Station locations and then 15 Transmitter locations - so in total
% 20 locations

clc;
clear;

GPSCoordinates = zeros(20,2);
baseStations = {'school','church','northfork','buffalo','hannah'};

for i = 1:5
    gpsfilename = ['./GPS_GND_TRUTH/gps_',baseStations{i},'.dat'];
    [lat_degree, long_degree] = gps_read(gpsfilename);
    GPSCoordinates(i,:) = [lat_degree, long_degree];
end

for i = 6:20
    tx1 = ['tx', num2str(i-5)]; % Remember to use correct GPS filenames for 12 and 13 - uncomment line 30
    if (i-5) == 12 || (i-5) == 13
        tx1 = [tx1, '_rep'];
    end
    gpsfilename = ['./GPS_GND_TRUTH/gps_', tx1,'.dat'];
    [lat_degree, long_degree] = gps_read(gpsfilename);
    GPSCoordinates(i,:) = [lat_degree, long_degree];
end

alt = 234;
origin = [GPSCoordinates(1,:) alt]; % Assigning Base Station 1 as Origin
[xEast,yNorth] = latlon2local(GPSCoordinates(:,1),GPSCoordinates(:,2),alt,origin);

function [lat_degree, long_degree] = gps_read(path)
%% Only works in North America

rawdata = readcell(path,'Delimiter',',');
a = cell2mat(rawdata(:,3)')';
lat_array = floor(a/100)+mod(a,100)/60;
lat_array_rad = lat_array * pi/180;
lat_array_exp = exp(1i*lat_array_rad);
mid = mean(lat_array_exp);
lat_degree = atan(imag(mid)/real(mid))*180/pi;

b = cell2mat(rawdata(:,5)')';
long_array = -1*(floor(b/100)+mod(b,100)/60);
long_array_rad = long_array * pi/180;
long_array_exp = exp(1i*long_array_rad);
mid = mean(long_array_exp);
long_degree = atan(imag(mid)/real(mid))*180/pi;

end